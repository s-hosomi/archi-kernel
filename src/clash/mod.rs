//! Interference (clash) checking — hard clashes and sleeve verification
//! (`DESIGN.md` §6-6, §10 Phase 7; `docs/research/06-domain.md` §6).
//!
//! The building model is a collection of placed members; clash checking finds
//! the pairs that physically overlap (a hard clash) and verifies that beam
//! penetrations (sleeves) obey the structural diameter / position rules. Both
//! run on the members' *evaluated* B-reps, so openings and priority clips are
//! already resolved.
//!
//! # Two-phase hard-clash detection
//!
//! 1. **Coarse phase** — every pair's axis-aligned bounding boxes are compared
//!    (`O(n²)`; see [`clash_check`] for why that is acceptable at building member
//!    counts). Only overlapping boxes proceed.
//! 2. **Fine phase** — the exact intersection *volume* of the pair is computed
//!    with the prismatic boolean ([`prismatic::intersection`]). A volume above a
//!    tolerance-derived threshold is a [`HardClash`](ClashKind::HardClash);
//!    a box overlap that yields (essentially) zero shared volume is a
//!    [`Touching`](ClashKind::Touching) contact.
//!
//! # Honest degeneracy — no silent "no clash"
//!
//! The exact volume is only available when the pair shares a common prismatic
//! direction (`DESIGN.md` §4.2). Two H-sections crossed at a right angle, or a
//! round column meeting a beam off-axis, have **no** common direction, so the
//! kernel cannot (in this phase) compute their shared volume. Rather than
//! silently report "no clash" — the exact failure mode the design forbids
//! (`DESIGN.md` §6-4, §6-6) — such a pair whose boxes overlap is returned as a
//! [`PotentialClash`](ClashKind::PotentialClash): an explicit, machine-readable
//! "the boxes overlap but the exact volume is undecidable here", carrying *no*
//! volume. The caller (or a later 3-D boolean phase) can then escalate it.
//!
//! # Intentional deductions are not clashes
//!
//! A girder that a column clips, or a wall pierced by its own opening, overlaps
//! the clipper *by design*. Pairs in a [`Clip`](crate::csg::CsgNode::Clip)
//! relationship (the base lists the other as a clipper) are therefore excluded
//! by default; set [`ClashOptions::include_clip_pairs`] to include them. This is
//! why clash checking takes the whole [`Model`] (which knows the clip DAG)
//! rather than a bare list of B-reps.
//!
//! # Local failure isolation
//!
//! A member that fails to evaluate does not abort the whole check: the pairs
//! involving it are reported with [`ClashError`] and every other pair is still
//! checked (`DESIGN.md` §2.3, §4.5).

mod aabb;
mod fem;
mod sleeve;

use std::collections::BTreeSet;

use crate::csg::{CsgNode, EvalError, StableId};
use crate::mass::signed_volume_checked;
use crate::model::Model;
use crate::tolerance::Tol;

use aabb::{aabb_of, Aabb};

pub use fem::{member_from_axis, AxisMemberError, ST_BRIDGE_MM_TO_M};
pub use sleeve::{
    check_sleeve, SleeveError, SleeveReport, SleeveRule, SleeveViolation, SleeveViolationKind,
};

/// What kind of interference a member pair exhibits.
///
/// `#[non_exhaustive]` so a future soft-clash (clearance) variant
/// (`docs/research/06-domain.md` §6.1) can be added in a semver-compatible way.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum ClashKind {
    /// The members share a positive volume — a true hard clash. The overlap
    /// volume (m³) is carried.
    HardClash {
        /// The exact shared volume in cubic metres.
        volume: f64,
    },
    /// The members touch (their boxes overlap and the exact shared volume is
    /// essentially zero) — a face/edge contact, not a volumetric clash.
    Touching,
    /// The members' bounding boxes overlap but no common prismatic direction
    /// exists, so the exact shared volume cannot be computed in this phase
    /// (`DESIGN.md` §4.2). Reported explicitly — never silently dropped — so the
    /// pair can be escalated to a 3-D boolean later (module docs).
    PotentialClash,
}

/// One detected interference between two members.
///
/// `a` and `b` are the members' [`StableId`]s with `a < b` (the pair is
/// unordered; the canonical ordering keeps results stable and de-duplicated).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub struct ClashResult {
    /// The lower-id member of the pair.
    pub a: StableId,
    /// The higher-id member of the pair.
    pub b: StableId,
    /// The kind of interference.
    pub kind: ClashKind,
}

/// A member pair whose interference could not be decided because one (or both)
/// of the members failed to evaluate. Reported alongside the successful results
/// so a single bad member does not hide the rest (local failure isolation).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub struct ClashError {
    /// The lower-id member of the pair.
    pub a: StableId,
    /// The higher-id member of the pair.
    pub b: StableId,
    /// The evaluation failure that prevented the pair from being checked.
    pub error: EvalError,
}

/// Options controlling the clash check.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub struct ClashOptions {
    /// Include pairs in a [`Clip`](crate::csg::CsgNode::Clip) relationship.
    ///
    /// `false` (the default) excludes intentional priority deductions and own
    /// openings from clash reporting, because a clipped girder *is meant* to
    /// overlap its column (module docs). Set `true` to audit those overlaps too.
    pub include_clip_pairs: bool,
    /// The minimum shared volume (m³) reported as a [`HardClash`].
    ///
    /// A volume at or below this is treated as a [`Touching`] contact. The
    /// default is derived from `tol` in [`ClashOptions::from_tol`]; this field
    /// lets a caller widen the dead band for noisy input.
    pub min_clash_volume: f64,
}

impl ClashOptions {
    /// Default options for a tolerance: exclude clip pairs, and set the
    /// hard-clash volume floor from `tol`.
    ///
    /// The floor is `tol.length³` scaled up by a small factor: a genuine clash in
    /// the building domain is at least millimetre-cubed scale, while numerical
    /// dust from a tolerant boolean of two touching solids is bounded by the
    /// length tolerance cubed. `1e3·tol.length³` (with the default `tol.length =
    /// 1e-6 m`, i.e. `1e-15 m³`) sits comfortably between the two without ever
    /// swallowing a real overlap.
    pub fn from_tol(tol: &Tol) -> Self {
        Self {
            include_clip_pairs: false,
            min_clash_volume: 1e3 * tol.length * tol.length * tol.length,
        }
    }
}

impl Default for ClashOptions {
    fn default() -> Self {
        Self::from_tol(&Tol::default())
    }
}

/// The full result of a clash check: the detected interferences and the pairs
/// that could not be evaluated.
#[derive(Debug, Clone, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub struct ClashCheck {
    /// The interferences found, in ascending `(a, b)` order.
    pub clashes: Vec<ClashResult>,
    /// The pairs that could not be checked because a member failed to evaluate.
    pub errors: Vec<ClashError>,
}

/// Check every member pair in `model` for interference (`DESIGN.md` §6-6).
///
/// Each member is evaluated once (so a failure is isolated to that member's
/// pairs), bounding boxes are computed, and every unordered pair is classified:
///
/// * boxes disjoint → not reported (coarse-phase reject);
/// * boxes overlap, common prismatic direction, volume > floor → [`HardClash`];
/// * boxes overlap, common direction, volume ≈ 0 → [`Touching`];
/// * boxes overlap, no common direction → [`PotentialClash`] (honest degeneracy);
/// * a [`Clip`](crate::csg::CsgNode::Clip) pair → excluded unless
///   [`ClashOptions::include_clip_pairs`].
///
/// # Complexity
///
/// The coarse phase compares all `O(n²)` pairs' boxes. A structural model has on
/// the order of hundreds to low-thousands of members per check scope (a frame,
/// a storey); `n²` box comparisons are microseconds at that scale and need no
/// spatial index. The expensive fine phase only runs on box-overlapping pairs,
/// which is near-linear in practice. A broad-phase tree is deferred until a
/// measured need (`docs/research/06-domain.md` §6: AABB pre-filter first).
///
/// The signature takes `&mut Model` to match the lazy-evaluation contract (a
/// check may trigger member re-evaluation), mirroring
/// [`takeoff`](crate::model::takeoff); the current implementation evaluates
/// through the immutable path and does not mutate the model.
///
/// # Errors
///
/// Never returns `Err` for the check as a whole — per-pair evaluation failures
/// are collected into [`ClashCheck::errors`] so the rest of the model is still
/// reported (local failure isolation).
pub fn clash_check(model: &mut Model, tol: &Tol, opts: &ClashOptions) -> ClashCheck {
    let ids: Vec<StableId> = model.ids().collect();

    // Evaluate each member once so a failure is isolated to its pairs. The fine
    // phase reads the occupancy leaf from the CSG, so only the AABB (built from
    // the evaluated B-rep) needs caching here — an `Ok(None)` is an empty member.
    let evaluated: Vec<(StableId, Result<Option<Aabb>, EvalError>)> = ids
        .iter()
        .map(|&id| (id, model.evaluate(id, tol).map(|brep| aabb_of(&brep))))
        .collect();

    let mut out = ClashCheck::default();

    for i in 0..evaluated.len() {
        for j in (i + 1)..evaluated.len() {
            let (a, ref res_a) = evaluated[i];
            let (b, ref res_b) = evaluated[j];
            // Ids are ascending in `evaluated`, so (a, b) is already canonical.

            // Skip intentional clip pairs unless asked to include them.
            if !opts.include_clip_pairs && is_clip_pair(model, a, b) {
                continue;
            }

            // A member that failed to evaluate makes the pair undecidable.
            let (box_a, box_b) = match (res_a, res_b) {
                (Err(e), _) | (Ok(_), Err(e)) => {
                    out.errors.push(ClashError {
                        a,
                        b,
                        error: e.clone(),
                    });
                    continue;
                }
                (Ok(box_a), Ok(box_b)) => (box_a, box_b),
            };

            // Both evaluated: compare boxes (coarse), then volume (fine). An empty
            // B-rep occupies no space — it never clashes.
            let (Some(box_a), Some(box_b)) = (box_a, box_b) else {
                continue;
            };
            if !box_a.overlaps(box_b, tol) {
                continue; // coarse reject
            }

            out.clashes.push(fine_classify(model, a, b, tol, opts));
        }
    }

    out
}

/// Classify a box-overlapping pair by its exact shared volume, or as a
/// [`PotentialClash`] when no common prismatic direction exists.
///
/// The exact shared volume of two members is the prismatic *intersection* of
/// their gross occupancy leaves — the single prisms they occupy, read straight
/// from their CSG trees. (Intentional clip overlaps are already filtered out
/// upstream, so the gross prism is the right operand here.) A member that is not
/// a single prism — a union, an oblique difference, anything outside the 2.5-D
/// fast path — has no single occupancy leaf, so the pair degrades to a
/// [`PotentialClash`] rather than a guessed volume.
fn fine_classify(
    model: &Model,
    a: StableId,
    b: StableId,
    tol: &Tol,
    opts: &ClashOptions,
) -> ClashResult {
    use crate::boolean::prismatic::{self, PrismError};

    let potential = ClashResult {
        a,
        b,
        kind: ClashKind::PotentialClash,
    };

    let (Some(ma), Some(mb)) = (model.get(a), model.get(b)) else {
        return potential;
    };
    let (Some(leaf_a), Some(leaf_b)) = (occupancy_leaf(ma.csg()), occupancy_leaf(mb.csg())) else {
        return potential;
    };

    match prismatic::intersection(&leaf_a, &leaf_b, tol) {
        Ok(inter) => {
            let vol = signed_volume_checked(&inter).map(f64::abs).unwrap_or(0.0);
            let kind = if vol > opts.min_clash_volume {
                ClashKind::HardClash { volume: vol }
            } else {
                ClashKind::Touching
            };
            ClashResult { a, b, kind }
        }
        // No common prismatic direction (H×H crossed, off-axis circle): the exact
        // volume is undecidable in this phase. Report honestly, never "no clash".
        Err(PrismError::NoCommonDirection)
        | Err(PrismError::CircularInvolved { .. })
        | Err(PrismError::ArcNotYetSupported) => ClashResult {
            a,
            b,
            kind: ClashKind::PotentialClash,
        },
        // Any other failure (complexity, internal) is also a degeneracy we cannot
        // resolve to a volume; surface it as potential rather than a wrong "clear".
        Err(_) => ClashResult {
            a,
            b,
            kind: ClashKind::PotentialClash,
        },
    }
}

/// `true` if `a` and `b` are in a clip relationship: one lists the other as a
/// [`Clip`](crate::csg::CsgNode::Clip) clipper anywhere in its CSG tree.
fn is_clip_pair(model: &Model, a: StableId, b: StableId) -> bool {
    clips_against(model, a, b) || clips_against(model, b, a)
}

/// `true` if member `base`'s CSG tree lists `clipper` as a clipper.
fn clips_against(model: &Model, base: StableId, clipper: StableId) -> bool {
    let Some(member) = model.get(base) else {
        return false;
    };
    let mut set = BTreeSet::new();
    collect_clippers(member.csg(), &mut set);
    set.contains(&clipper)
}

/// Collect every clipper id referenced anywhere in a CSG tree.
fn collect_clippers(node: &CsgNode, out: &mut BTreeSet<StableId>) {
    match node {
        CsgNode::Clip { base, clippers, .. } => {
            out.extend(clippers.iter().copied());
            collect_clippers(base, out);
        }
        CsgNode::OpeningSubtraction { base, .. } => collect_clippers(base, out),
        CsgNode::Union(nodes) => nodes.iter().for_each(|n| collect_clippers(n, out)),
        CsgNode::Difference { positive, negative } => {
            collect_clippers(positive, out);
            collect_clippers(negative, out);
        }
        CsgNode::Extrude { .. } => {}
    }
}

/// The gross occupancy [`ExtrudeLeaf`] of a member's CSG tree — the single prism
/// it occupies, ignoring its own openings and clips.
///
/// Mirrors the model layer's own occupancy descent: for clash purposes a member
/// displaces its full prism (a column with a sleeve still occupies its whole
/// cross-section against another member). Returns `None` when the member is not
/// ultimately a single extrusion (outside the 2.5-D fast path), so the fine
/// phase falls back to [`PotentialClash`].
fn occupancy_leaf(node: &CsgNode) -> Option<crate::boolean::prismatic::ExtrudeLeaf> {
    use crate::boolean::prismatic::ExtrudeLeaf;
    match node {
        CsgNode::Extrude {
            profile,
            origin,
            axis,
            length,
        } => Some(ExtrudeLeaf {
            profile: *profile,
            origin: *origin,
            axis: *axis,
            length: *length,
        }),
        CsgNode::OpeningSubtraction { base, .. } => occupancy_leaf(base),
        CsgNode::Clip { base, .. } => occupancy_leaf(base),
        _ => None,
    }
}
