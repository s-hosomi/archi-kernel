//! Prismatic (2.5-D) boolean reduction — the building domain's common path.
//!
//! The overwhelming majority of architectural booleans are a prismatic solid
//! minus / union / intersection another prismatic solid that share a common
//! prismatic direction (`DESIGN.md` §4.1–4.2). When such a direction exists the
//! 3-D boolean collapses to a 2-D polygon boolean per axial band, which the
//! proven [`poly2d`](crate::boolean::poly2d) engine solves robustly, plus a
//! one-shot watertight 3-D reassembly.
//!
//! ## Pipeline
//!
//! 1. The `detect` step finds a common prismatic direction of the two extruded
//!    leaves, builds the shared 2-D frame, and reprofiles both into a 2-D region
//!    plus an axial interval.
//! 2. The `arrange` step overlays every region into **one** global 2-D
//!    arrangement (atomic cells + edges), split once and reused for every band —
//!    the keystone that makes the result watertight without band stitching.
//! 3. The `build` step walks the `cell × band` voxels, emits wall and interface
//!    faces from the shared segmentation, splits disconnected pieces into
//!    separate solids, and validates at
//!    [`Full`](crate::topo::ValidateLevel::Full).
//!
//! ## Entry points
//!
//! [`difference`], [`union`] and [`intersection`] take the CSG extrude fields of
//! the operands directly. [`opening_subtraction`] removes many openings in one
//! multi-band pass (`DESIGN.md` §4.5, web-ifc's `BOOLEAN_UNION_THRESHOLD`
//! lesson). Failures are returned as a fine-grained [`PrismError`] that the CSG
//! layer maps onto [`EvalError`](crate::csg::EvalError).

mod arrange;
mod build;
pub(crate) mod detect;
pub(crate) mod error;

pub use error::{Operand, PrismError};

use crate::boolean::poly2d::{self, Op};
use crate::brep::Brep;
use crate::csg::Profile2d;
use crate::math::{Point3, Vec3};
use crate::tolerance::Tol;

use detect::{detect, Leaf};

/// Default boundary-complexity budget (`DESIGN.md` §4.5).
///
/// The build's cost is bounded by `edges·bands + cells·levels`. A single member
/// boolean in the building domain stays in the low hundreds (a few dozen profile
/// edges × a handful of bands); `100_000` isolates a pathological input (e.g. a
/// member accidentally given thousands of openings) long before it grinds, while
/// never tripping on a legitimate wall-with-windows. It is a guard, not a tuning
/// knob — chosen two orders of magnitude above any real member.
pub const DEFAULT_BUDGET: usize = 100_000;

/// One extruded operand described by its CSG [`Extrude`](crate::csg::CsgNode)
/// fields.
#[derive(Debug, Clone, Copy)]
pub struct ExtrudeLeaf {
    /// The cross-section profile.
    pub profile: Profile2d,
    /// World origin (bottom-cap centre).
    pub origin: Point3,
    /// Extrusion direction (need not be unit).
    pub axis: Vec3,
    /// Extrusion length in metres.
    pub length: f64,
}

impl ExtrudeLeaf {
    fn to_leaf(self) -> Leaf {
        Leaf::new(self.profile, self.origin, self.axis, self.length)
    }
}

/// `positive − negative` of two extruded leaves, as a watertight [`Brep`].
pub fn difference(
    positive: &ExtrudeLeaf,
    negative: &ExtrudeLeaf,
    tol: &Tol,
) -> Result<Brep, PrismError> {
    binary(positive, negative, Op::Difference, tol)
}

/// `positive − negative` with an explicit complexity budget.
///
/// Same as [`difference`] but the fail-safe budget (`DESIGN.md` §4.5) is
/// caller-supplied. Exceeding it yields [`PrismError::ComplexityLimit`]. The
/// default path uses [`DEFAULT_BUDGET`]; this entry point exists so callers (and
/// the complexity-limit regression test) can dial the guard.
pub fn difference_with_budget(
    positive: &ExtrudeLeaf,
    negative: &ExtrudeLeaf,
    tol: &Tol,
    budget: usize,
) -> Result<Brep, PrismError> {
    let la = positive.to_leaf();
    let lb = negative.to_leaf();
    let (frame, pa, pb) = detect(&la, &lb, tol)?;
    build::build(&frame, &pa, &pb, Op::Difference, tol, budget)
}

/// `a ∩ b` of two extruded leaves.
pub fn intersection(a: &ExtrudeLeaf, b: &ExtrudeLeaf, tol: &Tol) -> Result<Brep, PrismError> {
    binary(a, b, Op::Intersection, tol)
}

/// `a ∪ b` of two extruded leaves.
pub fn union_pair(a: &ExtrudeLeaf, b: &ExtrudeLeaf, tol: &Tol) -> Result<Brep, PrismError> {
    binary(a, b, Op::Union, tol)
}

/// `union(leaves)` of extruded leaves sharing a common prismatic direction.
///
/// An empty input is an empty B-rep and a single input is that leaf extruded.
/// The 2.5-D fast path covers the **two-leaf** union directly; three or more
/// leaves would require a `Brep × Brep` boolean (the accumulator is no longer a
/// single extruded leaf), which is outside this phase — reported as
/// [`PrismError::NoCommonDirection`] so the caller falls back cleanly.
pub fn union(leaves: &[ExtrudeLeaf], tol: &Tol) -> Result<Brep, PrismError> {
    match leaves {
        [] => Ok(Brep::new()),
        [single] => single_leaf(single, tol),
        [a, b] => union_pair(a, b, tol),
        _ => Err(PrismError::NoCommonDirection),
    }
}

/// Subtract many openings from a base in one pass (`DESIGN.md` §4.5).
///
/// All openings and the base are reprofiled onto the **same** shared frame and
/// fed to the single multi-band arrangement build, with the residency rule
/// `base ∧ ¬(opening₀ ∨ opening₁ ∨ …)`. Openings that overlap, touch, or sit at
/// different axial positions are handled uniformly — the union is taken in 3-D
/// (per band) rather than by collapsing the openings to one cross-section, which
/// would wrongly bridge openings separated along the prism axis. Every opening
/// must share the base's prismatic direction, else this is not 2.5-D and
/// [`PrismError::NoCommonDirection`] is returned.
pub fn opening_subtraction(
    base: &ExtrudeLeaf,
    openings: &[ExtrudeLeaf],
    tol: &Tol,
) -> Result<Brep, PrismError> {
    if openings.is_empty() {
        return single_leaf(base, tol);
    }
    if openings.len() == 1 {
        return difference(base, &openings[0], tol);
    }

    let base_leaf = base.to_leaf();
    let first = openings[0].to_leaf();
    let (frame, base_op, first_op) = detect(&base_leaf, &first, tol)?;

    // operands[0] = base; operands[1..] = each opening, all on the shared frame.
    let mut operands = Vec::with_capacity(openings.len() + 1);
    operands.push(base_op);
    operands.push(first_op);
    for opening in &openings[1..] {
        let leaf = opening.to_leaf();
        let (other_frame, _base_again, op) = detect(&base_leaf, &leaf, tol)?;
        // Every opening must reduce along the same direction as the first; if a
        // later opening picks a different common direction the set is not a
        // single 2.5-D problem.
        if frame.d.cross(other_frame.d).norm() > tol.angular {
            return Err(PrismError::NoCommonDirection);
        }
        operands.push(op);
    }

    // Residency: keep where the base is present and no opening covers the voxel.
    build::build_combined(
        &frame,
        &operands,
        |flags| flags[0] && !flags[1..].iter().any(|&f| f),
        tol,
        DEFAULT_BUDGET,
    )
}

/// Extrude a single leaf to a watertight B-rep (the degenerate "boolean" of one
/// operand). Reuses the standard extruder so the result matches Phase 2.
fn single_leaf(leaf: &ExtrudeLeaf, tol: &Tol) -> Result<Brep, PrismError> {
    use crate::build::extrude;
    use crate::primitives::Line3;
    let line = Line3::new(leaf.origin, leaf.axis).map_err(|_| {
        PrismError::Poly2(poly2d::Poly2Error::Internal {
            what: "degenerate extrusion axis",
        })
    })?;
    extrude(&leaf.profile, &line, leaf.length, tol)
        .map_err(|e| PrismError::Poly2(poly2d::Poly2Error::Internal { what: leak(&e) }))
}

fn leak(_e: &crate::error::KernelError) -> &'static str {
    "extrusion construction failed"
}

/// Shared binary driver: detect the common direction, then build for `op`.
fn binary(a: &ExtrudeLeaf, b: &ExtrudeLeaf, op: Op, tol: &Tol) -> Result<Brep, PrismError> {
    let la = a.to_leaf();
    let lb = b.to_leaf();
    let (frame, pa, pb) = detect(&la, &lb, tol)?;
    build::build(&frame, &pa, &pb, op, tol, DEFAULT_BUDGET)
}
