//! Sleeve (beam-penetration) verification (`DESIGN.md` §6-6,
//! `docs/research/06-domain.md` §6.2).
//!
//! A service pipe passing through a beam is a *sleeve*: a circular void cut
//! through the beam web. Whether it is structurally admissible is governed by
//! diameter and position rules from the AIJ RC standard (`docs/research/06-domain.md`
//! §6.2):
//!
//! * **diameter** ≤ beam depth (せい) × a ratio (RC: 1/3);
//! * **end distance** — the sleeve centre must sit at least a multiple of the
//!   beam depth away from the beam ends (the column faces), to clear the plastic
//!   hinge zone;
//! * **vertical position** — the sleeve must stay near mid-depth; its centre must
//!   keep at least a margin from the top and bottom edges.
//!
//! The check reads the sleeves straight from the member's **CSG tree** — the
//! source of truth — rather than the evaluated B-rep: openings are identified by
//! their stable [`OpeningId`](crate::csg::OpeningId) (`DESIGN.md` §5.1, "開口は
//! stable id で識別"), and a circular opening carries its radius and centre line
//! exactly, with no need to re-derive them from faces. Every violation is
//! returned machine-readably (which opening, which rule, the measured value and
//! the limit), never as a bare boolean.

use crate::csg::{CsgNode, Opening, OpeningId, Profile2d};
use crate::math::{Point3, Vec3};
use crate::primitives::plane_basis;
use crate::tolerance::Tol;

/// The structural rule set a sleeve must satisfy, as depth-relative ratios.
///
/// Ratios are relative to the beam depth (せい) so one rule set applies to any
/// beam size. Defaults follow the RC standard
/// (`docs/research/06-domain.md` §6.2); steel beams use a looser diameter ratio
/// (1/2) and are configured by changing [`max_diameter_ratio`](Self::max_diameter_ratio).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub struct SleeveRule {
    /// Maximum sleeve diameter as a fraction of the beam depth. RC: `1/3`.
    pub max_diameter_ratio: f64,
    /// Minimum distance of the sleeve centre from either beam end, as a fraction
    /// of the beam depth. RC: `1.5` (centre ≥ 1.5·depth from the column face, to
    /// clear the hinge region).
    pub min_end_distance_ratio: f64,
    /// Minimum distance of the sleeve centre from the top **and** bottom edge, as
    /// a fraction of the beam depth — keeps the sleeve near mid-depth. A centre
    /// closer than this to either flange is a violation. Mid-depth placement
    /// corresponds to `0.5`; the rule of thumb keeps the centre out of the outer
    /// third, i.e. `≥ 1/3` from each edge.
    pub min_edge_distance_ratio: f64,
}

impl Default for SleeveRule {
    /// The reinforced-concrete sleeve rule (`docs/research/06-domain.md` §6.2):
    /// diameter ≤ depth/3, centre ≥ 1.5·depth from each end, centre ≥ depth/3
    /// from each flange.
    fn default() -> Self {
        Self {
            max_diameter_ratio: 1.0 / 3.0,
            min_end_distance_ratio: 1.5,
            min_edge_distance_ratio: 1.0 / 3.0,
        }
    }
}

/// Which rule a sleeve violated.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum SleeveViolationKind {
    /// The sleeve diameter exceeds `depth × max_diameter_ratio`.
    DiameterTooLarge,
    /// The sleeve centre is closer to a beam end than `depth × min_end_distance_ratio`.
    EndDistanceTooSmall,
    /// The sleeve centre is closer to the top or bottom flange than
    /// `depth × min_edge_distance_ratio`.
    EdgeDistanceTooSmall,
}

/// One rule a sleeve breaks, with the measured value and the limit it failed.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub struct SleeveViolation {
    /// The opening that violated a rule.
    pub opening: OpeningId,
    /// Which rule was broken.
    pub kind: SleeveViolationKind,
    /// The measured value (diameter, or centre distance), in metres.
    pub measured: f64,
    /// The limiting value the measurement was compared against, in metres.
    pub limit: f64,
}

/// The machine-readable result of checking a beam's sleeves.
///
/// A report with an empty [`violations`](Self::violations) list means every
/// sleeve is admissible. The beam geometry that the limits were derived from
/// (depth and clear span) is echoed back for traceability.
#[derive(Debug, Clone, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub struct SleeveReport {
    /// The beam depth (せい) the ratios were applied to, in metres.
    pub beam_depth: f64,
    /// The number of sleeves (circular openings) examined.
    pub sleeves_checked: usize,
    /// Every rule violation found, in opening order.
    pub violations: Vec<SleeveViolation>,
}

impl SleeveReport {
    /// `true` if no sleeve violated any rule.
    pub fn is_admissible(&self) -> bool {
        self.violations.is_empty()
    }
}

/// A reason the sleeve check could not run on a member.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum SleeveError {
    /// The member is not a beam the sleeve check understands: it is not a
    /// rectangular extrusion with circular openings subtracted. The check is
    /// declined rather than producing a misleading empty report.
    NotARectangularBeam,
}

/// Verify the circular sleeves of a beam member against `rule` (`DESIGN.md`
/// §6-6).
///
/// `beam` is the member's CSG tree — a rectangular [`Extrude`](CsgNode::Extrude)
/// with an [`OpeningSubtraction`](CsgNode::OpeningSubtraction) carrying the
/// circular sleeves. The beam *depth* (せい) is the rectangle half-height ×2 in
/// the profile's local `y` (the extruder's convention), and the clear span is
/// the extrusion length. For each circular opening the check measures:
///
/// * the diameter against `depth × max_diameter_ratio`;
/// * the centre's axial distance from each beam end against
///   `depth × min_end_distance_ratio`;
/// * the centre's transverse distance from the top and bottom flange against
///   `depth × min_edge_distance_ratio`.
///
/// Non-circular openings are ignored (they are not sleeves). The geometry is
/// taken from the CSG, so the opening's stable id flows straight into each
/// [`SleeveViolation`].
///
/// # Errors
///
/// [`SleeveError::NotARectangularBeam`] if `beam` is not a rectangular
/// extrusion (possibly wrapped in an opening subtraction).
pub fn check_sleeve(
    beam: &CsgNode,
    rule: &SleeveRule,
    tol: &Tol,
) -> Result<SleeveReport, SleeveError> {
    let (base, openings) = split_beam(beam).ok_or(SleeveError::NotARectangularBeam)?;
    let BeamGeom {
        depth,
        length,
        origin,
        axis,
        u,
        v,
    } = beam_geom(base).ok_or(SleeveError::NotARectangularBeam)?;

    let mut report = SleeveReport {
        beam_depth: depth,
        sleeves_checked: 0,
        violations: Vec::new(),
    };

    let dia_limit = depth * rule.max_diameter_ratio;
    let end_limit = depth * rule.min_end_distance_ratio;
    let edge_limit = depth * rule.min_edge_distance_ratio;
    let half_depth = depth / 2.0;

    for (id, Opening { shape }) in openings {
        // Only circular openings are sleeves.
        let Some((radius, centre)) = circular_opening(shape) else {
            continue;
        };
        report.sleeves_checked += 1;
        let diameter = 2.0 * radius;

        // Diameter rule.
        if diameter > dia_limit + tol.length {
            report.violations.push(SleeveViolation {
                opening: *id,
                kind: SleeveViolationKind::DiameterTooLarge,
                measured: diameter,
                limit: dia_limit,
            });
        }

        // Axial position of the centre relative to the beam start (origin).
        let rel = centre - origin;
        let t = rel.dot(axis); // along the beam
        let end_dist = t.min(length - t); // nearer of the two ends
        if end_dist < end_limit - tol.length {
            report.violations.push(SleeveViolation {
                opening: *id,
                kind: SleeveViolationKind::EndDistanceTooSmall,
                measured: end_dist,
                limit: end_limit,
            });
        }

        // Transverse (depth-direction) position. The profile's local `y` (the
        // extruder's `v` for a `+z` extrusion is `−x`, etc.) carries the beam
        // depth; the half-height is `depth/2`. The centre's offset from mid-depth
        // is `rel·v_depth`; its distance to the nearer flange is
        // `half_depth − |offset|`.
        let offset = rel.dot(v).abs();
        let edge_dist = half_depth - offset;
        if edge_dist < edge_limit - tol.length {
            report.violations.push(SleeveViolation {
                opening: *id,
                kind: SleeveViolationKind::EdgeDistanceTooSmall,
                measured: edge_dist,
                limit: edge_limit,
            });
        }
        // `u` is the width direction; not constrained by the sleeve rules but
        // bound to silence the unused read.
        let _ = u;
    }

    Ok(report)
}

/// The beam's rectangular base node and its openings, when `beam` is a
/// rectangular extrusion optionally wrapped in an opening subtraction.
fn split_beam(beam: &CsgNode) -> Option<(&CsgNode, &[(OpeningId, Opening)])> {
    match beam {
        CsgNode::OpeningSubtraction { base, openings } => Some((base, openings.as_slice())),
        CsgNode::Extrude { .. } => Some((beam, &[])),
        _ => None,
    }
}

/// The geometric description of a rectangular beam, in world coordinates.
struct BeamGeom {
    /// Beam depth (せい), in metres — the rectangle's full local-`y` extent.
    depth: f64,
    /// Clear span (extrusion length), in metres.
    length: f64,
    /// Beam start (bottom-cap centre / one end), world coordinates.
    origin: Point3,
    /// Unit beam axis (direction of the span).
    axis: Vec3,
    /// In-plane width direction (unit), profile local `x`.
    u: Vec3,
    /// In-plane depth direction (unit), profile local `y` — the せい direction.
    v: Vec3,
}

/// Read a rectangular beam's geometry from its extrude node.
///
/// The profile is `Rect { half_w, half_h }`; the extruder places the profile
/// `u` axis along `half_w` (width) and `v` along `half_h` (depth/せい). The frame
/// `(u, v)` is the same [`plane_basis`] the extruder uses, so the depth direction
/// here matches the built geometry.
fn beam_geom(base: &CsgNode) -> Option<BeamGeom> {
    let CsgNode::Extrude {
        profile,
        origin,
        axis,
        length,
    } = base
    else {
        return None;
    };
    let Profile2d::Rect { half_w: _, half_h } = profile else {
        return None;
    };
    let unit = axis.try_unit()?;
    let (u, v) = plane_basis(unit);
    Some(BeamGeom {
        depth: 2.0 * half_h,
        length: *length,
        origin: *origin,
        axis: unit.as_vec(),
        u,
        v,
    })
}

/// The `(radius, world centre-line origin)` of a circular opening, or `None` if
/// the opening is not a circular extrusion.
///
/// The centre is the opening extrusion's own origin (the cylinder axis passes
/// through it). For the position checks we need a point on that axis closest to
/// the beam; the opening origin lies on the axis, and the beam is thin in the
/// sleeve-axis direction, so the origin is a faithful sleeve centre.
fn circular_opening(shape: &CsgNode) -> Option<(f64, Point3)> {
    let CsgNode::Extrude {
        profile, origin, ..
    } = shape
    else {
        return None;
    };
    let Profile2d::Circle { radius } = profile else {
        return None;
    };
    Some((*radius, *origin))
}
