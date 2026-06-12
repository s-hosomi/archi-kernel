//! Minimal section extraction (`DESIGN.md` §6-3, Phase 4 does the polishing).
//!
//! A section is the set of loops where a cutting plane meets a solid — the same
//! loops the half-space cut builds as cap faces. [`section`] reuses the cut
//! machinery (it runs an actual cut and reads back the cap loops), so the
//! outer / hole correspondence is exactly the one the cut produces: a through
//! hole gives an outer loop with one inner loop.
//!
//! The return type is intentionally minimal — full 2-D typing and tidy-up are
//! Phase 4. Each loop is returned as its ordered ring of 3-D points plus the
//! same ring projected to the cut plane's local 2-D frame, and holes are nested
//! under their outer loop so the correspondence is not lost.

use crate::boolean::{cut, CutResult, KeepSide};
use crate::brep::Brep;
use crate::csg::EvalError;
use crate::geom::SurfaceGeom;
use crate::math::Point3;
use crate::primitives::{plane_basis, Plane};
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::{Loop, Solid};

/// A single section loop: its 3-D point ring and the same ring in the cut
/// plane's local 2-D coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionLoop {
    /// The loop's vertices in order, in world 3-D coordinates (metres).
    pub points_3d: Vec<Point3>,
    /// The same vertices projected to the cut plane's local `(s, t)` frame.
    pub points_2d: Vec<[f64; 2]>,
}

/// A section outline: one outer loop and its hole loops.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionOutline {
    /// The outer boundary loop.
    pub outer: SectionLoop,
    /// The hole loops nested inside `outer`.
    pub holes: Vec<SectionLoop>,
}

/// The full section: every disjoint outline where the plane meets the solid.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SectionLoops {
    /// The disjoint outlines (each an outer loop plus its holes).
    pub outlines: Vec<SectionOutline>,
}

impl SectionLoops {
    /// Total number of loops (outer + holes) across all outlines.
    pub fn loop_count(&self) -> usize {
        self.outlines.iter().map(|o| 1 + o.holes.len()).sum()
    }
}

/// Extract the section of `solid` by `plane`.
///
/// Reuses the half-space cut: the cut's cap faces are precisely the section
/// outlines, with holes already nested. Returns the loops with their 3-D and
/// plane-local 2-D coordinates.
///
/// # Errors
///
/// Propagates [`EvalError`] from the underlying cut (the cut self-validates).
pub fn section(
    brep: &Brep,
    solid: Id<Solid>,
    plane: &Plane,
    tol: &Tol,
) -> Result<SectionLoops, EvalError> {
    let result = cut(brep, solid, plane, KeepSide::Below, tol)?;
    let CutResult::Cut {
        brep: cut_brep,
        caps,
    } = result
    else {
        // The plane missed the solid (all kept / empty): no section.
        return Ok(SectionLoops::default());
    };

    // 2-D frame on the cut plane.
    let (u, v) = plane_basis(plane.normal());
    let origin = plane.point();
    let project = |p: Point3| {
        let d = p - origin;
        [d.dot(u), d.dot(v)]
    };

    let mut outlines = Vec::new();
    for &cap in &caps {
        let Some(face) = cut_brep.topo.faces.get(cap) else {
            continue;
        };
        // Confirm this is a cap on the cut plane (planar surface).
        if !matches!(
            cut_brep.geom.surface(face.surface),
            Some(SurfaceGeom::Plane(_))
        ) {
            continue;
        }
        let outer = loop_to_section(&cut_brep, face.outer, &project);
        let holes: Vec<SectionLoop> = face
            .inners
            .iter()
            .filter_map(|&l| loop_to_section(&cut_brep, l, &project))
            .collect();
        if let Some(outer) = outer {
            outlines.push(SectionOutline { outer, holes });
        }
    }

    Ok(SectionLoops { outlines })
}

/// Convert a loop in the cut B-rep to a [`SectionLoop`].
fn loop_to_section(
    brep: &Brep,
    loop_id: Id<Loop>,
    project: &impl Fn(Point3) -> [f64; 2],
) -> Option<SectionLoop> {
    let lp = brep.topo.loops.get(loop_id)?;
    let mut points_3d = Vec::with_capacity(lp.half_edges.len());
    for &he_id in &lp.half_edges {
        let he = brep.topo.half_edges.get(he_id)?;
        let vert = brep.topo.vertices.get(he.start)?;
        let p = brep.geom.point(vert.point)?.as_point()?;
        points_3d.push(p);
    }
    if points_3d.len() < 3 {
        return None;
    }
    let points_2d = points_3d.iter().map(|&p| project(p)).collect();
    Some(SectionLoop {
        points_3d,
        points_2d,
    })
}
