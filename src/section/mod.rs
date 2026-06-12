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
use crate::geom::{CurveGeom, SurfaceGeom};
use crate::math::Point3;
use crate::primitives::{plane_basis, Plane};
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::{Loop, Solid};

/// One boundary edge of a section loop, in the cut plane's local 2-D frame.
///
/// A section boundary is in general a mix of straight segments (polygonal
/// cross-sections — walls, beams, H-sections) and circular arcs (round members:
/// a vertical section of a circular column yields a loop of two semicircular
/// arcs). Ellipse arcs (oblique cuts of a cylinder) are not produced by the
/// vertical/perpendicular section paths Phase 4 draws and are reported as their
/// chord segment here; full ellipse support is Phase 5.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum SectionEdge {
    /// A straight segment from `start` to `end` (2-D, cut-plane frame).
    Line {
        /// Segment start.
        start: [f64; 2],
        /// Segment end.
        end: [f64; 2],
    },
    /// A circular arc, centre/radius in the cut-plane frame, swept from
    /// `start_angle` to `end_angle` (radians, the curve's own parameterisation),
    /// with `start`/`end` the arc endpoints (2-D, cut-plane frame).
    Arc {
        /// Arc centre in the cut-plane frame.
        center: [f64; 2],
        /// Arc radius (metres).
        radius: f64,
        /// Start angle (radians).
        start_angle: f64,
        /// End angle (radians).
        end_angle: f64,
        /// Arc start point.
        start: [f64; 2],
        /// Arc end point.
        end: [f64; 2],
    },
}

/// A single section loop, described both as a typed edge list (line segments and
/// circular arcs) and — for straightforward polygonal consumption — its ordered
/// vertex ring in 3-D and in the cut plane's local 2-D `(s, t)` frame.
///
/// For a polygonal loop the vertex ring has one point per edge. For a round
/// member the loop is two arcs over two seam vertices, so the vertex ring has two
/// points while `edges` carries the two arcs.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionLoop {
    /// The loop's vertices in order, in world 3-D coordinates (metres).
    pub points_3d: Vec<Point3>,
    /// The same vertices projected to the cut plane's local `(s, t)` frame.
    pub points_2d: Vec<[f64; 2]>,
    /// The loop's boundary as a typed edge list (segments + arcs) in the 2-D
    /// frame, in loop order.
    pub edges: Vec<SectionEdge>,
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
/// plane-local 2-D coordinates, plus a typed [`SectionEdge`] list (straight
/// segments and circular arcs), so a round member sections to an arc loop.
///
/// # Coplanar-coincidence convention
///
/// When the cutting plane is **coplanar with a face of the solid** (the everyday
/// "section on a grid line / slab top" case), that coincident face is the kept
/// material's lid and is *not* re-emitted as a separate section outline — the
/// section reports the cross-section of the material the plane passes *through*,
/// not the faces it lies *in*. (Whether Phase 4 should instead include a
/// coincident face's outline as a section profile is a Phase-4 output-spec
/// decision per `DESIGN.md` §6-3 item 3; the current behaviour is "a coincident
/// face is not part of the section".)
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

/// Convert a loop in the cut B-rep to a [`SectionLoop`], reading each half-edge's
/// curve so straight edges become [`SectionEdge::Line`] and circular arcs become
/// [`SectionEdge::Arc`]. A round member's section (two semicircular arcs over two
/// seam vertices) is therefore preserved rather than dropped.
fn loop_to_section(
    brep: &Brep,
    loop_id: Id<Loop>,
    project: &impl Fn(Point3) -> [f64; 2],
) -> Option<SectionLoop> {
    let lp = brep.topo.loops.get(loop_id)?;
    let n = lp.half_edges.len();
    if n == 0 {
        return None;
    }
    let mut points_3d = Vec::with_capacity(n);
    let mut edges: Vec<SectionEdge> = Vec::with_capacity(n);
    let mut has_arc = false;
    for i in 0..n {
        let he = brep.topo.half_edges.get(lp.half_edges[i])?;
        let next = brep.topo.half_edges.get(lp.half_edges[(i + 1) % n])?;
        let start = brep
            .geom
            .point(brep.topo.vertices.get(he.start)?.point)?
            .as_point()?;
        let end = brep
            .geom
            .point(brep.topo.vertices.get(next.start)?.point)?
            .as_point()?;
        points_3d.push(start);
        let (s2, e2) = (project(start), project(end));
        match brep.geom.curve(he.curve)? {
            CurveGeom::Circle(c) => {
                has_arc = true;
                let center = project(c.center());
                let ang = |p: [f64; 2]| (p[1] - center[1]).atan2(p[0] - center[0]);
                edges.push(SectionEdge::Arc {
                    center,
                    radius: c.radius(),
                    start_angle: ang(s2),
                    end_angle: ang(e2),
                    start: s2,
                    end: e2,
                });
            }
            // Ellipse arcs (oblique cylinder cuts) are reported as their chord
            // until Phase 5; their vertical/perpendicular sections never reach
            // this path.
            _ => edges.push(SectionEdge::Line { start: s2, end: e2 }),
        }
    }
    // A polygonal loop needs at least three vertices; an arc loop (a round
    // member) can legitimately have only two seam vertices, so allow it when an
    // arc edge is present.
    if points_3d.len() < 3 && !has_arc {
        return None;
    }
    let points_2d = points_3d.iter().map(|&p| project(p)).collect();
    Some(SectionLoop {
        points_3d,
        points_2d,
        edges,
    })
}
