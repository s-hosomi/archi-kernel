//! Section drawing — the cut plane's cross-section as tidy holed profiles
//! (`DESIGN.md` §6-3, Phase 4).
//!
//! A section is the cross-section a cutting plane makes through a solid: the
//! plan (伏図) at a column's mid-height, the elevation (軸組図) through a grid
//! line, the section (断面図) through a slab. The output is a list of
//! [`SectionProfile`]s — each a CCW outer loop plus its CW holes — expressed in
//! the cut plane's own 2-D `(u, v)` frame, together with the [`SectionFrame`]
//! that maps those 2-D coordinates back to world 3-D.
//!
//! # How it is computed
//!
//! [`section`] reuses the half-space [`cut`](crate::boolean::cut): the cut's cap
//! faces are precisely the cross-section of the material the plane passes
//! *through*. Every loop of every cap (plus the coplanar lid loops, see below)
//! is collected into one flat pool and **re-nested globally** by exact
//! orient2d-based containment and signed area (with a circular-segment
//! correction for arcs), so a profile's outer / hole / island-in-hole structure
//! is correct even when the cut produced several disjoint cap faces (the fused
//! sleeve / compound-void case that the per-cap grouping got wrong, the Phase 3c
//! known limitation).
//!
//! # Coplanar-coincidence convention
//!
//! When the cutting plane is **coplanar with a face of the solid** — the
//! everyday "plan on a slab top" or "elevation on a grid line" case — the rule
//! is:
//!
//! > A face whose outward normal **agrees with the cut plane's `+normal`** is
//! > included in the section (it is the lid of the material on the kept side);
//! > a coincident face whose normal opposes `+normal` is not.
//!
//! Concretely: a horizontal section taken **exactly on a slab's top face**
//! (slab top outward normal `+z`, cut normal `+z`) reports the slab's plan
//! (openings included as holes). The same plane taken on the slab's **bottom**
//! face (outward normal `−z`) reports nothing — the material is on the other
//! side. This draws "the lid of the material that is there", which is what a
//! plan/elevation on a coincident plane is expected to show. The convention is
//! pinned by `tests/section_drawing.rs`.

use crate::boolean::support::{loop_signed_area_2d, point_in_loop_2d};
use crate::boolean::{cut, CutResult, KeepSide};
use crate::brep::Brep;
use crate::csg::{EvalError, Member};
use crate::geom::{CurveGeom, SurfaceGeom};
use crate::math::{Point3, Vec3};
use crate::primitives::{plane_basis, Plane};
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::{Face, Loop, Sense, Solid};

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
/// vertex ring in 3-D and in the cut plane's local 2-D `(u, v)` frame.
///
/// For a polygonal loop the vertex ring has one point per edge. For a round
/// member the loop is two arcs over two seam vertices, so the vertex ring has two
/// points while `edges` carries the two arcs.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionLoop {
    /// The loop's vertices in order, in world 3-D coordinates (metres).
    pub points_3d: Vec<Point3>,
    /// The same vertices projected to the cut plane's local `(u, v)` frame.
    pub points_2d: Vec<[f64; 2]>,
    /// The loop's boundary as a typed edge list (segments + arcs) in the 2-D
    /// frame, in loop order.
    pub edges: Vec<SectionEdge>,
}

/// The cut plane's local 2-D coordinate frame.
///
/// `(u, v)` is an orthonormal basis of the plane with `u × v = normal`, and
/// `origin` is the plane's reference point. A world point `p` on the plane has
/// 2-D coordinates `(d·u, d·v)` with `d = p − origin`; conversely a 2-D point
/// `(s, t)` maps back to 3-D via [`to_3d`](Self::to_3d) as `origin + s·u + t·v`.
/// The basis comes from the kernel's deterministic plane-basis seed rule (the
/// same one `Circle3` and the extruder use), so the frame is reproducible from
/// the plane alone — a caller
/// can re-derive it and round-trip every 2-D coordinate back to its world point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SectionFrame {
    /// The plane's reference point (world 3-D).
    pub origin: Point3,
    /// In-plane `u` axis (unit, world 3-D).
    pub u: Vec3,
    /// In-plane `v` axis (unit, world 3-D), with `u × v = plane normal`.
    pub v: Vec3,
}

impl SectionFrame {
    /// Build the frame from a cutting plane.
    pub fn new(plane: &Plane) -> Self {
        let (u, v) = plane_basis(plane.normal());
        Self {
            origin: plane.point(),
            u,
            v,
        }
    }

    /// Project a world 3-D point onto the frame's 2-D `(u, v)` coordinates.
    #[inline]
    pub fn project(&self, p: Point3) -> [f64; 2] {
        let d = p - self.origin;
        [d.dot(self.u), d.dot(self.v)]
    }

    /// Map a 2-D `(s, t)` coordinate back to its world 3-D point on the plane.
    #[inline]
    pub fn to_3d(&self, p: [f64; 2]) -> Point3 {
        self.origin + self.u * p[0] + self.v * p[1]
    }
}

/// One section profile: a CCW outer boundary loop and its CW hole loops.
///
/// A hole that itself contains an island (material inside a hole) is reported as
/// a *separate* [`SectionProfile`], so every profile here is a simple
/// ring-with-holes — exactly the holed `Polygon2d` a drawing consumes.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionProfile {
    /// The outer boundary loop (CCW in the frame's `(u, v)`).
    pub outer: SectionLoop,
    /// The hole loops directly nested inside `outer` (CW in the frame).
    pub holes: Vec<SectionLoop>,
}

/// The full result of a section: the 2-D frame plus every disjoint profile.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SectionResult {
    /// The cut plane's 2-D frame (for restoring 2-D coordinates to world 3-D).
    pub frame: Option<SectionFrame>,
    /// The disjoint section profiles (each an outer loop plus its holes).
    pub profiles: Vec<SectionProfile>,
}

impl SectionResult {
    /// Total number of loops (outer + holes) across all profiles.
    pub fn loop_count(&self) -> usize {
        self.profiles.iter().map(|p| 1 + p.holes.len()).sum()
    }
}

/// Extract the section of `solid` of `brep` by `plane`.
///
/// Reuses the half-space cut and re-nests all section loops globally (see the
/// module docs). The coplanar-coincidence convention there governs faces that
/// lie exactly in the cut plane. The result's profiles are in the cut plane's
/// local 2-D frame; [`SectionResult::frame`] maps them back to world 3-D.
///
/// # Errors
///
/// Propagates [`EvalError`] from the underlying cut (the cut self-validates).
pub fn section(
    brep: &Brep,
    solid: Id<Solid>,
    plane: &Plane,
    tol: &Tol,
) -> Result<SectionResult, EvalError> {
    let frame = SectionFrame::new(plane);
    let project = |p: Point3| frame.project(p);

    // Run the cut. Its caps are the cross-section of the through-material; its
    // kept coplanar lid faces are the coincident-face contribution.
    let result = cut(brep, solid, plane, KeepSide::Below, tol)?;
    let (cut_brep, cap_faces): (Brep, Vec<Id<Face>>) = match result {
        CutResult::Cut { brep, caps } => (brep, caps),
        // The plane missed (AllKept) or removed everything (Empty): no through
        // cross-section. A coincident lid can still exist on an all-kept solid,
        // so fall through with the *input* brep and no cap faces; the coplanar
        // scan below picks up any lid.
        CutResult::AllKept { brep } => (brep, Vec::new()),
        CutResult::Empty => return Ok(SectionResult::default()),
    };

    // Collect every loop that belongs to the section: cap-face loops plus the
    // loops of any kept coplanar lid face whose outward normal agrees with the
    // cut plane's +normal (the coplanar convention).
    let mut pool: Vec<SectionLoop> = Vec::new();
    let cap_set: std::collections::HashSet<Id<Face>> = cap_faces.iter().copied().collect();

    for &cap in &cap_faces {
        collect_face_loops(&cut_brep, cap, &project, &mut pool);
    }
    collect_coplanar_lids(&cut_brep, plane, &cap_set, tol, &project, &mut pool);

    let profiles = nest_profiles(pool);
    Ok(SectionResult {
        frame: Some(frame),
        profiles,
    })
}

/// Section every member of `members`, returning one [`SectionResult`] per member
/// (伏図 API: the plan is the overlay of each member's section).
///
/// Each member is evaluated through its cache ([`Member::brep`]) and sectioned by
/// the shared `plane`. A member whose evaluation **fails** is *not* skipped: its
/// slot carries the member-local [`EvalError`] (local failure isolation,
/// `DESIGN.md` §4.5), so the caller sees exactly which member could not be drawn
/// while the rest still produce their profiles. The returned vector is parallel
/// to `members` (one entry each, in order).
pub fn section_members(
    members: &mut [Member],
    plane: &Plane,
    tol: &Tol,
) -> Vec<Result<SectionResult, EvalError>> {
    members
        .iter_mut()
        .map(|m| {
            let brep = m.brep(tol)?.clone();
            let Some(&solid) = brep.solids.first() else {
                return Ok(SectionResult::default());
            };
            section(&brep, solid, plane, tol)
        })
        .collect()
}

/// Push each loop of `face` (outer + inners) into the pool as a [`SectionLoop`].
fn collect_face_loops(
    brep: &Brep,
    face_id: Id<Face>,
    project: &impl Fn(Point3) -> [f64; 2],
    pool: &mut Vec<SectionLoop>,
) {
    let Some(face) = brep.topo.faces.get(face_id) else {
        return;
    };
    // Caps and coplanar lids are planar; a non-planar face is never a section.
    if !matches!(brep.geom.surface(face.surface), Some(SurfaceGeom::Plane(_))) {
        return;
    }
    let mut loops = vec![face.outer];
    loops.extend(face.inners.iter().copied());
    for lid in loops {
        if let Some(sl) = loop_to_section(brep, lid, project) {
            pool.push(sl);
        }
    }
}

/// Scan the cut result for kept coplanar lid faces and add their loops.
///
/// A lid face lies in the cut plane with outward normal along `+normal`; the cut
/// keeps exactly those (for `KeepSide::Below`), so any planar face of the kept
/// brep that is coincident with the cut plane *and* faces `+normal` is a lid. We
/// skip faces already counted as caps so a face is never emitted twice.
fn collect_coplanar_lids(
    brep: &Brep,
    plane: &Plane,
    cap_set: &std::collections::HashSet<Id<Face>>,
    tol: &Tol,
    project: &impl Fn(Point3) -> [f64; 2],
    pool: &mut Vec<SectionLoop>,
) {
    let plane_n = plane.normal().as_vec();
    for &solid_id in &brep.solids {
        let Some(solid) = brep.topo.solids.get(solid_id) else {
            continue;
        };
        for &shell_id in &solid.shells {
            let Some(shell) = brep.topo.shells.get(shell_id) else {
                continue;
            };
            for &face_id in &shell.faces {
                if cap_set.contains(&face_id) {
                    continue;
                }
                let Some(face) = brep.topo.faces.get(face_id) else {
                    continue;
                };
                let Some(SurfaceGeom::Plane(face_plane)) = brep.geom.surface(face.surface) else {
                    continue;
                };
                // Coplanar with the cut plane?
                let cross = face_plane.normal().as_vec().cross(plane_n);
                let coplanar = cross.norm() <= tol.angular
                    && face_plane.signed_distance(plane.point()).abs() <= tol.length;
                if !coplanar {
                    continue;
                }
                // Outward normal of the face (folding its sense).
                let outward = match face.sense {
                    Sense::Same => face_plane.normal().as_vec(),
                    Sense::Reversed => -face_plane.normal().as_vec(),
                };
                // Convention: include the face only if it agrees with +normal.
                if outward.dot(plane_n) <= tol.angular {
                    continue;
                }
                collect_face_loops(brep, face_id, project, pool);
            }
        }
    }
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

/// Globally nest a flat pool of section loops into holed profiles.
///
/// Each loop's signed area (arc-segment corrected) and an interior representative
/// point drive an exact orient2d containment test. A loop's **parent** is the
/// smallest-area loop that strictly contains it; the *nesting depth* (number of
/// ancestors) decides role by parity:
///
/// * even depth (0, 2, …) → an **outer** boundary → starts a new profile,
/// * odd depth (1, 3, …) → a **hole** of its nearest even-depth ancestor.
///
/// So material-inside-a-hole (depth 2) becomes its own profile, and a hole's
/// hole is correctly a hole of the inner island. Output loops are oriented CCW
/// (outers) / CW (holes) in the frame for a well-formed holed polygon.
fn nest_profiles(loops: Vec<SectionLoop>) -> Vec<SectionProfile> {
    let n = loops.len();
    if n == 0 {
        return Vec::new();
    }
    let areas: Vec<f64> = loops.iter().map(|l| loop_signed_area(l).abs()).collect();
    let reps: Vec<[f64; 2]> = loops.iter().map(representative_point).collect();

    // parent[i] = smallest-area loop strictly containing loop i (or None).
    let mut parent: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let mut best: Option<usize> = None;
        for j in 0..n {
            if i == j {
                continue;
            }
            if areas[j] > areas[i] && point_in_loop_2d(reps[i], &loops[j].points_2d) {
                match best {
                    Some(b) if areas[j] < areas[b] => best = Some(j),
                    None => best = Some(j),
                    _ => {}
                }
            }
        }
        parent[i] = best;
    }

    // Depth (number of ancestors) of each loop.
    let depth = |mut i: usize| -> usize {
        let mut d = 0usize;
        let mut guard = 0usize;
        while let Some(p) = parent[i] {
            d += 1;
            i = p;
            guard += 1;
            if guard > n {
                break;
            }
        }
        d
    };

    // Nearest even-depth ancestor (the outer this hole belongs to).
    let nearest_outer = |i: usize| -> Option<usize> { parent[i] };

    // Build profiles keyed by even-depth outer loop. The index `i` keys
    // `profile_of` / `parent` (via `depth`/`nearest_outer`) *and* selects
    // `loops[i]`, so an index loop is clearer than an enumerate here.
    let mut profile_of: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut profiles: Vec<SectionProfile> = Vec::new();
    #[allow(clippy::needless_range_loop)]
    for i in 0..n {
        if depth(i) % 2 == 0 {
            let mut outer = loops[i].clone();
            orient_loop(&mut outer, true);
            profile_of.insert(i, profiles.len());
            profiles.push(SectionProfile {
                outer,
                holes: Vec::new(),
            });
        }
    }
    #[allow(clippy::needless_range_loop)]
    for i in 0..n {
        if depth(i) % 2 == 1 {
            if let Some(p) = nearest_outer(i) {
                if let Some(&pi) = profile_of.get(&p) {
                    let mut hole = loops[i].clone();
                    orient_loop(&mut hole, false);
                    profiles[pi].holes.push(hole);
                }
            }
        }
    }
    profiles
}

/// Signed area of a section loop in the 2-D frame, with the arc circular-segment
/// correction (same formula as [`crate::boolean::poly2d`] `signed_area`).
fn loop_signed_area(l: &SectionLoop) -> f64 {
    let arcs: Vec<(f64, f64)> = l
        .edges
        .iter()
        .filter_map(|e| match *e {
            SectionEdge::Arc {
                radius,
                start_angle,
                end_angle,
                ..
            } => Some((radius, signed_sweep(start_angle, end_angle))),
            SectionEdge::Line { .. } => None,
        })
        .collect();
    loop_signed_area_2d(&l.points_2d, &arcs)
}

/// The signed sweep `end − start` of an arc, brought into `(−π, π]` (a section
/// arc spans at most a semicircle, so the shorter signed sweep is correct).
fn signed_sweep(start_angle: f64, end_angle: f64) -> f64 {
    let mut d = end_angle - start_angle;
    let two_pi = std::f64::consts::TAU;
    while d <= -std::f64::consts::PI {
        d += two_pi;
    }
    while d > std::f64::consts::PI {
        d -= two_pi;
    }
    d
}

/// Reverse a section loop so its 2-D winding matches `ccw` (outer = CCW, hole =
/// CW). Reverses the vertex rings and the typed edge list, swapping each edge's
/// endpoints / arc sweep so the edge list stays consistent with the ring.
fn orient_loop(l: &mut SectionLoop, ccw: bool) {
    let positive = loop_signed_area(l) > 0.0;
    if positive == ccw {
        return;
    }
    l.points_3d.reverse();
    l.points_2d.reverse();
    let mut rev: Vec<SectionEdge> = l
        .edges
        .iter()
        .rev()
        .map(|e| match *e {
            SectionEdge::Line { start, end } => SectionEdge::Line {
                start: end,
                end: start,
            },
            SectionEdge::Arc {
                center,
                radius,
                start_angle,
                end_angle,
                start,
                end,
            } => SectionEdge::Arc {
                center,
                radius,
                start_angle: end_angle,
                end_angle: start_angle,
                start: end,
                end: start,
            },
        })
        .collect();
    std::mem::swap(&mut l.edges, &mut rev);
}

/// An interior representative point of a section loop in the 2-D frame.
///
/// For a polygonal loop the centroid of the first triangle of the ring is
/// interior; for a two-vertex arc loop (a round member) the centre of the two
/// seam points plus a nudge toward the arc bulge is used.
fn representative_point(l: &SectionLoop) -> [f64; 2] {
    let p = &l.points_2d;
    if p.len() >= 3 {
        return [
            (p[0][0] + p[1][0] + p[2][0]) / 3.0,
            (p[0][1] + p[1][1] + p[2][1]) / 3.0,
        ];
    }
    // Arc loop: use the first arc's centre, which is interior to a full disk /
    // annulus section.
    for e in &l.edges {
        if let SectionEdge::Arc { center, .. } = e {
            return *center;
        }
    }
    p.first().copied().unwrap_or([0.0, 0.0])
}
