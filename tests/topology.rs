//! Hand-built topology fixtures exercising `validate()`.
//!
//! The fixtures are assembled with a small in-test builder ([`SolidBuilder`])
//! that, given each face as an ordered ring of explicit corner coordinates,
//! produces shared vertices, shared curves and correctly reversed sibling
//! boundaries. Every numeric literal carries an `f64` annotation and explicit
//! tolerances, per `DESIGN.md` §12.

use std::collections::HashMap;

use archi_kernel::brep::Brep;
use archi_kernel::csg::{CsgNode, EvalError, Member, Profile2d};
use archi_kernel::geom::{CurveGeom, SurfaceGeom, VertexGeom};
use archi_kernel::topo::arena::Id;
use archi_kernel::topo::validate::{validate_topology, Defect};
use archi_kernel::topo::{Face, HalfEdge, Loop, Sense, Shell, Solid, Vertex};
use archi_kernel::{Line3, Plane, Point3, Tol, ValidateLevel, Vec3};

const EPS: f64 = 1e-9;

// ── Builder ──────────────────────────────────────────────────────────────────

/// Integer-keyed coordinate so that "the same corner" deduplicates exactly,
/// despite `f64` not being `Hash`/`Eq`. Fixtures use coordinates on a coarse
/// grid, so scaling by 1e6 and rounding is exact here.
type CoordKey = (i64, i64, i64);

fn key(p: Point3) -> CoordKey {
    let q = |v: f64| (v * 1.0e6_f64).round() as i64;
    (q(p.x), q(p.y), q(p.z))
}

/// Builds a single solid (one shell) face by face, sharing vertices and curves.
struct SolidBuilder {
    brep: Brep,
    vertices: HashMap<CoordKey, Id<Vertex>>,
    coords: HashMap<CoordKey, Point3>,
    /// Edge (unordered vertex pair) → (curve, origin coord key). The origin is
    /// the endpoint chosen as the line's parameter origin, so both half-edges
    /// agree on the parameterisation and siblings get reversed boundaries.
    edges: HashMap<(CoordKey, CoordKey), (archi_kernel::CurveId, CoordKey)>,
    faces: Vec<Id<Face>>,
    tol: Tol,
}

impl SolidBuilder {
    fn new() -> Self {
        Self {
            brep: Brep::new(),
            vertices: HashMap::new(),
            coords: HashMap::new(),
            edges: HashMap::new(),
            faces: Vec::new(),
            tol: Tol::default(),
        }
    }

    fn vertex(&mut self, p: Point3) -> Id<Vertex> {
        let k = key(p);
        if let Some(&v) = self.vertices.get(&k) {
            return v;
        }
        let pid = self.brep.geom.insert_point(VertexGeom::Explicit(p));
        let v = self.brep.topo.add_vertex(Vertex { point: pid });
        self.vertices.insert(k, v);
        self.coords.insert(k, p);
        v
    }

    /// Get or create the shared line for an edge between two corners, returning
    /// the curve and the origin coord key for its parameterisation.
    fn edge_curve(&mut self, a: Point3, b: Point3) -> (archi_kernel::CurveId, CoordKey) {
        let (ka, kb) = (key(a), key(b));
        let unordered = if ka <= kb { (ka, kb) } else { (kb, ka) };
        if let Some(&entry) = self.edges.get(&unordered) {
            return entry;
        }
        // Parameter origin is the lexicographically smaller endpoint.
        let (origin_key, origin_pt, other_pt) = if ka <= kb { (ka, a, b) } else { (kb, b, a) };
        let dir = other_pt - origin_pt;
        let line = Line3::new(origin_pt, dir).expect("non-degenerate edge");
        let cid = self.brep.geom.insert_curve(CurveGeom::Line(line));
        let entry = (cid, origin_key);
        self.edges.insert(unordered, entry);
        entry
    }

    /// Add a face from an ordered ring of corners (CCW seen from outside) on the
    /// given plane. Returns the face id.
    fn face(&mut self, ring: &[Point3], plane: Plane) -> Id<Face> {
        let (surface, flipped) = self.brep.geom.insert_plane(plane, &self.tol);
        let sense = if flipped {
            Sense::Reversed
        } else {
            Sense::Same
        };
        let outer = self.loop_from_ring(ring);
        self.brep.topo.add_face(Face {
            surface,
            sense,
            outer,
            inners: Vec::new(),
        })
    }

    /// Add a face with interior (hole) loops.
    fn face_with_holes(&mut self, ring: &[Point3], holes: &[&[Point3]], plane: Plane) -> Id<Face> {
        let (surface, flipped) = self.brep.geom.insert_plane(plane, &self.tol);
        let sense = if flipped {
            Sense::Reversed
        } else {
            Sense::Same
        };
        let outer = self.loop_from_ring(ring);
        let inners: Vec<Id<Loop>> = holes.iter().map(|h| self.loop_from_ring(h)).collect();
        self.brep.topo.add_face(Face {
            surface,
            sense,
            outer,
            inners,
        })
    }

    /// Build a loop from an ordered ring; half-edges share vertices/curves and
    /// boundaries are oriented along the shared line, so the sibling on the
    /// adjacent face gets the reversed boundary automatically.
    fn loop_from_ring(&mut self, ring: &[Point3]) -> Id<Loop> {
        let n = ring.len();
        let mut hes = Vec::with_capacity(n);
        for i in 0..n {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            let start = self.vertex(a);
            let _ = self.vertex(b);
            let (curve, origin_key) = self.edge_curve(a, b);
            // Parameter = signed arc length from the curve origin. Since dir is
            // (other - origin) un-normalised? No: Line3 stores a *unit* dir, so
            // the parameter is true arc length. Compute both endpoints' params.
            let origin_pt = self.coords[&origin_key];
            let line = match self.brep.geom.curve(curve).expect("curve") {
                CurveGeom::Line(l) => *l,
                _ => unreachable!("edges are lines"),
            };
            let ta = (a - origin_pt).dot(line.dir().as_vec());
            let tb = (b - origin_pt).dot(line.dir().as_vec());
            let he = self.brep.topo.add_half_edge(HalfEdge {
                start,
                curve,
                boundary: [ta, tb],
            });
            hes.push(he);
        }
        self.brep.topo.add_loop(Loop { half_edges: hes })
    }

    /// Finish: wrap the accumulated faces in a shell and a solid.
    fn finish(mut self) -> (Brep, Id<Solid>) {
        let shell = self.brep.topo.add_shell(Shell {
            faces: self.faces.clone(),
        });
        let solid = self.brep.topo.add_solid(Solid {
            shells: vec![shell],
        });
        self.brep.solids = vec![solid];
        (self.brep, solid)
    }

    fn push_face(&mut self, f: Id<Face>) {
        self.faces.push(f);
    }
}

/// Convenience: a plane through `p` with outward normal `n`.
fn plane(p: Point3, n: Vec3) -> Plane {
    Plane::new(p, n).expect("valid plane")
}

// ── Cube ─────────────────────────────────────────────────────────────────────

/// Build the unit cube `[0,1]³` with outward-facing loops.
fn build_cube() -> (Brep, Id<Solid>) {
    let mut b = SolidBuilder::new();
    // 8 corners.
    let v = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let p000 = v(0.0, 0.0, 0.0);
    let p100 = v(1.0, 0.0, 0.0);
    let p110 = v(1.0, 1.0, 0.0);
    let p010 = v(0.0, 1.0, 0.0);
    let p001 = v(0.0, 0.0, 1.0);
    let p101 = v(1.0, 0.0, 1.0);
    let p111 = v(1.0, 1.0, 1.0);
    let p011 = v(0.0, 1.0, 1.0);

    // bottom z=0, outward normal -z, CCW seen from below.
    let f = b.face(
        &[p000, p010, p110, p100],
        plane(p000, Vec3::new(0.0, 0.0, -1.0)),
    );
    b.push_face(f);
    // top z=1, outward +z, CCW seen from above.
    let f = b.face(
        &[p001, p101, p111, p011],
        plane(p001, Vec3::new(0.0, 0.0, 1.0)),
    );
    b.push_face(f);
    // front y=0, outward -y.
    let f = b.face(
        &[p000, p100, p101, p001],
        plane(p000, Vec3::new(0.0, -1.0, 0.0)),
    );
    b.push_face(f);
    // back y=1, outward +y.
    let f = b.face(
        &[p010, p011, p111, p110],
        plane(p010, Vec3::new(0.0, 1.0, 0.0)),
    );
    b.push_face(f);
    // left x=0, outward -x.
    let f = b.face(
        &[p000, p001, p011, p010],
        plane(p000, Vec3::new(-1.0, 0.0, 0.0)),
    );
    b.push_face(f);
    // right x=1, outward +x.
    let f = b.face(
        &[p100, p110, p111, p101],
        plane(p100, Vec3::new(1.0, 0.0, 0.0)),
    );
    b.push_face(f);

    b.finish()
}

#[test]
fn cube_validates_light_and_full() {
    let tol = Tol::default();
    let (brep, solid) = build_cube();

    // Counts: 8V, 24HE (12 edges × 2), 6F, 6L, 1S, genus 0.
    assert_eq!(brep.topo.vertices.len(), 8usize);
    assert_eq!(brep.topo.half_edges.len(), 24usize);
    assert_eq!(brep.topo.faces.len(), 6usize);
    assert_eq!(brep.topo.loops.len(), 6usize);

    validate_topology(&brep.topo, &[solid], &tol, Some(0u32))
        .expect("cube must validate with genus 0");
    brep.validate(&tol, ValidateLevel::Light)
        .expect("cube light");
    brep.validate(&tol, ValidateLevel::Full).expect("cube full");
}

// ── L-shaped solid (concave, genus 0) ───────────────────────────────────────

/// An extruded L-shaped prism: footprint is an L in the xy-plane, extruded in z.
/// Footprint corners (CCW): (0,0)→(2,0)→(2,1)→(1,1)→(1,2)→(0,2).
fn build_l_solid() -> (Brep, Id<Solid>) {
    let mut b = SolidBuilder::new();
    let foot = [
        (0.0_f64, 0.0_f64),
        (2.0_f64, 0.0_f64),
        (2.0_f64, 1.0_f64),
        (1.0_f64, 1.0_f64),
        (1.0_f64, 2.0_f64),
        (0.0_f64, 2.0_f64),
    ];
    let lo = |x: f64, y: f64| Point3::new(x, y, 0.0_f64);
    let hi = |x: f64, y: f64| Point3::new(x, y, 1.0_f64);

    // bottom (z=0), outward -z: CCW seen from below = reverse the footprint.
    let bottom: Vec<Point3> = foot.iter().rev().map(|&(x, y)| lo(x, y)).collect();
    let f = b.face(&bottom, plane(lo(0.0, 0.0), Vec3::new(0.0, 0.0, -1.0)));
    b.push_face(f);

    // top (z=1), outward +z: CCW seen from above = footprint order.
    let top: Vec<Point3> = foot.iter().map(|&(x, y)| hi(x, y)).collect();
    let f = b.face(&top, plane(hi(0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)));
    b.push_face(f);

    // side walls, one per footprint edge.
    let n = foot.len();
    for i in 0..n {
        let (x0, y0) = foot[i];
        let (x1, y1) = foot[(i + 1) % n];
        // outward normal = edge direction rotated -90° in xy (right-hand outside
        // for a CCW footprint).
        let ex = x1 - x0;
        let ey = y1 - y0;
        let nrm = Vec3::new(ey, -ex, 0.0_f64);
        // ring CCW as seen from outside: lo0 → lo1 → hi1 → hi0.
        let ring = [lo(x0, y0), lo(x1, y1), hi(x1, y1), hi(x0, y0)];
        let f = b.face(&ring, plane(lo(x0, y0), nrm));
        b.push_face(f);
    }

    b.finish()
}

#[test]
fn l_solid_validates() {
    let tol = Tol::default();
    let (brep, solid) = build_l_solid();
    // 6 footprint corners → 12 vertices, 8 faces (top+bottom+6 walls), genus 0.
    assert_eq!(brep.topo.vertices.len(), 12usize);
    assert_eq!(brep.topo.faces.len(), 8usize);
    validate_topology(&brep.topo, &[solid], &tol, Some(0u32)).expect("L solid genus 0");
    brep.validate(&tol, ValidateLevel::Full).expect("L full");
}

// ── Box with a through hole (torus topology, genus 1) ────────────────────────

/// A `[0,3]×[0,3]×[0,1]` box pierced by a `[1,2]×[1,2]` square hole through z.
/// Top and bottom faces each carry one inner loop (the hole rim); the hole adds
/// 4 inner side walls. Genus 1.
fn build_box_with_hole() -> (Brep, Id<Solid>) {
    let mut b = SolidBuilder::new();
    let lo = |x: f64, y: f64| Point3::new(x, y, 0.0_f64);
    let hi = |x: f64, y: f64| Point3::new(x, y, 1.0_f64);

    // Outer footprint (CCW): (0,0)(3,0)(3,3)(0,3). Hole footprint (1,1)(2,1)(2,2)(1,2).
    let outer = [
        (0.0_f64, 0.0_f64),
        (3.0_f64, 0.0_f64),
        (3.0_f64, 3.0_f64),
        (0.0_f64, 3.0_f64),
    ];
    let hole = [
        (1.0_f64, 1.0_f64),
        (2.0_f64, 1.0_f64),
        (2.0_f64, 2.0_f64),
        (1.0_f64, 2.0_f64),
    ];

    // Bottom z=0, outward -z. Outer ring CW-from-below; hole ring opposite.
    let bottom_outer: Vec<Point3> = outer.iter().rev().map(|&(x, y)| lo(x, y)).collect();
    let bottom_hole: Vec<Point3> = hole.iter().map(|&(x, y)| lo(x, y)).collect();
    let f = b.face_with_holes(
        &bottom_outer,
        &[&bottom_hole],
        plane(lo(0.0, 0.0), Vec3::new(0.0, 0.0, -1.0)),
    );
    b.push_face(f);

    // Top z=1, outward +z. Outer ring CCW-from-above; hole ring opposite.
    let top_outer: Vec<Point3> = outer.iter().map(|&(x, y)| hi(x, y)).collect();
    let top_hole: Vec<Point3> = hole.iter().rev().map(|&(x, y)| hi(x, y)).collect();
    let f = b.face_with_holes(
        &top_outer,
        &[&top_hole],
        plane(hi(0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
    );
    b.push_face(f);

    // Outer side walls (4), outward.
    for i in 0..4 {
        let (x0, y0) = outer[i];
        let (x1, y1) = outer[(i + 1) % 4];
        let nrm = Vec3::new(y1 - y0, -(x1 - x0), 0.0_f64);
        let ring = [lo(x0, y0), lo(x1, y1), hi(x1, y1), hi(x0, y0)];
        let f = b.face(&ring, plane(lo(x0, y0), nrm));
        b.push_face(f);
    }

    // Inner hole walls (4), normal points *into* the hole (toward the void),
    // i.e. opposite the outer-wall convention. The hole footprint is CCW, but
    // the inner surface faces inward, so we traverse the wall the other way.
    for i in 0..4 {
        let (x0, y0) = hole[i];
        let (x1, y1) = hole[(i + 1) % 4];
        // Inward-facing normal for a CCW hole footprint.
        let nrm = Vec3::new(-(y1 - y0), x1 - x0, 0.0_f64);
        // Reverse traversal so the wall's outward (into-void) normal is consistent.
        let ring = [lo(x1, y1), lo(x0, y0), hi(x0, y0), hi(x1, y1)];
        let f = b.face(&ring, plane(lo(x0, y0), nrm));
        b.push_face(f);
    }

    b.finish()
}

#[test]
fn box_with_hole_has_genus_one() {
    let tol = Tol::default();
    let (brep, solid) = build_box_with_hole();
    // V: 8 outer + 8 hole = 16. F: 2 caps + 4 outer walls + 4 inner walls = 10.
    // L: 2 caps each with outer+inner = 4 loops, + 8 wall loops = 12.
    assert_eq!(brep.topo.vertices.len(), 16usize);
    assert_eq!(brep.topo.faces.len(), 10usize);
    assert_eq!(brep.topo.loops.len(), 12usize);

    // The ring-term Euler must derive genus 1.
    validate_topology(&brep.topo, &[solid], &tol, Some(1u32))
        .expect("box with hole must be genus 1");

    // The ring-free formula would mis-derive the genus; assert genus 0 is
    // rejected (this is exactly why the ring term exists).
    let err = validate_topology(&brep.topo, &[solid], &tol, Some(0u32))
        .expect_err("genus 0 must be rejected for a torus");
    assert!(err
        .iter()
        .any(|d| matches!(d, Defect::GenusMismatch { .. })));

    brep.validate(&tol, ValidateLevel::Full)
        .expect("box with hole full");
}

// ── Defect detection ─────────────────────────────────────────────────────────

#[test]
fn missing_face_breaks_watertightness() {
    let tol = Tol::default();
    let (mut brep, solid) = build_cube();
    // Drop one face from the shell.
    let shell_id = brep.topo.solids.get(solid).unwrap().shells[0];
    let dropped = brep.topo.shells.get(shell_id).unwrap().faces[0];
    brep.topo
        .shells
        .get_mut(shell_id)
        .unwrap()
        .faces
        .retain(|&f| f != dropped);

    let err = validate_topology(&brep.topo, &[solid], &tol, None)
        .expect_err("removing a face must be detected");
    // The four edges of the dropped face now lack siblings.
    assert!(
        err.iter()
            .any(|d| matches!(d, Defect::MissingSibling { .. })),
        "expected MissingSibling, got {err:?}"
    );
}

#[test]
fn located_defect_attaches_a_coordinate() {
    // The B-rep layer augments a topological defect with a representative world
    // coordinate (`docs/design/progress.md`: Defect への座標付加). A missing
    // sibling must come back with the offending half-edge's start position, and
    // that position must be a real corner of the unit cube.
    let tol = Tol::default();
    let (mut brep, solid) = build_cube();
    let shell_id = brep.topo.solids.get(solid).unwrap().shells[0];
    let dropped = brep.topo.shells.get(shell_id).unwrap().faces[0];
    brep.topo
        .shells
        .get_mut(shell_id)
        .unwrap()
        .faces
        .retain(|&f| f != dropped);
    brep.solids = vec![solid];

    let located = brep
        .validate_located(&tol, ValidateLevel::Full)
        .expect_err("missing sibling must be located");
    let with_coord = located
        .iter()
        .find(|d| matches!(d.defect, Defect::MissingSibling { .. }))
        .expect("a MissingSibling located defect");
    let p = with_coord
        .location
        .expect("MissingSibling carries a representative coordinate");
    // Every coordinate of a unit-cube corner is 0 or 1.
    for c in p {
        assert!(
            c.abs() <= tol.length || (c - 1.0).abs() <= tol.length,
            "located coordinate {c} is not a unit-cube corner"
        );
    }
}

#[test]
fn unreversed_boundary_breaks_sibling_pairing() {
    let tol = Tol::default();
    let (mut brep, solid) = build_cube();
    // Pick any half-edge and corrupt its boundary so its sibling no longer
    // matches (forget to reverse: make boundary equal on both endpoints).
    let some_he = brep.topo.half_edges.ids().next().expect("has half-edges");
    let b = brep.topo.half_edges.get(some_he).unwrap().boundary;
    // Set boundary to a non-reversed value vs its sibling: shift it.
    brep.topo.half_edges.get_mut(some_he).unwrap().boundary = [b[0] + 10.0_f64, b[1] + 10.0_f64];

    let err = validate_topology(&brep.topo, &[solid], &tol, None)
        .expect_err("corrupted boundary must be detected");
    assert!(
        err.iter()
            .any(|d| matches!(d, Defect::MissingSibling { .. })),
        "expected MissingSibling from broken reversal, got {err:?}"
    );
}

#[test]
fn dangling_handle_is_detected() {
    let tol = Tol::default();
    let (mut brep, solid) = build_cube();
    // Remove a vertex that a half-edge still references.
    let some_he = brep.topo.half_edges.ids().next().expect("has half-edges");
    let v = brep.topo.half_edges.get(some_he).unwrap().start;
    brep.topo.vertices.remove(v);

    let err = validate_topology(&brep.topo, &[solid], &tol, None)
        .expect_err("dangling vertex handle must be detected");
    assert!(
        err.iter()
            .any(|d| matches!(d, Defect::DanglingReference { .. })),
        "expected DanglingReference, got {err:?}"
    );
}

// ── Plane canonicalisation ───────────────────────────────────────────────────

#[test]
fn same_plane_twice_yields_same_id() {
    let mut brep = Brep::new();
    let tol = Tol::default();
    let p = plane(Point3::origin(), Vec3::Z);
    let (id1, f1) = brep.geom.insert_plane(p, &tol);
    let (id2, f2) = brep.geom.insert_plane(p, &tol);
    assert_eq!(id1, id2);
    assert!(!f1);
    assert!(!f2);
    assert_eq!(brep.geom.surface_count(), 1usize);
}

#[test]
fn opposite_normal_is_canonicalised_and_flagged() {
    let mut brep = Brep::new();
    let tol = Tol::default();
    let up = plane(Point3::origin(), Vec3::Z);
    let down = plane(Point3::origin(), Vec3::new(0.0_f64, 0.0_f64, -1.0_f64));
    let (id_up, f_up) = brep.geom.insert_plane(up, &tol);
    let (id_down, f_down) = brep.geom.insert_plane(down, &tol);
    assert_eq!(id_up, id_down, "opposite normals are the same plane");
    assert!(!f_up);
    assert!(f_down, "the down-facing normal must be flagged flipped");
    assert_eq!(brep.geom.surface_count(), 1usize);
}

#[test]
fn distinct_planes_get_distinct_ids() {
    let mut brep = Brep::new();
    let tol = Tol::default();
    let z0 = plane(Point3::origin(), Vec3::Z);
    let z1 = plane(Point3::new(0.0_f64, 0.0_f64, 1.0_f64), Vec3::Z);
    let (id0, _) = brep.geom.insert_plane(z0, &tol);
    let (id1, _) = brep.geom.insert_plane(z1, &tol);
    assert_ne!(id0, id1, "planes 1 m apart are distinct (beyond Tol)");
    assert_eq!(brep.geom.surface_count(), 2usize);
}

// ── CSG container smoke test ─────────────────────────────────────────────────

#[test]
fn member_extrude_evaluates_and_non_extrude_is_not_yet_implemented() {
    let tol = Tol::default();
    // An Extrude leaf now evaluates for real (Phase 2).
    let profile = Profile2d::rect(0.15_f64, 0.3_f64).expect("valid rect");
    let node = CsgNode::Extrude {
        origin: Point3::origin(),
        profile,
        axis: Vec3::Z,
        length: 3.0_f64,
    };
    let mut member = Member::new(node);
    member.brep(&tol).expect("extrude member must evaluate");
    assert!(member.last_valid().is_some());

    // The priority-clip node is still unimplemented in this phase (booleans and
    // openings are wired in Phase 3b, but Clip is not).
    let mut other = Member::new(CsgNode::Clip {
        base: Box::new(CsgNode::Union(Vec::new())),
        clippers: Vec::new(),
        rule: archi_kernel::csg::ClipRule::Priority,
    });
    assert!(matches!(
        other.brep(&tol),
        Err(EvalError::NotYetImplemented)
    ));
    assert!(other.last_valid().is_none());
}

// Keep the surface module reachable from tests without an unused import warning.
#[test]
fn surface_geom_signed_distance_is_exposed() {
    let s = SurfaceGeom::Plane(plane(Point3::origin(), Vec3::Z));
    assert!((s.signed_distance(Point3::new(0.0_f64, 0.0_f64, 2.0_f64)) - 2.0_f64).abs() < EPS);
}
