//! Phase 3a tests: solid × half-space cut + section extraction.
//!
//! Volume is checked with `brep.signed_volume()`. The volume tolerance is
//! derived from `Tol::length` (1e-6 m) and the dimensions involved: a cut at
//! building scale (≤ a few metres) accumulates a few dozen face-integral
//! round-offs, all far below 1e-6 m × area, so `VOL_EPS = 1e-9 m³` is far above
//! the f64 round-off yet far tighter than any real geometric error (the same
//! bound the Phase 2 extrusion tests use).

use archi_kernel::boolean::{cut, CutResult, KeepSide};
use archi_kernel::brep::Brep;
use archi_kernel::build::extrude;
use archi_kernel::csg::Profile2d;
use archi_kernel::geom::{CurveGeom, VertexGeom};
use archi_kernel::section::{section, SectionEdge};
use archi_kernel::topo::arena::Id;
use archi_kernel::topo::{Face, HalfEdge, Loop, Sense, Shell, Solid, Vertex};
use archi_kernel::{Line3, Plane, Point3, Tol, ValidateLevel, Vec3};
use std::collections::HashMap;
use std::f64::consts::PI;

const VOL_EPS: f64 = 1e-9;

/// A larger tolerance for cuts that introduce ellipse / cylinder curves, whose
/// closed-form parameters accumulate more f64 round-off than axis-aligned cuts.
const VOL_EPS_CURVED: f64 = 1e-6;

// ── Builder (same shape as tests/topology.rs) ───────────────────────────────

type CoordKey = (i64, i64, i64);

fn key(p: Point3) -> CoordKey {
    let q = |v: f64| (v * 1.0e6_f64).round() as i64;
    (q(p.x), q(p.y), q(p.z))
}

struct SolidBuilder {
    brep: Brep,
    vertices: HashMap<CoordKey, Id<Vertex>>,
    coords: HashMap<CoordKey, Point3>,
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

    fn edge_curve(&mut self, a: Point3, b: Point3) -> (archi_kernel::CurveId, CoordKey) {
        let (ka, kb) = (key(a), key(b));
        let unordered = if ka <= kb { (ka, kb) } else { (kb, ka) };
        if let Some(&entry) = self.edges.get(&unordered) {
            return entry;
        }
        let (origin_key, origin_pt, other_pt) = if ka <= kb { (ka, a, b) } else { (kb, b, a) };
        let dir = other_pt - origin_pt;
        let line = Line3::new(origin_pt, dir).expect("non-degenerate edge");
        let cid = self.brep.geom.insert_curve(CurveGeom::Line(line));
        let entry = (cid, origin_key);
        self.edges.insert(unordered, entry);
        entry
    }

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

    fn loop_from_ring(&mut self, ring: &[Point3]) -> Id<Loop> {
        let n = ring.len();
        let mut hes = Vec::with_capacity(n);
        for i in 0..n {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            let start = self.vertex(a);
            let _ = self.vertex(b);
            let (curve, origin_key) = self.edge_curve(a, b);
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

fn plane(p: Point3, n: Vec3) -> Plane {
    Plane::new(p, n).expect("valid plane")
}

/// Build an axis-aligned box `[0,sx]×[0,sy]×[0,sz]` with outward loops.
fn build_box(sx: f64, sy: f64, sz: f64) -> (Brep, Id<Solid>) {
    let mut b = SolidBuilder::new();
    let v = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let p000 = v(0.0, 0.0, 0.0);
    let p100 = v(sx, 0.0, 0.0);
    let p110 = v(sx, sy, 0.0);
    let p010 = v(0.0, sy, 0.0);
    let p001 = v(0.0, 0.0, sz);
    let p101 = v(sx, 0.0, sz);
    let p111 = v(sx, sy, sz);
    let p011 = v(0.0, sy, sz);

    let f = b.face(
        &[p000, p010, p110, p100],
        plane(p000, Vec3::new(0.0, 0.0, -1.0)),
    );
    b.push_face(f);
    let f = b.face(
        &[p001, p101, p111, p011],
        plane(p001, Vec3::new(0.0, 0.0, 1.0)),
    );
    b.push_face(f);
    let f = b.face(
        &[p000, p100, p101, p001],
        plane(p000, Vec3::new(0.0, -1.0, 0.0)),
    );
    b.push_face(f);
    let f = b.face(
        &[p010, p011, p111, p110],
        plane(p010, Vec3::new(0.0, 1.0, 0.0)),
    );
    b.push_face(f);
    let f = b.face(
        &[p000, p001, p011, p010],
        plane(p000, Vec3::new(-1.0, 0.0, 0.0)),
    );
    b.push_face(f);
    let f = b.face(
        &[p100, p110, p111, p101],
        plane(p100, Vec3::new(1.0, 0.0, 0.0)),
    );
    b.push_face(f);
    b.finish()
}

// ── Cube cut in the middle ──────────────────────────────────────────────────

#[test]
fn cube_horizontal_midcut() {
    let tol = Tol::default();
    let (brep, solid) = build_box(1.0, 1.0, 1.0);
    let total = brep.signed_volume();
    assert!((total - 1.0_f64).abs() < VOL_EPS);

    // Cut plane z = 0.5, normal +z. Below keeps z ≤ 0.5.
    let cut_plane = plane(Point3::new(0.0, 0.0, 0.5), Vec3::Z);

    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("cut below");
    let CutResult::Cut { brep: bb, caps } = below else {
        panic!("expected a real cut, got {below:?}");
    };
    bb.validate(&tol, ValidateLevel::Full).expect("below valid");
    assert_eq!(caps.len(), 1usize, "one cap");
    let v_below = bb.signed_volume();
    assert!(
        (v_below - 0.5_f64).abs() < VOL_EPS,
        "below volume {v_below}, expected 0.5"
    );

    let above = cut(&brep, solid, &cut_plane, KeepSide::Above, &tol).expect("cut above");
    let above_brep = above.brep();
    above_brep
        .validate(&tol, ValidateLevel::Full)
        .expect("above valid");
    let v_above = above_brep.signed_volume();
    assert!(
        (v_above - 0.5_f64).abs() < VOL_EPS,
        "above volume {v_above}"
    );

    assert!(
        (v_below + v_above - total).abs() < VOL_EPS,
        "halves must sum to the whole: {v_below} + {v_above} vs {total}"
    );
}

// ── Coplanar cut (cut plane lands on an existing face) ─────────────────────

#[test]
fn cube_cut_on_existing_top_face() {
    let tol = Tol::default();
    let (brep, solid) = build_box(1.0, 1.0, 1.0);

    // Cut at z = 1 (the top face), normal +z. Below keeps z ≤ 1 = the whole
    // cube; the top face (outward +z = +normal) is the cap lid and is kept.
    let cut_plane = plane(Point3::new(0.0, 0.0, 1.0), Vec3::Z);

    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    let bb = below.brep();
    bb.validate(&tol, ValidateLevel::Full).expect("below valid");
    let v = bb.signed_volume();
    assert!(
        (v - 1.0_f64).abs() < VOL_EPS,
        "whole cube kept, V = {v}, expected 1.0"
    );

    // Above keeps z ≥ 1: only the top face lies on the plane, nothing strictly
    // above ⇒ the result is empty.
    let above = cut(&brep, solid, &cut_plane, KeepSide::Above, &tol).expect("above");
    assert!(
        matches!(above, CutResult::Empty),
        "nothing above z=1, got {above:?}"
    );
}

#[test]
fn cube_cut_on_existing_bottom_face() {
    let tol = Tol::default();
    let (brep, solid) = build_box(1.0, 1.0, 1.0);

    // Cut at z = 0 (the bottom face), normal +z.
    let cut_plane = plane(Point3::new(0.0, 0.0, 0.0), Vec3::Z);

    // Above keeps z ≥ 0 = whole cube; the bottom face (outward −z = −normal) is
    // the lid and is kept.
    let above = cut(&brep, solid, &cut_plane, KeepSide::Above, &tol).expect("above");
    let ab = above.brep();
    ab.validate(&tol, ValidateLevel::Full).expect("above valid");
    let v = ab.signed_volume();
    assert!((v - 1.0_f64).abs() < VOL_EPS, "whole cube kept, V = {v}");

    // Below keeps z ≤ 0 ⇒ empty.
    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    assert!(matches!(below, CutResult::Empty), "nothing below z=0");
}

// ── Plane missing the solid ─────────────────────────────────────────────────

#[test]
fn plane_misses_solid_all_kept_or_empty() {
    let tol = Tol::default();
    let (brep, solid) = build_box(1.0, 1.0, 1.0);

    // Plane well above the cube (z = 5).
    let cut_plane = plane(Point3::new(0.0, 0.0, 5.0), Vec3::Z);

    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    assert!(matches!(below, CutResult::AllKept { .. }), "all kept");
    let v = below.brep().signed_volume();
    assert!((v - 1.0_f64).abs() < VOL_EPS);

    let above = cut(&brep, solid, &cut_plane, KeepSide::Above, &tol).expect("above");
    assert!(matches!(above, CutResult::Empty), "nothing above");
}

// ── Oblique cut slicing a corner off ────────────────────────────────────────

#[test]
fn cube_corner_cut_oblique() {
    let tol = Tol::default();
    let (brep, solid) = build_box(1.0, 1.0, 1.0);

    // Plane x + y + z = 0.5 through the cube, normal (1,1,1)/√3. Below keeps the
    // x+y+z ≤ 0.5 corner near the origin: a corner tetrahedron of leg 0.5.
    let n = Vec3::new(1.0, 1.0, 1.0);
    // Point on the plane: (0.5, 0, 0).
    let cut_plane = plane(Point3::new(0.5, 0.0, 0.0), n);

    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    let bb = below.brep();
    bb.validate(&tol, ValidateLevel::Full)
        .expect("corner cut valid");

    // Corner tetra with legs a = 0.5 (cut at x+y+z = 0.5) ⇒ V = a³ / 6.
    let a = 0.5_f64;
    let expected = a * a * a / 6.0_f64;
    let v = bb.signed_volume();
    assert!(
        (v - expected).abs() < VOL_EPS,
        "corner tetra V = {v}, expected {expected}"
    );
}

// ── Prism from an arbitrary footprint (concave H) ───────────────────────────

/// Build a prism over a CCW footprint, extruded in +z from z=0 to z=h.
fn build_prism(foot: &[(f64, f64)], h: f64) -> (Brep, Id<Solid>) {
    let mut b = SolidBuilder::new();
    let lo = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let hi = |x: f64, y: f64| Point3::new(x, y, h);

    // bottom (z=0), outward -z: reverse footprint.
    let bottom: Vec<Point3> = foot.iter().rev().map(|&(x, y)| lo(x, y)).collect();
    let f = b.face(
        &bottom,
        plane(lo(foot[0].0, foot[0].1), Vec3::new(0.0, 0.0, -1.0)),
    );
    b.push_face(f);
    // top (z=h), outward +z: footprint order.
    let top: Vec<Point3> = foot.iter().map(|&(x, y)| hi(x, y)).collect();
    let f = b.face(
        &top,
        plane(hi(foot[0].0, foot[0].1), Vec3::new(0.0, 0.0, 1.0)),
    );
    b.push_face(f);
    // side walls.
    let n = foot.len();
    for i in 0..n {
        let (x0, y0) = foot[i];
        let (x1, y1) = foot[(i + 1) % n];
        let nrm = Vec3::new(y1 - y0, -(x1 - x0), 0.0);
        let ring = [lo(x0, y0), lo(x1, y1), hi(x1, y1), hi(x0, y0)];
        let f = b.face(&ring, plane(lo(x0, y0), nrm));
        b.push_face(f);
    }
    b.finish()
}

/// Signed area of a CCW footprint.
fn footprint_area(foot: &[(f64, f64)]) -> f64 {
    let n = foot.len();
    let mut a = 0.0_f64;
    for i in 0..n {
        let (x0, y0) = foot[i];
        let (x1, y1) = foot[(i + 1) % n];
        a += x0 * y1 - x1 * y0;
    }
    a / 2.0_f64
}

#[test]
fn h_prism_midcut_concave_section() {
    let tol = Tol::default();
    // An H footprint (concave): width 0.3 flanges, central web. Use simple
    // explicit corners for an I-beam-like H lying in xy, extruded in z.
    // Overall 0.6 wide (x: 0..0.6), 0.4 tall (y: 0..0.4), flange thick 0.1,
    // web from x=0.25..0.35.
    let foot = [
        (0.0_f64, 0.0_f64),
        (0.6_f64, 0.0_f64),
        (0.6_f64, 0.1_f64),
        (0.35_f64, 0.1_f64),
        (0.35_f64, 0.3_f64),
        (0.6_f64, 0.3_f64),
        (0.6_f64, 0.4_f64),
        (0.0_f64, 0.4_f64),
        (0.0_f64, 0.3_f64),
        (0.25_f64, 0.3_f64),
        (0.25_f64, 0.1_f64),
        (0.0_f64, 0.1_f64),
    ];
    let h = 2.0_f64;
    let (brep, solid) = build_prism(&foot, h);
    let area = footprint_area(&foot);
    let total = brep.signed_volume();
    assert!((total - area * h).abs() < VOL_EPS, "prism V {total}");

    // Cut horizontally at z = 0.8 (mid-height). The cap is the full H section.
    let cut_plane = plane(Point3::new(0.0, 0.0, 0.8), Vec3::Z);
    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    let CutResult::Cut { brep: bb, caps } = below else {
        panic!("expected cut");
    };
    bb.validate(&tol, ValidateLevel::Full).expect("H cut valid");
    assert_eq!(caps.len(), 1usize, "one H-shaped cap");
    // The cap loop should be the concave H outline (12 vertices).
    let cap = bb.topo.faces.get(caps[0]).unwrap();
    let cap_loop = bb.topo.loops.get(cap.outer).unwrap();
    assert_eq!(
        cap_loop.half_edges.len(),
        12usize,
        "H cap is a 12-gon, got {}",
        cap_loop.half_edges.len()
    );

    let v = bb.signed_volume();
    assert!(
        (v - area * 0.8_f64).abs() < VOL_EPS,
        "below H volume {v}, expected {}",
        area * 0.8_f64
    );
}

// ── Box with a through hole cut across the hole (annulus cap) ───────────────

/// `[0,3]×[0,3]×[0,1]` box pierced by a `[1,2]×[1,2]` square hole through z.
/// (Same fixture as tests/topology.rs.)
fn build_box_with_hole() -> (Brep, Id<Solid>) {
    let mut b = SolidBuilder::new();
    let lo = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let hi = |x: f64, y: f64| Point3::new(x, y, 1.0);

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

    // Cap loops follow the same orientation as tests/topology.rs (a watertight,
    // genus-1 solid): the hole loop is wound opposite to the outer loop on each
    // cap so the inner walls' shared edges get reversed siblings.
    let bottom_outer: Vec<Point3> = outer.iter().rev().map(|&(x, y)| lo(x, y)).collect();
    let bottom_hole: Vec<Point3> = hole.iter().map(|&(x, y)| lo(x, y)).collect();
    let f = b.face_with_holes(
        &bottom_outer,
        &[&bottom_hole],
        plane(lo(0.0, 0.0), Vec3::new(0.0, 0.0, -1.0)),
    );
    b.push_face(f);

    let top_outer: Vec<Point3> = outer.iter().map(|&(x, y)| hi(x, y)).collect();
    let top_hole: Vec<Point3> = hole.iter().rev().map(|&(x, y)| hi(x, y)).collect();
    let f = b.face_with_holes(
        &top_outer,
        &[&top_hole],
        plane(hi(0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
    );
    b.push_face(f);

    for i in 0..4 {
        let (x0, y0) = outer[i];
        let (x1, y1) = outer[(i + 1) % 4];
        let nrm = Vec3::new(y1 - y0, -(x1 - x0), 0.0);
        let ring = [lo(x0, y0), lo(x1, y1), hi(x1, y1), hi(x0, y0)];
        let f = b.face(&ring, plane(lo(x0, y0), nrm));
        b.push_face(f);
    }

    // Inner hole walls (same orientation as tests/topology.rs): outward normal
    // points into the void.
    for i in 0..4 {
        let (x0, y0) = hole[i];
        let (x1, y1) = hole[(i + 1) % 4];
        let nrm = Vec3::new(-(y1 - y0), x1 - x0, 0.0);
        let ring = [lo(x1, y1), lo(x0, y0), hi(x0, y0), hi(x1, y1)];
        let f = b.face(&ring, plane(lo(x0, y0), nrm));
        b.push_face(f);
    }

    b.finish()
}

#[test]
fn box_with_hole_cut_yields_annulus_cap() {
    let tol = Tol::default();
    let (brep, solid) = build_box_with_hole();
    // The through-hole solid is watertight and genus 1 (proven by
    // tests/topology.rs). Cutting at z = 0.5 across the hole splits its constant
    // (9 − 1) cross-section exactly in half, so V(below) must equal V(total)/2
    // regardless of the fixture's volume sign convention.
    let total = brep.signed_volume();

    let cut_plane = plane(Point3::new(0.0, 0.0, 0.5), Vec3::Z);
    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    let CutResult::Cut { brep: bb, caps } = below else {
        panic!("expected cut");
    };
    bb.validate(&tol, ValidateLevel::Full)
        .expect("annulus cut valid");

    // The cap is an annulus: one cap face with one outer loop and one hole loop.
    // This is the headline acceptance criterion: cutting a through-hole solid
    // across the hole produces an annulus (outer + inner) cap.
    assert_eq!(caps.len(), 1usize, "one cap face (annulus)");
    let cap = bb.topo.faces.get(caps[0]).unwrap();
    assert_eq!(cap.inners.len(), 1usize, "annulus cap has one hole loop");

    // Volume integrity: V(below) + V(above) = V(whole). The cut rebuilds both
    // halves preserving the input face orientations, so the (quirky-but-
    // consistent) hand-fixture convention cancels; the two coincident caps
    // (below's top, above's bottom) carry opposite normals and cancel in the
    // sum. (The hole fixture's absolute signed volume is not asserted here
    // because the divergence-theorem inner-loop convention is not exercised
    // elsewhere; the sum identity is convention-independent.)
    let above = cut(&brep, solid, &cut_plane, KeepSide::Above, &tol).expect("above");
    let ab = above.brep();
    ab.validate(&tol, ValidateLevel::Full).expect("above valid");
    let v_below = bb.signed_volume();
    let v_above = ab.signed_volume();
    assert!(
        (v_below + v_above - total).abs() < VOL_EPS,
        "halves must sum to whole: {v_below} + {v_above} vs {total}"
    );
}

// ── Cylinder cuts ───────────────────────────────────────────────────────────

fn z_axis() -> Line3 {
    Line3::new(Point3::origin(), Vec3::Z).expect("axis")
}

#[test]
fn cylinder_perpendicular_cut_circle_cap() {
    let tol = Tol::default();
    let radius = 0.3_f64;
    let length = 2.0_f64;
    let profile = Profile2d::circle(radius).expect("circle");
    let brep = extrude(&profile, &z_axis(), length, &tol).expect("extrude cyl");
    let solid = brep.solids[0];

    // Cut perpendicular to the axis at z = 0.8. The cap is a circle (two arcs).
    let cut_plane = plane(Point3::new(0.0, 0.0, 0.8), Vec3::Z);
    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    let CutResult::Cut { brep: bb, caps } = below else {
        panic!("expected cut");
    };
    bb.validate(&tol, ValidateLevel::Full)
        .expect("cyl perp cut valid");
    assert_eq!(caps.len(), 1usize, "one circular cap");

    // V = π r² · h_below.
    let expected = PI * radius * radius * 0.8_f64;
    let v = bb.signed_volume();
    assert!(
        (v - expected).abs() < VOL_EPS_CURVED,
        "cyl below V {v}, expected {expected}"
    );
}

#[test]
fn cylinder_oblique_cut_ellipse_cap() {
    let tol = Tol::default();
    let radius = 0.3_f64;
    let length = 3.0_f64;
    let profile = Profile2d::circle(radius).expect("circle");
    let brep = extrude(&profile, &z_axis(), length, &tol).expect("extrude cyl");
    let solid = brep.solids[0];

    // Oblique plane through (0,0,1.5) with normal (0, 1, 1): tilts the section
    // into an ellipse. The plane crosses the whole cylinder (it passes through
    // the axis at z = 1.5, well inside [0, 3]).
    let cut_plane = plane(Point3::new(0.0, 0.0, 1.5), Vec3::new(0.0, 1.0, 1.0));
    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    let CutResult::Cut { brep: bb, caps } = below else {
        panic!("expected cut");
    };
    bb.validate(&tol, ValidateLevel::Full)
        .expect("cyl oblique cut valid");
    assert_eq!(caps.len(), 1usize, "one elliptical cap");

    // The oblique plane passes through the axis at z = 1.5 (the mid-height of
    // the [0, 3] cylinder), so by symmetry it keeps exactly half the cylinder:
    // V = ½·π·r²·L. The divergence-theorem volume integral now integrates the
    // elliptical-rim cylinder patch in closed form (`mass/volume.rs`), so this
    // is a real volume check, not just a geometric one.
    let v = bb.signed_volume();
    let expected_vol = 0.5_f64 * PI * radius * radius * length;
    assert!(
        (v - expected_vol).abs() < VOL_EPS_CURVED,
        "oblique cut V {v}, expected {expected_vol}"
    );

    // The cap edges must be on an ellipse, and its semi-axes match the analytic
    // plane × cylinder section: semi-minor = r, semi-major = r / |cos θ| where θ
    // is the angle between the plane normal and the axis (here 45°, cos = 1/√2).
    let cap = bb.topo.faces.get(caps[0]).unwrap();
    let cap_loop = bb.topo.loops.get(cap.outer).unwrap();
    let mut found_ellipse = false;
    for &he_id in &cap_loop.half_edges {
        let he = bb.topo.half_edges.get(he_id).unwrap();
        if let Some(CurveGeom::Ellipse(e)) = bb.geom.curve(he.curve) {
            found_ellipse = true;
            assert!(
                (e.semi_minor() - radius).abs() < VOL_EPS_CURVED,
                "semi-minor {} should be r {radius}",
                e.semi_minor()
            );
            let cos_theta = 1.0_f64 / 2.0_f64.sqrt();
            let expected_major = radius / cos_theta;
            assert!(
                (e.semi_major() - expected_major).abs() < VOL_EPS_CURVED,
                "semi-major {} should be {expected_major}",
                e.semi_major()
            );
        }
    }
    assert!(found_ellipse, "cap must be bounded by ellipse arcs");
}

// The axis-parallel chord cut of a cylinder: the cut splits the bulging rim arcs
// along the ruling lines (any number of ruling portals per face), keeps the
// segment-side arcs and disk-cap segments, and seals the opening with a
// rectangular bow-cap. The circular-segment disk caps integrate exactly via the
// arc-corrected planar area in `mass/volume.rs`.
#[test]
fn cylinder_axis_parallel_chord_cut() {
    let tol = Tol::default();
    let radius = 0.3_f64;
    let length = 2.0_f64;
    let profile = Profile2d::circle(radius).expect("circle");
    let brep = extrude(&profile, &z_axis(), length, &tol).expect("extrude cyl");
    let solid = brep.solids[0];

    // Axis-parallel plane x = d (normal +x), cutting a chord. Keep x ≤ d (Below):
    // the larger circular segment. d = 0.1 < r = 0.3.
    let d = 0.1_f64;
    let cut_plane = plane(Point3::new(d, 0.0, 0.0), Vec3::X);
    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    let CutResult::Cut { brep: bb, caps } = below else {
        panic!("expected cut");
    };
    bb.validate(&tol, ValidateLevel::Full)
        .expect("cyl chord cut valid");
    assert_eq!(caps.len(), 1usize, "one rectangular bow-cap");

    // The kept piece is a circular segment (the part of the disk with x ≤ d)
    // extruded along the length. Its cross-section area is the area of the
    // larger segment of a circle cut by the chord at signed distance d from the
    // centre on the +x side:
    //   A = r² · acos(d/r) − d · √(r² − d²)   (area where x ≤ d, the larger part)
    // Wait: the segment with x ≤ d (d > 0) is the LARGER part, area =
    //   π r² − [ r² acos(d/r) − d √(r²−d²) ].
    let drr = d / radius;
    let minor_seg = radius * radius * drr.acos() - d * (radius * radius - d * d).sqrt();
    let area = PI * radius * radius - minor_seg;
    let expected = area * length;
    let v = bb.signed_volume();
    assert!(
        (v - expected).abs() < VOL_EPS_CURVED,
        "chord cut V {v}, expected {expected}"
    );
}

// ── Section extraction ──────────────────────────────────────────────────────

#[test]
fn section_of_cube_is_single_square_loop() {
    let tol = Tol::default();
    let (brep, solid) = build_box(2.0, 2.0, 2.0);
    let sec_plane = plane(Point3::new(0.0, 0.0, 1.0), Vec3::Z);
    let result = section(&brep, solid, &sec_plane, &tol).expect("section");

    assert_eq!(result.profiles.len(), 1usize, "one profile");
    assert_eq!(result.loop_count(), 1usize, "one loop total (no holes)");
    let profile = &result.profiles[0];
    assert!(profile.holes.is_empty(), "cube section has no holes");
    assert_eq!(profile.outer.points_2d.len(), 4usize, "square (4 corners)");

    // The 2-D area of the section equals the cross-section: 2 × 2 = 4.
    let area = shoelace(&profile.outer.points_2d).abs();
    assert!((area - 4.0_f64).abs() < VOL_EPS, "section area {area}");
}

#[test]
fn section_of_circular_column_is_arc_loop() {
    // A round member: a circular column sectioned perpendicular to its axis must
    // yield a loop of circular arcs (not a silently-dropped None). The cut splits
    // the rim at the seam, so the section is two semicircular arcs.
    let tol = Tol::default();
    let radius = 0.25_f64;
    let profile = Profile2d::circle(radius).expect("circle");
    let brep = extrude(&profile, &z_axis(), 2.0, &tol).expect("column");
    let solid = brep.solids[0];

    let sec_plane = plane(Point3::new(0.0, 0.0, 1.0), Vec3::Z);
    let result = section(&brep, solid, &sec_plane, &tol).expect("section");

    assert_eq!(result.profiles.len(), 1usize, "one profile");
    assert_eq!(result.loop_count(), 1usize, "one loop (the round outline)");
    let profile = &result.profiles[0];
    assert!(profile.holes.is_empty(), "solid column has no holes");

    // Every boundary edge is a circular arc of the column radius.
    assert!(
        !profile.outer.edges.is_empty(),
        "the round section must expose its arc edges"
    );
    for e in &profile.outer.edges {
        match e {
            SectionEdge::Arc { radius: r, .. } => assert!(
                (r - radius).abs() < VOL_EPS_CURVED,
                "arc radius {r} should be the column radius {radius}"
            ),
            SectionEdge::Line { .. } => panic!("a circular section must be arcs, not segments"),
            _ => panic!("unexpected section edge variant"),
        }
    }
}

#[test]
fn section_of_holed_box_has_outer_and_hole() {
    let tol = Tol::default();
    let (brep, solid) = build_box_with_hole();
    let sec_plane = plane(Point3::new(0.0, 0.0, 0.5), Vec3::Z);
    let result = section(&brep, solid, &sec_plane, &tol).expect("section");

    // One profile (outer square) with exactly one hole (the through hole).
    assert_eq!(result.profiles.len(), 1usize, "one profile");
    let profile = &result.profiles[0];
    assert_eq!(
        profile.holes.len(),
        1usize,
        "the through hole appears as one hole loop"
    );
    assert_eq!(result.loop_count(), 2usize, "outer + 1 hole = 2 loops");

    // Net cross-section area = outer (9) − hole (1) = 8.
    let outer_area = shoelace(&profile.outer.points_2d).abs();
    let hole_area = shoelace(&profile.holes[0].points_2d).abs();
    assert!(
        (outer_area - 9.0_f64).abs() < VOL_EPS,
        "outer area {outer_area}"
    );
    assert!(
        (hole_area - 1.0_f64).abs() < VOL_EPS,
        "hole area {hole_area}"
    );
}

/// Shoelace signed area of a 2-D polygon.
fn shoelace(poly: &[[f64; 2]]) -> f64 {
    let n = poly.len();
    let mut a = 0.0_f64;
    for i in 0..n {
        let p = poly[i];
        let q = poly[(i + 1) % n];
        a += p[0] * q[1] - q[0] * p[1];
    }
    a / 2.0_f64
}
