//! Regression tests for the coplanar / arc-bearing cut and section bugs in
//! `boolean::half_space`.
//!
//! Two root causes are covered:
//!
//! 1. `nest_cap_cycles` measured cap-loop area and containment from the *chord*
//!    polygon only, so a circular cap (whose chord ring is just two points — two
//!    semicircular arcs) collapsed to a degenerate 2-gon: area 0, containment
//!    always false. A hollow tube cut across the bore therefore failed to nest
//!    the inner bore circle as a hole and produced two independent disk caps
//!    instead of one annular cap. The fix nests with the arc-corrected loop area
//!    and an arc-sampled containment polygon.
//!
//! 2. `register_on_boundary_caps` registered only *straight* On–On boundary edges
//!    as cap edges and silently dropped circular / elliptical ones. A horizontal
//!    section coincident with an arc-bearing coplanar lid (the rim of a round
//!    void / sleeve / clip column) therefore left the lid's arc rim unsealed,
//!    producing a non-watertight brep and an `InvalidResult` (`MissingSibling`).
//!    The fix registers arc On–On edges on the shared section conic with the same
//!    signed-multiplicity cancellation the straight edges use.

use std::collections::HashMap;
use std::f64::consts::PI;

use archi_kernel::boolean::prismatic::{self, ExtrudeLeaf};
use archi_kernel::boolean::{cut, CutResult, KeepSide};
use archi_kernel::csg::Profile2d;
use archi_kernel::geom::{CurveGeom, CurveId, SurfaceGeom, VertexGeom};
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::{Circle3, Cylinder, Line3, Plane};
use archi_kernel::section::section;
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::arena::Id;
use archi_kernel::topo::{Face, HalfEdge, Loop, Sense, Shell, Solid, ValidateLevel, Vertex};
use archi_kernel::Brep;

/// Volume tolerance for curved members (the tube is faceted only in area
/// integrals, but its analytic volume is exact via the divergence form).
const VOL_EPS_CURVED: f64 = 1.0e-6_f64;

fn horizontal(z: f64) -> Plane {
    Plane::new(Point3::new(0.0_f64, 0.0_f64, z), Vec3::Z).expect("horizontal plane")
}

// ── Fix 1: annular cap nesting for a circular outer rim ────────────────────

type CoordKey = (i64, i64, i64);

fn ckey(p: Point3) -> CoordKey {
    let q = |x: f64| (x * 1.0e9_f64).round() as i64;
    (q(p.x), q(p.y), q(p.z))
}

fn z_axis() -> Line3 {
    Line3::new(Point3::origin(), Vec3::Z).expect("z axis")
}

fn plane(point: Point3, normal: Vec3) -> Plane {
    Plane::new(point, normal).expect("plane")
}

/// A watertight genus-1 hollow cylinder, hand-built so the inner *and* outer cap
/// rims are circular arcs (the degenerate 2-point chord ring that the old
/// `nest_cap_cycles` collapsed to area 0). Outer radius `r_out`, bore `r_in`,
/// axis +Z, z ∈ `[0, length]`. The annular caps each carry an inner hole loop.
///
/// Ported from the verifier's reproduction fixture (validated at
/// `ValidateLevel::Full`).
fn build_hollow_cylinder(r_out: f64, r_in: f64, length: f64) -> (Brep, Id<Solid>) {
    let tol = Tol::default();
    let mut brep = Brep::new();
    let two_pi = 2.0 * PI;
    let axis = z_axis();
    let dir = Vec3::Z;

    let mut vmap: HashMap<CoordKey, Id<Vertex>> = HashMap::new();
    let mut lmap: HashMap<(CoordKey, CoordKey), (CurveId, Point3)> = HashMap::new();
    let vertex =
        |brep: &mut Brep, vmap: &mut HashMap<CoordKey, Id<Vertex>>, p: Point3| -> Id<Vertex> {
            let k = ckey(p);
            if let Some(&v) = vmap.get(&k) {
                return v;
            }
            let pid = brep.geom.insert_point(VertexGeom::Explicit(p));
            let v = brep.topo.add_vertex(Vertex { point: pid });
            vmap.insert(k, v);
            v
        };

    let mut faces: Vec<Id<Face>> = Vec::new();

    let make_circles = |brep: &mut Brep, radius: f64| -> (Circle3, Circle3, CurveId, CurveId) {
        let bottom_circle = Circle3::new(Point3::new(0.0, 0.0, 0.0), dir, radius).expect("circle");
        let top_circle = Circle3::new(Point3::new(0.0, 0.0, length), dir, radius).expect("circle");
        let cb = brep.geom.insert_curve(CurveGeom::Circle(bottom_circle));
        let ct = brep.geom.insert_curve(CurveGeom::Circle(top_circle));
        (bottom_circle, top_circle, cb, ct)
    };

    let arc = |brep: &mut Brep, curve: CurveId, start: Id<Vertex>, a: f64, b: f64| {
        brep.topo.add_half_edge(HalfEdge {
            start,
            curve,
            boundary: [a, b],
        })
    };
    let line_he = |brep: &mut Brep,
                   vmap: &mut HashMap<CoordKey, Id<Vertex>>,
                   lmap: &mut HashMap<(CoordKey, CoordKey), (CurveId, Point3)>,
                   a: Point3,
                   b: Point3|
     -> Id<HalfEdge> {
        let (ka, kb) = (ckey(a), ckey(b));
        let unordered = if ka <= kb { (ka, kb) } else { (kb, ka) };
        let (cid, origin) = if let Some(&entry) = lmap.get(&unordered) {
            entry
        } else {
            let origin = if ka <= kb { a } else { b };
            let other = if ka <= kb { b } else { a };
            let line = Line3::new(origin, other - origin).expect("line");
            let cid = brep.geom.insert_curve(CurveGeom::Line(line));
            lmap.insert(unordered, (cid, origin));
            (cid, origin)
        };
        let line = match brep.geom.curve(cid).expect("line curve") {
            CurveGeom::Line(l) => *l,
            _ => unreachable!(),
        };
        let va = vertex(brep, vmap, a);
        let _ = vertex(brep, vmap, b);
        let ta = (a - origin).dot(line.dir().as_vec());
        let tb = (b - origin).dot(line.dir().as_vec());
        brep.topo.add_half_edge(HalfEdge {
            start: va,
            curve: cid,
            boundary: [ta, tb],
        })
    };

    // Outer wall (r_out): outward normal (sense Same), two half-cylinders.
    let (obc, otc, ocb, oct) = make_circles(&mut brep, r_out);
    let a_bot = obc.point_at(0.0);
    let b_bot = obc.point_at(PI);
    let a_top = otc.point_at(0.0);
    let b_top = otc.point_at(PI);
    let va_bot = vertex(&mut brep, &mut vmap, a_bot);
    let vb_bot = vertex(&mut brep, &mut vmap, b_bot);
    let va_top = vertex(&mut brep, &mut vmap, a_top);
    let vb_top = vertex(&mut brep, &mut vmap, b_top);
    let cyl_out = Cylinder::new(axis, r_out).expect("cyl");
    {
        let h0 = arc(&mut brep, ocb, va_bot, 0.0, PI);
        let h1 = line_he(&mut brep, &mut vmap, &mut lmap, b_bot, b_top);
        let h2 = arc(&mut brep, oct, vb_top, PI, 0.0);
        let h3 = line_he(&mut brep, &mut vmap, &mut lmap, a_top, a_bot);
        let lp = brep.topo.add_loop(Loop {
            half_edges: vec![h0, h1, h2, h3],
        });
        let surface = brep.geom.insert_surface(SurfaceGeom::Cylinder(cyl_out));
        faces.push(brep.topo.add_face(Face {
            surface,
            sense: Sense::Same,
            outer: lp,
            inners: Vec::new(),
        }));
    }
    {
        let h0 = arc(&mut brep, ocb, vb_bot, PI, two_pi);
        let h1 = line_he(&mut brep, &mut vmap, &mut lmap, a_bot, a_top);
        let h2 = arc(&mut brep, oct, va_top, two_pi, PI);
        let h3 = line_he(&mut brep, &mut vmap, &mut lmap, b_top, b_bot);
        let lp = brep.topo.add_loop(Loop {
            half_edges: vec![h0, h1, h2, h3],
        });
        let surface = brep.geom.insert_surface(SurfaceGeom::Cylinder(cyl_out));
        faces.push(brep.topo.add_face(Face {
            surface,
            sense: Sense::Same,
            outer: lp,
            inners: Vec::new(),
        }));
    }

    // Inner wall (r_in): normal points toward the axis (sense Reversed).
    let (ibc, itc, icb, ict) = make_circles(&mut brep, r_in);
    let ia_bot = ibc.point_at(0.0);
    let ib_bot = ibc.point_at(PI);
    let ia_top = itc.point_at(0.0);
    let ib_top = itc.point_at(PI);
    let iva_bot = vertex(&mut brep, &mut vmap, ia_bot);
    let ivb_bot = vertex(&mut brep, &mut vmap, ib_bot);
    let iva_top = vertex(&mut brep, &mut vmap, ia_top);
    let ivb_top = vertex(&mut brep, &mut vmap, ib_top);
    let cyl_in = Cylinder::new(axis, r_in).expect("cyl");
    {
        let h0 = arc(&mut brep, icb, iva_bot, 0.0, PI);
        let h1 = line_he(&mut brep, &mut vmap, &mut lmap, ib_bot, ib_top);
        let h2 = arc(&mut brep, ict, ivb_top, PI, 0.0);
        let h3 = line_he(&mut brep, &mut vmap, &mut lmap, ia_top, ia_bot);
        let lp = brep.topo.add_loop(Loop {
            half_edges: vec![h0, h1, h2, h3],
        });
        let surface = brep.geom.insert_surface(SurfaceGeom::Cylinder(cyl_in));
        faces.push(brep.topo.add_face(Face {
            surface,
            sense: Sense::Reversed,
            outer: lp,
            inners: Vec::new(),
        }));
    }
    {
        let h0 = arc(&mut brep, icb, ivb_bot, PI, two_pi);
        let h1 = line_he(&mut brep, &mut vmap, &mut lmap, ia_bot, ia_top);
        let h2 = arc(&mut brep, ict, iva_top, two_pi, PI);
        let h3 = line_he(&mut brep, &mut vmap, &mut lmap, ib_top, ib_bot);
        let lp = brep.topo.add_loop(Loop {
            half_edges: vec![h0, h1, h2, h3],
        });
        let surface = brep.geom.insert_surface(SurfaceGeom::Cylinder(cyl_in));
        faces.push(brep.topo.add_face(Face {
            surface,
            sense: Sense::Reversed,
            outer: lp,
            inners: Vec::new(),
        }));
    }

    // Bottom annular cap (z=0, outward −Z).
    {
        let h0 = arc(&mut brep, ocb, va_bot, two_pi, PI);
        let h1 = arc(&mut brep, ocb, vb_bot, PI, 0.0);
        let outer = brep.topo.add_loop(Loop {
            half_edges: vec![h0, h1],
        });
        let hi0 = arc(&mut brep, icb, ivb_bot, PI, 0.0);
        let hi1 = arc(&mut brep, icb, iva_bot, two_pi, PI);
        let inner = brep.topo.add_loop(Loop {
            half_edges: vec![hi0, hi1],
        });
        let (surface, flipped) = brep
            .geom
            .insert_plane(plane(Point3::new(0.0, 0.0, 0.0), -dir), &tol);
        let sense = if flipped {
            Sense::Reversed
        } else {
            Sense::Same
        };
        faces.push(brep.topo.add_face(Face {
            surface,
            sense,
            outer,
            inners: vec![inner],
        }));
    }
    // Top annular cap (z=length, outward +Z).
    {
        let h0 = arc(&mut brep, oct, va_top, 0.0, PI);
        let h1 = arc(&mut brep, oct, vb_top, PI, two_pi);
        let outer = brep.topo.add_loop(Loop {
            half_edges: vec![h0, h1],
        });
        let hi0 = arc(&mut brep, ict, iva_top, 0.0, PI);
        let hi1 = arc(&mut brep, ict, ivb_top, PI, two_pi);
        let inner = brep.topo.add_loop(Loop {
            half_edges: vec![hi0, hi1],
        });
        let (surface, flipped) = brep
            .geom
            .insert_plane(plane(Point3::new(0.0, 0.0, length), dir), &tol);
        let sense = if flipped {
            Sense::Reversed
        } else {
            Sense::Same
        };
        faces.push(brep.topo.add_face(Face {
            surface,
            sense,
            outer,
            inners: vec![inner],
        }));
    }

    let shell = brep.topo.add_shell(Shell {
        faces: faces.clone(),
    });
    let solid = brep.topo.add_solid(Solid {
        shells: vec![shell],
    });
    brep.solids = vec![solid];
    (brep, solid)
}

#[test]
fn hollow_tube_cut_yields_annulus_cap() {
    let tol = Tol::default();
    let r_out = 0.5_f64;
    let r_in = 0.25_f64;
    let length = 2.0_f64;
    let (brep, solid) = build_hollow_cylinder(r_out, r_in, length);
    brep.validate(&tol, ValidateLevel::Full)
        .expect("tube fixture valid");
    let total = brep.signed_volume().abs();

    // Cut horizontally through the middle of the tube.
    let cut_plane = horizontal(1.0_f64);
    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below");
    let CutResult::Cut { brep: bb, caps } = below else {
        panic!("expected a real cut");
    };
    bb.validate(&tol, ValidateLevel::Full)
        .expect("annulus tube cut valid");

    // Headline: the cut yields exactly *one* annular cap — one outer circle rim
    // plus one nested bore-hole loop. Before the fix the bore circle failed to
    // nest (its 2-point chord ring had area 0 and contained nothing), so it was
    // emitted as a second independent disk cap (caps.len()==2, both inners==0).
    assert_eq!(caps.len(), 1_usize, "one annular cap face");
    let cap = bb.topo.faces.get(caps[0]).expect("cap face");
    assert_eq!(
        cap.inners.len(),
        1_usize,
        "annulus cap has exactly one bore hole loop"
    );

    // Volume integrity, convention-independent: V(below) + V(above) = V(whole).
    // The two halves preserve the input face orientations and the two coincident
    // caps (below's top, above's bottom) carry opposite normals and cancel in the
    // sum, so the identity holds regardless of the fixture's inner-loop sign
    // convention. (Absolute signed volume is *not* asserted because the
    // divergence-theorem inner-loop convention of a hand-built reversed-sense bore
    // is the same one the existing `box_with_hole_cut_yields_annulus_cap` test
    // deliberately leaves out of its absolute check.)
    let v_below = bb.signed_volume();
    let above = cut(&brep, solid, &cut_plane, KeepSide::Above, &tol).expect("above");
    let ab = above.brep();
    ab.validate(&tol, ValidateLevel::Full)
        .expect("above half valid");
    let v_above = ab.signed_volume();
    assert!(
        (v_below + v_above - total).abs() < VOL_EPS_CURVED,
        "halves must sum to whole: {v_below} + {v_above} vs {total}"
    );
}

// ── Fix 2: section coincident with an arc-bearing coplanar lid ─────────────

/// A vertical wall box, extruded along +z.
fn wall() -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(6.4_f64, 0.09_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 7.8_f64),
        axis: Vec3::Z,
        length: 2.9_f64,
    }
}

/// Rectangular window void: cx=3.0, z0=8.65, 1.35 wide (x), 0.38 deep (y).
fn window_void() -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(0.675_f64, 0.19_f64).expect("rect"),
        origin: Point3::new(3.0_f64, 0.0_f64, 8.65_f64),
        axis: Vec3::Z,
        length: 1.6_f64,
    }
}

/// Round column clipper: cx=0,cy=0, r=0.325, z0=7.5, h=3.5.
fn round_column() -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::circle(0.325_f64).expect("circle"),
        origin: Point3::new(0.0_f64, 0.0_f64, 7.5_f64),
        axis: Vec3::Z,
        length: 3.5_f64,
    }
}

/// Rectangular column clipper (control): same footprint span, no arc.
fn rect_column() -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(0.325_f64, 0.325_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 7.5_f64),
        axis: Vec3::Z,
        length: 3.5_f64,
    }
}

#[test]
fn arc_lid_section_at_sill_round_column() {
    let tol = Tol::default();
    let brep = prismatic::clip(&wall(), &[window_void()], &[round_column()], &tol)
        .expect("clip round column");
    let solid = brep.solids[0];

    // Exactly on the window sill z = 8.65 (coincident with the void's bottom
    // lid, whose rim now carries the column's circular arc). This used to fail
    // with InvalidResult(MissingSibling); after the fix it must seal cleanly.
    let exact = section(&brep, solid, &horizontal(8.65_f64), &tol);
    assert!(
        exact.is_ok(),
        "exact-coincident arc-lid section must seal: {exact:?}"
    );

    // Nudged ±1mm must also succeed (the non-coincident chord-cut path).
    let below = section(&brep, solid, &horizontal(8.649_f64), &tol);
    let above = section(&brep, solid, &horizontal(8.651_f64), &tol);
    assert!(below.is_ok(), "z=8.649 nudge should succeed: {below:?}");
    assert!(above.is_ok(), "z=8.651 nudge should succeed: {above:?}");
}

#[test]
fn arc_lid_section_at_sill_rect_column_control() {
    let tol = Tol::default();
    let brep = prismatic::clip(&wall(), &[window_void()], &[rect_column()], &tol)
        .expect("clip rect column");
    let solid = brep.solids[0];
    let exact = section(&brep, solid, &horizontal(8.65_f64), &tol);
    assert!(
        exact.is_ok(),
        "control: rect clip section at sill should succeed: {exact:?}"
    );
}

#[test]
fn arc_cap_section_c2_column_cap_z() {
    // C2: no window; subtract a 1.0m-tall round column from inside the wall.
    // The column's own top/bottom caps land at z = 8.0 / 9.0, producing
    // arc-bearing coplanar lids inside the member.
    let tol = Tol::default();
    let col = ExtrudeLeaf {
        profile: Profile2d::circle(0.325_f64).expect("circle"),
        origin: Point3::new(0.0_f64, 0.0_f64, 8.0_f64),
        axis: Vec3::Z,
        length: 1.0_f64,
    };
    let brep = prismatic::clip(&wall(), &[], &[col], &tol).expect("clip column");
    let solid = brep.solids[0];

    let exact = section(&brep, solid, &horizontal(8.0_f64), &tol);
    assert!(
        exact.is_ok(),
        "C2 exact column-cap-z section must seal: {exact:?}"
    );
    let nudged = section(&brep, solid, &horizontal(8.01_f64), &tol);
    assert!(
        nudged.is_ok(),
        "C2 nudged z=8.01 should succeed: {nudged:?}"
    );
}

// ── Fix 2 sibling: blind circular sleeve in a slab ─────────────────────────

/// slab box centred at (2,2,0), 4x4 footprint, 0.22 tall.
fn slab() -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(2.0_f64, 2.0_f64).expect("rect"),
        origin: Point3::new(2.0_f64, 2.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 0.22_f64,
    }
}

/// blind sleeve cyl(cx=2,cy=2,z0=0.1,r=0.3,h=0.2): bottom cap at z=0.10.
fn blind_sleeve() -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::circle(0.3_f64).expect("circle"),
        origin: Point3::new(2.0_f64, 2.0_f64, 0.1_f64),
        axis: Vec3::Z,
        length: 0.2_f64,
    }
}

#[test]
fn blind_sleeve_bottom_cap_section_seals() {
    let tol = Tol::default();
    let b = prismatic::opening_subtraction(&slab(), &[blind_sleeve()], &tol)
        .expect("opening_subtraction builds");
    b.validate(&tol, ValidateLevel::Full)
        .expect("input solid is watertight");
    let solid = b.solids[0];

    let r_below = section(&b, solid, &horizontal(0.099_f64), &tol);
    let r_exact = section(&b, solid, &horizontal(0.10_f64), &tol);
    let r_above = section(&b, solid, &horizontal(0.101_f64), &tol);

    assert!(r_below.is_ok(), "z=0.099 must be Ok: {r_below:?}");
    assert!(r_above.is_ok(), "z=0.101 must be Ok: {r_above:?}");
    assert!(
        r_exact.is_ok(),
        "z=0.100 (exact sleeve bottom cap) must now seal: {r_exact:?}"
    );
}
