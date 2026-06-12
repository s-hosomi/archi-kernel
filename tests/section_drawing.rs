//! Phase 4 — section drawing E2E (`DESIGN.md` §6-3, §10 Phase 4).
//!
//! A mini building: two rectangular columns + one circular column (all along
//! `+z`), one girder (along `+x`), and a slab with a rectangular opening and a
//! circular sleeve (an `OpeningSubtraction`). The tests draw its plans
//! (伏図), a coincident-plane section on the slab top, and an elevation (軸組図),
//! and check the holed-profile structure, the arc edges of round members, the
//! coplanar convention, and the 2-D→3-D frame round-trip.
//!
//! Every literal carries an `f64` annotation and an explicit tolerance
//! (`DESIGN.md` §12).

use archi_kernel::boolean::prismatic::{self, ExtrudeLeaf};
use archi_kernel::brep::Brep;
use archi_kernel::csg::{CsgNode, Member, Profile2d};
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::Plane;
use archi_kernel::section::{section, section_members, SectionEdge, SectionLoop};
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::ValidateLevel;

const EPS: f64 = 1e-9;
const EPS_CURVED: f64 = 1e-6;

// ── building members ────────────────────────────────────────────────────────

/// A rectangular column along `+z`: footprint `0.6 × 0.4` (world `x × y`), from
/// `z0` up `h`, centred at `(cx, cy)`. (`Profile2d::rect(half_w, half_h)` maps
/// `half_w → world Y`, `half_h → world X` for a `+z` extrusion.)
fn rect_column(cx: f64, cy: f64, z0: f64, h: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(0.2_f64, 0.3_f64).expect("rect"),
        origin: Point3::new(cx, cy, z0),
        axis: Vec3::Z,
        length: h,
    }
}

/// A circular column along `+z`, radius `r`, centred at `(cx, cy)`, from `z0` up `h`.
fn circ_column(cx: f64, cy: f64, z0: f64, h: f64, r: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::circle(r).expect("circle"),
        origin: Point3::new(cx, cy, z0),
        axis: Vec3::Z,
        length: h,
    }
}

/// A girder along `+x`: cross-section `0.3` (world Y) × `0.5` (world Z), running
/// from `x0` for `len`, its section centred at `(cy, cz)`.
fn girder(x0: f64, cy: f64, cz: f64, len: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(0.15_f64, 0.25_f64).expect("rect"),
        origin: Point3::new(x0, cy, cz),
        axis: Vec3::X,
        length: len,
    }
}

/// The slab base: a flat box `4 (x) × 4 (y) × 0.2 (z)` from `z0`, corner at origin.
fn slab_base(z0: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(2.0_f64, 2.0_f64).expect("rect"),
        origin: Point3::new(2.0_f64, 2.0_f64, z0),
        axis: Vec3::Z,
        length: 0.2_f64,
    }
}

/// A rectangular vertical opening through the slab, `0.8 (x) × 0.6 (y)`, centred
/// at `(cx, cy)`.
fn rect_opening(cx: f64, cy: f64, z0: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(0.3_f64, 0.4_f64).expect("rect"),
        origin: Point3::new(cx, cy, z0 - 0.1_f64),
        axis: Vec3::Z,
        length: 0.4_f64,
    }
}

/// A circular vertical sleeve through the slab, radius `r`, centred at `(cx, cy)`.
fn circ_sleeve(cx: f64, cy: f64, z0: f64, r: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::circle(r).expect("circle"),
        origin: Point3::new(cx, cy, z0 - 0.1_f64),
        axis: Vec3::Z,
        length: 0.4_f64,
    }
}

/// Build a single-solid B-rep from one extruded leaf (for direct `section`).
fn leaf_brep(leaf: &ExtrudeLeaf, tol: &Tol) -> Brep {
    use archi_kernel::build::extrude;
    use archi_kernel::Line3;
    let line = Line3::new(leaf.origin, leaf.axis).expect("axis");
    let b = extrude(&leaf.profile, &line, leaf.length, tol).expect("extrude");
    b.validate(tol, ValidateLevel::Full).expect("valid leaf");
    b
}

/// A horizontal cut plane at height `z`, normal `+z`.
fn horizontal(z: f64) -> Plane {
    Plane::new(Point3::new(0.0_f64, 0.0_f64, z), Vec3::Z).expect("plane")
}

/// Absolute shoelace area of a 2-D ring.
fn ring_area(p: &[[f64; 2]]) -> f64 {
    let n = p.len();
    let mut a = 0.0_f64;
    for i in 0..n {
        let q = p[(i + 1) % n];
        a += p[i][0] * q[1] - q[0] * p[i][1];
    }
    (a / 2.0_f64).abs()
}

/// `true` if a loop has any arc edge.
fn has_arc(l: &SectionLoop) -> bool {
    l.edges.iter().any(|e| matches!(e, SectionEdge::Arc { .. }))
}

// ── 伏図: column plan ─────────────────────────────────────────────────────────

#[test]
fn plan_through_columns_shows_each_cross_section() {
    let tol = Tol::default();
    // Two rectangular columns + one circular column, all standing on z∈[0,3].
    let c1 = rect_column(1.0_f64, 1.0_f64, 0.0_f64, 3.0_f64);
    let c2 = rect_column(3.0_f64, 1.0_f64, 0.0_f64, 3.0_f64);
    let c3 = circ_column(2.0_f64, 3.0_f64, 0.0_f64, 3.0_f64, 0.25_f64);

    let plane = horizontal(1.5_f64); // mid-height

    // Each column drawn independently (伏図 = overlay of member sections).
    for col in [c1, c2] {
        let brep = leaf_brep(&col, &tol);
        let res = section(&brep, brep.solids[0], &plane, &tol).expect("section");
        assert_eq!(res.profiles.len(), 1usize, "one profile per column");
        let prof = &res.profiles[0];
        assert!(prof.holes.is_empty(), "a solid column has no holes");
        assert!(!has_arc(&prof.outer), "rectangular column → straight edges");
        // Cross-section area 0.6 × 0.4 = 0.24.
        let a = ring_area(&prof.outer.points_2d);
        assert!((a - 0.24_f64).abs() < EPS, "rect column area {a}");
    }

    let brep = leaf_brep(&c3, &tol);
    let res = section(&brep, brep.solids[0], &plane, &tol).expect("section");
    assert_eq!(res.profiles.len(), 1usize, "one circular profile");
    let prof = &res.profiles[0];
    assert!(prof.holes.is_empty());
    assert!(has_arc(&prof.outer), "circular column → arc edges");
    for e in &prof.outer.edges {
        if let SectionEdge::Arc { radius, .. } = e {
            assert!(
                (radius - 0.25_f64).abs() < EPS_CURVED,
                "arc radius {radius} == column radius"
            );
        }
    }
}

// ── スラブ伏図: holed plan with correct nesting ───────────────────────────────

/// The slab as an `OpeningSubtraction` of a rectangular opening and a circular
/// sleeve. Evaluated via the prismatic engine.
fn slab_brep(tol: &Tol) -> Brep {
    let base = slab_base(0.0_f64);
    let rect = rect_opening(1.0_f64, 1.0_f64, 0.0_f64);
    let circ = circ_sleeve(3.0_f64, 3.0_f64, 0.0_f64, 0.3_f64);
    let b =
        prismatic::opening_subtraction(&base, &[rect, circ], tol).expect("slab with two openings");
    b.validate(tol, ValidateLevel::Full).expect("valid slab");
    b
}

#[test]
fn slab_plan_has_two_holes_correctly_nested() {
    let tol = Tol::default();
    let slab = slab_brep(&tol);
    let plane = horizontal(0.1_f64); // mid-slab

    let res = section(&slab, slab.solids[0], &plane, &tol).expect("section");
    assert_eq!(res.profiles.len(), 1usize, "one slab outline");
    let prof = &res.profiles[0];
    assert_eq!(prof.holes.len(), 2usize, "rectangular + circular hole");

    // Outer area 4×4 = 16; holes: rect 0.8×0.6 = 0.48, circle πr² (r=0.3).
    let outer = ring_area(&prof.outer.points_2d);
    assert!((outer - 16.0_f64).abs() < EPS, "slab outer area {outer}");

    let rect_hole = prof
        .holes
        .iter()
        .find(|h| !has_arc(h))
        .expect("a straight (rectangular) hole");
    let circ_hole = prof
        .holes
        .iter()
        .find(|h| has_arc(h))
        .expect("an arc (circular) hole");

    let rect_a = ring_area(&rect_hole.points_2d);
    assert!((rect_a - 0.48_f64).abs() < EPS, "rect hole area {rect_a}");
    for e in &circ_hole.edges {
        if let SectionEdge::Arc { radius, .. } = e {
            assert!(
                (radius - 0.3_f64).abs() < EPS_CURVED,
                "sleeve arc radius {radius}"
            );
        }
    }
}

// ── スラブ天端一致断面: the coplanar convention ───────────────────────────────

#[test]
fn slab_top_coincident_section_draws_the_plan() {
    let tol = Tol::default();
    let slab = slab_brep(&tol);

    // Cut exactly on the slab top (z = 0.2), normal +z. The top face's outward
    // normal is +z = +normal ⇒ included (the slab's plan, holes and all).
    let top = horizontal(0.2_f64);
    let res = section(&slab, slab.solids[0], &top, &tol).expect("section on top");
    assert_eq!(res.profiles.len(), 1usize, "the slab plan at the top face");
    let prof = &res.profiles[0];
    assert_eq!(prof.holes.len(), 2usize, "openings show as holes");
    let outer = ring_area(&prof.outer.points_2d);
    assert!(
        (outer - 16.0_f64).abs() < EPS,
        "top plan outer area {outer}"
    );

    // Cut exactly on the slab bottom (z = 0), normal +z. The bottom face's
    // outward normal is −z, which opposes +normal ⇒ NOT drawn.
    let bottom = horizontal(0.0_f64);
    let res = section(&slab, slab.solids[0], &bottom, &tol).expect("section on bottom");
    assert!(
        res.profiles.is_empty(),
        "bottom face (normal −z) opposes +normal: nothing drawn, got {:?}",
        res.profiles.len()
    );
}

// ── 軸組図: vertical elevation through columns and the girder ──────────────────

#[test]
fn elevation_through_columns_and_girder() {
    let tol = Tol::default();
    // A rectangular column on z∈[0,3] at (1,1) and a girder along +x at height
    // z=2.75 crossing it. A vertical plane y = 1 (normal +y) cuts both.
    let col = rect_column(1.0_f64, 1.0_f64, 0.0_f64, 3.0_f64);
    let gird = girder(0.0_f64, 1.0_f64, 2.75_f64, 4.0_f64);

    let plane = Plane::new(Point3::new(0.0_f64, 1.0_f64, 0.0_f64), Vec3::Y).expect("plane");

    // Column elevation: a tall rectangle (its width × full height 3.0).
    let cb = leaf_brep(&col, &tol);
    let res = section(&cb, cb.solids[0], &plane, &tol).expect("column elevation");
    assert_eq!(res.profiles.len(), 1usize);
    let prof = &res.profiles[0];
    assert!(
        !has_arc(&prof.outer),
        "rect column elevation is a rectangle"
    );
    // Column cross-section at y=1: world-X extent 0.6, full height 3.0 ⇒ area 1.8.
    let a = ring_area(&prof.outer.points_2d);
    assert!((a - 1.8_f64).abs() < EPS, "column elevation area {a}");

    // Girder elevation: the y = 1 plane runs along the girder's full length 4.0
    // and through its height (world-Z extent 0.3 for a +x extrusion of
    // rect(0.15, 0.25)) ⇒ a 4.0 × 0.3 rectangle, area 1.2.
    let gb = leaf_brep(&gird, &tol);
    let res = section(&gb, gb.solids[0], &plane, &tol).expect("girder elevation");
    assert_eq!(res.profiles.len(), 1usize);
    let prof = &res.profiles[0];
    let a = ring_area(&prof.outer.points_2d);
    assert!((a - 1.2_f64).abs() < EPS, "girder elevation area {a}");

    // A vertical plane through the circular column yields a rectangle (full
    // width chord through the centre): width 2r = 0.5, height 3.0 ⇒ area 1.5.
    let circ = circ_column(2.0_f64, 3.0_f64, 0.0_f64, 3.0_f64, 0.25_f64);
    let circ_plane = Plane::new(Point3::new(0.0_f64, 3.0_f64, 0.0_f64), Vec3::Y).expect("plane");
    let cb = leaf_brep(&circ, &tol);
    let res = section(&cb, cb.solids[0], &circ_plane, &tol).expect("circ column elevation");
    assert_eq!(res.profiles.len(), 1usize);
    let a = ring_area(&res.profiles[0].outer.points_2d);
    assert!(
        (a - 1.5_f64).abs() < EPS_CURVED,
        "circular column full-width elevation area {a}"
    );
}

// ── frame 3-D restoration ─────────────────────────────────────────────────────

#[test]
fn frame_restores_2d_to_original_3d_points() {
    let tol = Tol::default();
    let col = rect_column(1.0_f64, 1.0_f64, 0.0_f64, 3.0_f64);
    let brep = leaf_brep(&col, &tol);
    let plane = horizontal(1.5_f64);
    let res = section(&brep, brep.solids[0], &plane, &tol).expect("section");
    let frame = res.frame.expect("a frame");
    let prof = &res.profiles[0];
    // Every 2-D point maps back to its world 3-D point.
    for (p2, p3) in prof.outer.points_2d.iter().zip(&prof.outer.points_3d) {
        let back = frame.to_3d(*p2);
        let d = (back - *p3).norm();
        assert!(d < EPS, "2-D→3-D round-trip error {d}");
        // And projecting the 3-D point reproduces the 2-D coordinate.
        let fwd = frame.project(*p3);
        assert!(
            (fwd[0] - p2[0]).abs() < EPS && (fwd[1] - p2[1]).abs() < EPS,
            "3-D→2-D projection mismatch"
        );
    }
}

// ── section_members: per-member isolation ─────────────────────────────────────

#[test]
fn section_members_returns_per_member_results_and_isolates_failure() {
    let tol = Tol::default();
    let plane = horizontal(1.5_f64);

    // One good member (a column) and one that fails to evaluate (an H×H general
    // boolean with no common prismatic direction).
    let good = Member::new(CsgNode::Extrude {
        profile: Profile2d::rect(0.2_f64, 0.3_f64).expect("rect"),
        origin: Point3::new(1.0_f64, 1.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 3.0_f64,
    });
    // A difference of two H-sections crossed at a right angle has no common
    // prismatic direction ⇒ EvalError, isolated per member (not skipped).
    let bad = Member::new(CsgNode::Difference {
        positive: Box::new(CsgNode::Extrude {
            profile: Profile2d::h_section(0.3_f64, 0.2_f64, 0.02_f64, 0.02_f64).expect("h"),
            origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
            axis: Vec3::Z,
            length: 1.0_f64,
        }),
        negative: Box::new(CsgNode::Extrude {
            profile: Profile2d::h_section(0.3_f64, 0.2_f64, 0.02_f64, 0.02_f64).expect("h"),
            origin: Point3::new(0.0_f64, 0.0_f64, 0.5_f64),
            axis: Vec3::X,
            length: 1.0_f64,
        }),
    });

    let mut members = vec![good, bad];
    let results = section_members(&mut members, &plane, &tol);
    assert_eq!(results.len(), 2usize, "one result slot per member");
    assert!(results[0].is_ok(), "the good column sections");
    assert_eq!(
        results[0].as_ref().unwrap().profiles.len(),
        1usize,
        "good column has one profile"
    );
    assert!(
        results[1].is_err(),
        "the failing member surfaces its EvalError, not a skip"
    );
}
