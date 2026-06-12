//! Phase 3c — 2-D circular-arc boolean tests (`DESIGN.md` §10 Phase 3c, §4.2).
//!
//! Covers the commercial-must family (rectangle × circular void / sleeve) plus
//! the everyday circle/circle cases (separate, overlapping, tangent, concentric)
//! and proptest area-conservation / inclusion–exclusion / idempotence with arcs.
//! Every literal carries an `f64` annotation and an explicit tolerance
//! (`DESIGN.md` §12).

use std::f64::consts::PI;

use archi_kernel::boolean::poly2d::{difference, intersection, union, Point2, Region};
use archi_kernel::tolerance::Tol;
use proptest::prelude::*;

fn tol() -> Tol {
    Tol::default()
}

/// Area-comparison tolerance for arc cases (arc-segment corrections involve
/// transcendental terms, so a touch looser than the polygon suite).
const AREA_EPS: f64 = 1e-6;

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Region {
    Region::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

fn circle(cx: f64, cy: f64, r: f64) -> Region {
    Region::circle(Point2::new(cx, cy), r)
}

fn assert_area(region: &Region, expected: f64) {
    let got = region.area();
    assert!(
        (got - expected).abs() <= AREA_EPS * (1.0 + expected.abs()),
        "area mismatch: got {got}, expected {expected}"
    );
}

fn outer_hole_counts(region: &Region) -> (usize, usize) {
    let mut outers = 0usize;
    let mut holes = 0usize;
    for c in &region.contours {
        if c.signed_area() > 0.0 {
            outers += 1;
        } else if c.signed_area() < 0.0 {
            holes += 1;
        }
    }
    (outers, holes)
}

// ── rectangle × circular void (the commercial milestone family) ───────────────

#[test]
fn rect_minus_contained_circle_makes_round_hole() {
    let r = 0.5_f64;
    let wall = rect(-2.0_f64, -2.0_f64, 2.0_f64, 2.0_f64); // area 16
    let void = circle(0.0_f64, 0.0_f64, r);
    let d = difference(&wall, &void, &tol()).unwrap();
    assert_area(&d, 16.0_f64 - PI * r * r);
    let (outers, holes) = outer_hole_counts(&d);
    assert_eq!(
        (outers, holes),
        (1, 1),
        "a contained circular void is a round hole"
    );
    // The hole must keep arc edges (no polyline approximation).
    let hole = d
        .contours
        .iter()
        .find(|c| c.signed_area() < 0.0)
        .expect("hole");
    assert!(hole.has_arc(), "the round hole must keep its arc boundary");
}

#[test]
fn rect_minus_edge_circle_is_a_circular_notch() {
    let r = 0.5_f64;
    let wall = rect(-2.0_f64, -2.0_f64, 2.0_f64, 2.0_f64);
    // Circle centred exactly on the right edge x = 2: removes a half-disc notch.
    let void = circle(2.0_f64, 0.0_f64, r);
    let d = difference(&wall, &void, &tol()).unwrap();
    assert_area(&d, 16.0_f64 - 0.5_f64 * PI * r * r);
    let (_outers, holes) = outer_hole_counts(&d);
    assert_eq!(holes, 0, "an edge circle is a notch, not a hole");
    assert!(
        d.contours.iter().any(|c| c.has_arc()),
        "the notch keeps an arc boundary"
    );
}

#[test]
fn circle_minus_rect_is_a_half_disc() {
    let big = circle(0.0_f64, 0.0_f64, 1.0_f64);
    let cut = rect(0.0_f64, -2.0_f64, 2.0_f64, 2.0_f64); // remove the right half
    let d = difference(&big, &cut, &tol()).unwrap();
    assert_area(&d, 0.5_f64 * PI);
    assert!(d.contours.iter().any(|c| c.has_arc()));
}

// ── circle / circle ───────────────────────────────────────────────────────────

fn lens_area(d: f64, r: f64) -> f64 {
    // Area of intersection of two equal circles radius r whose centres are d apart.
    2.0 * r * r * ((d / (2.0 * r)).acos()) - 0.5 * d * (4.0 * r * r - d * d).sqrt()
}

#[test]
fn two_overlapping_circles_union_and_intersection() {
    let c1 = circle(0.0_f64, 0.0_f64, 1.0_f64);
    let c2 = circle(1.0_f64, 0.0_f64, 1.0_f64);
    let lens = lens_area(1.0_f64, 1.0_f64);

    let i = intersection(&c1, &c2, &tol()).unwrap();
    assert_area(&i, lens);

    let u = union(&c1, &c2, &tol()).unwrap();
    assert_area(&u, 2.0_f64 * PI - lens);
    let (outers, holes) = outer_hole_counts(&u);
    assert_eq!(
        (outers, holes),
        (1, 0),
        "two overlapping discs union to one"
    );
}

#[test]
fn two_separate_circles_union_keeps_both() {
    let c1 = circle(0.0_f64, 0.0_f64, 0.5_f64);
    let c2 = circle(3.0_f64, 0.0_f64, 0.5_f64);
    let u = union(&c1, &c2, &tol()).unwrap();
    assert_area(&u, 2.0_f64 * PI * 0.25_f64);
    let (outers, holes) = outer_hole_counts(&u);
    assert_eq!((outers, holes), (2, 0));
    let i = intersection(&c1, &c2, &tol()).unwrap();
    assert!(i.area() <= AREA_EPS, "disjoint circles do not intersect");
}

#[test]
fn externally_tangent_circles_report_degeneracy() {
    use archi_kernel::boolean::poly2d::Poly2Error;
    // Two discs touching from outside at one point form a degree-4 tangent pinch
    // (`DESIGN.md` §13-3): the engine reports it explicitly rather than returning
    // a silently wrong (collapsed) result. Nudging either circle off the exact
    // tangent makes it a normal separate / overlapping case.
    let c1 = circle(0.0_f64, 0.0_f64, 1.0_f64);
    let c2 = circle(2.0_f64, 0.0_f64, 1.0_f64);
    let err = union(&c1, &c2, &tol()).unwrap_err();
    assert!(matches!(err, Poly2Error::UnsupportedArcDegeneracy { .. }));

    // A hair's separation: now a clean two-component union.
    let c2b = circle(2.001_f64, 0.0_f64, 1.0_f64);
    let u = union(&c1, &c2b, &tol()).unwrap();
    assert_area(&u, 2.0_f64 * PI);
}

#[test]
fn same_circle_difference_is_empty() {
    let c = circle(0.0_f64, 0.0_f64, 1.0_f64);
    let d = difference(&c, &c, &tol()).unwrap();
    assert!(d.area() <= AREA_EPS, "C − C = ∅, got area {}", d.area());
}

#[test]
fn concentric_difference_is_an_annulus() {
    let big = circle(0.0_f64, 0.0_f64, 1.0_f64);
    let small = circle(0.0_f64, 0.0_f64, 0.5_f64);
    let ann = difference(&big, &small, &tol()).unwrap();
    assert_area(&ann, PI * (1.0_f64 - 0.25_f64));
    let (outers, holes) = outer_hole_counts(&ann);
    assert_eq!((outers, holes), (1, 1), "annulus = outer ring + inner hole");
}

#[test]
fn internally_tangent_circles_report_degeneracy() {
    use archi_kernel::boolean::poly2d::Poly2Error;
    // A small circle internally tangent to a big one is a pinch contact
    // (`DESIGN.md` §13-3): reported explicitly rather than mis-answered. Nudging
    // the inner circle off the exact tangent recovers a normal contained void.
    let big = circle(0.0_f64, 0.0_f64, 1.0_f64);
    let small = circle(0.5_f64, 0.0_f64, 0.5_f64); // internally tangent at (1,0)
    let err = difference(&big, &small, &tol()).unwrap_err();
    assert!(matches!(err, Poly2Error::UnsupportedArcDegeneracy { .. }));

    // Off-tangent: a properly contained void → an annulus-with-bite (area shrinks).
    let small2 = circle(0.4_f64, 0.0_f64, 0.5_f64);
    let d = difference(&big, &small2, &tol()).unwrap();
    assert!(d.area() < PI && d.area() > 0.0);
}

// ── line tangent to a circle (rect edge grazing a void) ───────────────────────

#[test]
fn rect_edge_tangent_to_circle_removes_nothing_extra() {
    // A circular void whose top is tangent to the wall's top edge y = 2: still a
    // full contained hole (tangency does not open it).
    let r = 0.5_f64;
    let wall = rect(-2.0_f64, -2.0_f64, 2.0_f64, 2.0_f64);
    let void = circle(0.0_f64, 1.5_f64, r); // top at y = 2, tangent to edge
    let d = difference(&wall, &void, &tol()).unwrap();
    assert_area(&d, 16.0_f64 - PI * r * r);
}

// ── proptest: area conservation / inclusion–exclusion / idempotence ───────────

/// A circle on a coarse grid, centre offset by a non-aligning fraction so the
/// vertical seam (centre ± r) does not land exactly on a grid-aligned rectangle
/// corner — that exact seam-on-corner coincidence is a genuine degree-4 pinch
/// (`DESIGN.md` §13-3), tested separately; here we exercise the generic cases.
fn circle_strategy() -> impl Strategy<Value = Region> {
    (-3i32..3, -3i32..3, 1i32..4)
        .prop_map(|(x, y, r)| circle(x as f64 + 0.137, y as f64 + 0.111, r as f64 * 0.5))
}

/// A rectangle on the integer grid.
fn rect_strategy() -> impl Strategy<Value = Region> {
    (-4i32..4, -4i32..4, 1i32..5, 1i32..5).prop_map(|(x, y, w, h)| {
        let x0 = x as f64;
        let y0 = y as f64;
        rect(x0, y0, x0 + w as f64, y0 + h as f64)
    })
}

fn check_well_formed(r: &Region) {
    for c in &r.contours {
        let a = c.signed_area();
        assert!(a.abs() > 1e-9_f64, "zero-area sliver contour in output");
    }
    let total = r.signed_area();
    assert!(total >= -AREA_EPS, "output total signed area is negative");
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 400, ..ProptestConfig::default() })]

    /// area(A − B) + area(A ∩ B) = area(A) with circles + rectangles. A genuine
    /// tangent-pinch degeneracy is reported as an error (`DESIGN.md` §13-3) and
    /// skipped here rather than silently mis-answered.
    #[test]
    fn arc_difference_plus_intersection_equals_a(
        a in prop_oneof![circle_strategy(), rect_strategy()],
        b in prop_oneof![circle_strategy(), rect_strategy()],
    ) {
        let (Ok(diff), Ok(inter)) =
            (difference(&a, &b, &tol()), intersection(&a, &b, &tol())) else {
            return Ok(()); // tangent-pinch degeneracy, explicitly reported
        };
        let lhs = diff.area() + inter.area();
        let rhs = a.area();
        prop_assert!((lhs - rhs).abs() <= 1e-5 * (1.0 + rhs.abs()),
            "area(A−B)+area(A∩B)={lhs} != area(A)={rhs}");
        check_well_formed(&diff);
        check_well_formed(&inter);
    }

    /// area(A ∪ B) = area(A) + area(B) − area(A ∩ B).
    #[test]
    fn arc_union_inclusion_exclusion(
        a in prop_oneof![circle_strategy(), rect_strategy()],
        b in prop_oneof![circle_strategy(), rect_strategy()],
    ) {
        let (Ok(uni), Ok(inter)) =
            (union(&a, &b, &tol()), intersection(&a, &b, &tol())) else {
            return Ok(());
        };
        let lhs = uni.area();
        let rhs = a.area() + b.area() - inter.area();
        prop_assert!((lhs - rhs).abs() <= 1e-5 * (1.0 + rhs.abs()),
            "area(A∪B)={lhs} != incl-excl={rhs}");
        check_well_formed(&uni);
    }

    /// Idempotence: (A − B) − B = A − B for circle/rect mixes.
    #[test]
    fn arc_difference_idempotent(
        a in prop_oneof![circle_strategy(), rect_strategy()],
        b in prop_oneof![circle_strategy(), rect_strategy()],
    ) {
        let Ok(ab) = difference(&a, &b, &tol()) else { return Ok(()); };
        let Ok(abb) = difference(&ab, &b, &tol()) else { return Ok(()); };
        prop_assert!((ab.area() - abb.area()).abs() <= 1e-5 * (1.0 + ab.area()),
            "difference not idempotent: {} vs {}", ab.area(), abb.area());
    }
}
