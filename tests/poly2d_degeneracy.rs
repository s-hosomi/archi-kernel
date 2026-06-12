//! Degeneracy suite — the building-domain "exact coincidence is the common
//! case" scenarios. If any of these regress after integration, the engine is
//! broken in the way that matters most for archi-kernel.

use archi_kernel::boolean::poly2d::{difference, intersection, union, Point2, Poly2Error, Region};
use archi_kernel::tolerance::Tol;

fn tol() -> Tol {
    Tol::default()
}

/// Build an axis-aligned rectangle region (CCW) from corner `(x0,y0)` to `(x1,y1)`.
fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Region {
    Region::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

/// Build a region from an explicit CCW point list (single contour).
fn poly(points: &[(f64, f64)]) -> Region {
    use archi_kernel::boolean::poly2d::Contour;
    let pts: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    Region::new(vec![Contour::from_points(&pts)])
}

/// Assert two areas agree to `1e-7`.
fn assert_area(region: &Region, expected: f64) {
    let got = region.area();
    assert!(
        (got - expected).abs() <= 1e-7_f64,
        "area mismatch: got {got}, expected {expected}"
    );
}

/// Validate structural invariants every output region must satisfy.
fn assert_well_formed(region: &Region) {
    for (i, c) in region.contours.iter().enumerate() {
        let a = c.signed_area();
        assert!(
            a.abs() > 1e-12_f64,
            "contour {i} is a zero-area sliver (signed area {a})"
        );
    }
    let total = region.signed_area();
    assert!(
        total >= -1e-7_f64,
        "total signed area is negative ({total}); orientation is inconsistent"
    );
}

/// Count contours with positive (outer) and negative (hole) signed area.
fn outer_hole_counts(region: &Region) -> (usize, usize) {
    let mut outers = 0usize;
    let mut holes = 0usize;
    for c in &region.contours {
        if c.signed_area() > 0.0_f64 {
            outers += 1;
        } else if c.signed_area() < 0.0_f64 {
            holes += 1;
        }
    }
    (outers, holes)
}

#[test]
fn identical_square_difference_is_empty() {
    // A − A = ∅.
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let r = difference(&a, &a, &tol()).unwrap();
    assert!(r.is_empty(), "A − A must be empty, got area {}", r.area());
}

#[test]
fn identical_square_intersection_is_self() {
    // A ∩ A = A.
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let r = intersection(&a, &a, &tol()).unwrap();
    assert_area(&r, 1.0_f64);
    assert_well_formed(&r);
}

#[test]
fn identical_square_union_is_self() {
    // A ∪ A = A (idempotent), and the result is a single square (no doubled
    // boundary).
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let r = union(&a, &a, &tol()).unwrap();
    assert_area(&r, 1.0_f64);
    assert_well_formed(&r);
    let (outers, holes) = outer_hole_counts(&r);
    assert_eq!((outers, holes), (1, 0));
}

#[test]
fn shared_edge_union_is_one_rectangle() {
    // Two unit squares sharing the edge x=1 → a single 2×1 rectangle; the
    // internal shared edge must vanish (one outer contour, no holes).
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let b = rect(1.0_f64, 0.0_f64, 2.0_f64, 1.0_f64);
    let r = union(&a, &b, &tol()).unwrap();
    assert_area(&r, 2.0_f64);
    assert_well_formed(&r);
    let (outers, holes) = outer_hole_counts(&r);
    assert_eq!(
        (outers, holes),
        (1, 0),
        "shared edge should merge into one rectangle"
    );
}

#[test]
fn shared_edge_difference_is_unchanged() {
    // Two squares touching only along the edge x=1. A − B removes nothing (they
    // share no area), so the result is A.
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let b = rect(1.0_f64, 0.0_f64, 2.0_f64, 1.0_f64);
    let r = difference(&a, &b, &tol()).unwrap();
    assert_area(&r, 1.0_f64);
    assert_well_formed(&r);
}

#[test]
fn shared_edge_intersection_is_empty_or_zero_area() {
    // Edge-only contact has zero area intersection → regularized result is empty
    // (no dangling 1-D sliver).
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let b = rect(1.0_f64, 0.0_f64, 2.0_f64, 1.0_f64);
    let r = intersection(&a, &b, &tol()).unwrap();
    assert!(
        r.area() <= 1e-9_f64,
        "edge contact has no interior; got area {}",
        r.area()
    );
}

#[test]
fn corner_touch_union_keeps_both() {
    // Two squares touching only at the corner (1,1). Union keeps both; area 2.
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let b = rect(1.0_f64, 1.0_f64, 2.0_f64, 2.0_f64);
    let r = union(&a, &b, &tol()).unwrap();
    assert_area(&r, 2.0_f64);
    assert_well_formed(&r);
}

#[test]
fn corner_touch_intersection_empty() {
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let b = rect(1.0_f64, 1.0_f64, 2.0_f64, 2.0_f64);
    let r = intersection(&a, &b, &tol()).unwrap();
    assert!(r.area() <= 1e-9_f64);
}

#[test]
fn vertex_on_edge_difference() {
    // B's left edge midpoint touches A's interior: B = [1,3]×[0,1] partly inside
    // A = [0,2]×[0,1], and a third square C starts a vertex exactly on A's edge.
    // Here: A=[0,2]×[0,1], B=[1,2]×[0.25,0.75] is a notch fully inside A's right
    // half whose corners (2,0.25),(2,0.75) land ON A's right edge x=2.
    let a = rect(0.0_f64, 0.0_f64, 2.0_f64, 1.0_f64);
    let b = rect(1.0_f64, 0.25_f64, 2.0_f64, 0.75_f64);
    let r = difference(&a, &b, &tol()).unwrap();
    // Area = 2 − 0.5 = 1.5.
    assert_area(&r, 1.5_f64);
    assert_well_formed(&r);
}

#[test]
fn opening_flush_with_outer_edge() {
    // An opening whose right edge is flush with the wall's right edge: the
    // result is a C-shape (open on the right), not a donut. No hole, one outer.
    let wall = rect(0.0_f64, 0.0_f64, 4.0_f64, 3.0_f64);
    let opening = rect(2.0_f64, 1.0_f64, 4.0_f64, 2.0_f64); // flush with x=4
    let r = difference(&wall, &opening, &tol()).unwrap();
    assert_area(&r, 12.0_f64 - 2.0_f64);
    assert_well_formed(&r);
    let (outers, holes) = outer_hole_counts(&r);
    assert_eq!(holes, 0, "flush opening should not create a hole");
    assert_eq!(outers, 1);
}

#[test]
fn interior_opening_makes_a_hole() {
    // An opening strictly inside the wall makes a donut: one outer, one hole.
    let wall = rect(0.0_f64, 0.0_f64, 4.0_f64, 3.0_f64);
    let opening = rect(1.0_f64, 1.0_f64, 2.0_f64, 2.0_f64);
    let r = difference(&wall, &opening, &tol()).unwrap();
    assert_area(&r, 12.0_f64 - 1.0_f64);
    assert_well_formed(&r);
    let (outers, holes) = outer_hole_counts(&r);
    assert_eq!(
        (outers, holes),
        (1, 1),
        "interior opening must create one hole"
    );
}

#[test]
fn hole_touching_outer_edge_from_inside() {
    // Opening touching the outer boundary at a single point (its corner on the
    // wall edge). Should remain a valid region with the expected area.
    let wall = rect(0.0_f64, 0.0_f64, 4.0_f64, 3.0_f64);
    let opening = rect(0.0_f64, 1.0_f64, 1.0_f64, 2.0_f64); // left edge flush x=0
    let r = difference(&wall, &opening, &tol()).unwrap();
    assert_area(&r, 12.0_f64 - 1.0_f64);
    assert_well_formed(&r);
}

#[test]
fn nested_holes_island_in_hole() {
    // Wall with a big opening, and inside the opening a smaller solid island
    // (achieved by unioning the wall-with-hole result with an inner square).
    let wall = rect(0.0_f64, 0.0_f64, 6.0_f64, 6.0_f64);
    let big_opening = rect(1.0_f64, 1.0_f64, 5.0_f64, 5.0_f64);
    let donut = difference(&wall, &big_opening, &tol()).unwrap();
    let island = rect(2.0_f64, 2.0_f64, 4.0_f64, 4.0_f64);
    let r = union(&donut, &island, &tol()).unwrap();
    // Area = (36 − 16) + 4 = 24.
    assert_area(&r, 24.0_f64);
    assert_well_formed(&r);
    let (outers, holes) = outer_hole_counts(&r);
    // Two outer rings (wall frame + island) and one hole (the ring gap).
    assert_eq!(outers, 2, "expected wall frame outer + island outer");
    assert_eq!(holes, 1, "expected the ring-shaped gap as a hole");
}

#[test]
fn l_shape_minus_rect() {
    // An L-shape minus a rectangle that clips its inner corner.
    let l = poly(&[
        (0.0_f64, 0.0_f64),
        (3.0_f64, 0.0_f64),
        (3.0_f64, 1.0_f64),
        (1.0_f64, 1.0_f64),
        (1.0_f64, 3.0_f64),
        (0.0_f64, 3.0_f64),
    ]); // area = 3*1 + 1*2 = 5
    let cut = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64); // removes the corner unit square
    let r = difference(&l, &cut, &tol()).unwrap();
    assert_area(&r, 4.0_f64);
    assert_well_formed(&r);
}

#[test]
fn zero_width_bridge_does_not_spawn_sliver() {
    // Two squares connected by a coincident shared edge segment of partial
    // length (the classic "zero-width bridge" where boundaries overlap on a
    // sub-segment). Union must produce a clean single region with no zero-area
    // fragment.
    let a = rect(0.0_f64, 0.0_f64, 2.0_f64, 2.0_f64);
    let b = rect(2.0_f64, 0.5_f64, 4.0_f64, 1.5_f64); // shares part of edge x=2
    let r = union(&a, &b, &tol()).unwrap();
    assert_area(&r, 4.0_f64 + 2.0_f64);
    assert_well_formed(&r); // assert_well_formed rejects sliver contours
}

#[test]
fn fully_contained_difference_makes_hole() {
    let a = rect(0.0_f64, 0.0_f64, 10.0_f64, 10.0_f64);
    let b = rect(3.0_f64, 3.0_f64, 7.0_f64, 7.0_f64);
    let r = difference(&a, &b, &tol()).unwrap();
    assert_area(&r, 100.0_f64 - 16.0_f64);
    let (outers, holes) = outer_hole_counts(&r);
    assert_eq!((outers, holes), (1, 1));
    assert_well_formed(&r);
}

#[test]
fn disjoint_union_two_components() {
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let b = rect(5.0_f64, 5.0_f64, 6.0_f64, 6.0_f64);
    let r = union(&a, &b, &tol()).unwrap();
    assert_area(&r, 2.0_f64);
    let (outers, holes) = outer_hole_counts(&r);
    assert_eq!((outers, holes), (2, 0));
}

#[test]
fn disjoint_difference_unchanged() {
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let b = rect(5.0_f64, 5.0_f64, 6.0_f64, 6.0_f64);
    let r = difference(&a, &b, &tol()).unwrap();
    assert_area(&r, 1.0_f64);
}

#[test]
fn b_fully_covers_a_difference_empty() {
    let a = rect(2.0_f64, 2.0_f64, 3.0_f64, 3.0_f64);
    let b = rect(0.0_f64, 0.0_f64, 5.0_f64, 5.0_f64);
    let r = difference(&a, &b, &tol()).unwrap();
    assert!(
        r.area() <= 1e-9_f64,
        "A ⊂ B ⇒ A − B = ∅, got area {}",
        r.area()
    );
}

#[test]
fn near_coincident_within_eps_treated_as_coincident() {
    // B is offset from A by 1e-7 (< eps = 1e-6); they should be treated as the
    // identical square, so A − B ≈ ∅.
    let a = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let b = Region::from_points(&[
        Point2::new(1e-7_f64, 0.0_f64),
        Point2::new(1.0_f64 + 1e-7_f64, 0.0_f64),
        Point2::new(1.0_f64 + 1e-7_f64, 1.0_f64),
        Point2::new(1e-7_f64, 1.0_f64),
    ]);
    let r = difference(&a, &b, &tol()).unwrap();
    assert!(
        r.area() <= 1e-5_f64,
        "near-coincident squares: A−B should be ~empty, got {}",
        r.area()
    );
}

#[test]
fn arc_input_is_rejected() {
    use archi_kernel::boolean::poly2d::{Arc, Contour, Edge2};
    let arc_contour = Contour::new(vec![
        Edge2::seg(Point2::new(0.0_f64, 0.0_f64), Point2::new(1.0_f64, 0.0_f64)),
        Edge2::Arc(Arc::new(
            Point2::new(0.5_f64, 0.0_f64),
            0.5_f64,
            0.0_f64,
            std::f64::consts::PI,
        )),
        Edge2::seg(Point2::new(0.0_f64, 0.0_f64), Point2::new(0.0_f64, 0.0_f64)),
    ]);
    let a = Region::new(vec![arc_contour]);
    let b = rect(0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64);
    let err = difference(&a, &b, &tol()).unwrap_err();
    assert!(matches!(err, Poly2Error::ArcNotYetSupported));
}
