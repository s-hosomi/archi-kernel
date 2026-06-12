//! Property tests — the statistical robustness net.
//!
//! These run hundreds of cases each (>1000 total) and check the algebraic
//! identities that must hold for any correct boolean engine, plus structural
//! validity of every output. Coordinates are drawn on a coarse integer-ish grid
//! so that "exact coincidence" (the building-domain common case) is hit *often*,
//! not avoided.

use archi_kernel::boolean::poly2d::{difference, intersection, union, Contour, Point2, Region};
use archi_kernel::tolerance::Tol;
use proptest::prelude::*;

fn tol() -> Tol {
    Tol::default()
}

/// Area-comparison tolerance: generous relative to coordinate magnitudes (grid
/// up to ~20), accounting for accumulated floating-point error.
const AREA_EPS: f64 = 1e-6;

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
    let pts: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    Region::new(vec![Contour::from_points(&pts)])
}

/// A strategy for axis-aligned rectangles on a coarse grid (lots of shared
/// coordinates → frequent coincidence).
fn rect_strategy() -> impl Strategy<Value = Region> {
    (0i32..10, 0i32..10, 1i32..6, 1i32..6).prop_map(|(x, y, w, h)| {
        let x0 = x as f64;
        let y0 = y as f64;
        rect(x0, y0, x0 + w as f64, y0 + h as f64)
    })
}

/// A strategy for axis-aligned L / rectangle / step shapes, still grid-aligned.
fn shape_strategy() -> impl Strategy<Value = Region> {
    prop_oneof![
        rect_strategy(),
        // An L-shape from two grid params.
        (0i32..8, 0i32..8, 2i32..6, 2i32..6).prop_map(|(x, y, a, b)| {
            let x0 = x as f64;
            let y0 = y as f64;
            let a = a as f64;
            let b = b as f64;
            poly(&[
                (x0, y0),
                (x0 + a, y0),
                (x0 + a, y0 + b / 2.0),
                (x0 + a / 2.0, y0 + b / 2.0),
                (x0 + a / 2.0, y0 + b),
                (x0, y0 + b),
            ])
        }),
    ]
}

/// Assert the structural invariants of an output region:
/// * no zero-area sliver contour,
/// * each contour is simple (no self-intersection),
/// * total signed area non-negative (outer CCW, holes CW; outers ≥ holes).
fn check_well_formed(r: &Region) {
    for c in &r.contours {
        let a = c.signed_area();
        prop_assert_ok(a.abs() > 1e-9_f64, "zero-area sliver contour in output");
        prop_assert_ok(contour_is_simple(c), "self-intersecting contour in output");
    }
    let total = r.signed_area();
    prop_assert_ok(total >= -AREA_EPS, "output total signed area is negative");
}

/// `true` if a contour's segments do not cross except at shared endpoints.
fn contour_is_simple(c: &Contour) -> bool {
    let v = c.vertices();
    let n = v.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        let a0 = v[i];
        let a1 = v[(i + 1) % n];
        for j in (i + 1)..n {
            // Skip adjacent edges (they legitimately share an endpoint).
            if j == i || (j + 1) % n == i || (i + 1) % n == j {
                continue;
            }
            let b0 = v[j];
            let b1 = v[(j + 1) % n];
            if segments_properly_cross(a0, a1, b0, b1) {
                return false;
            }
        }
    }
    true
}

/// Proper crossing test (interiors intersect) using orientation signs.
fn segments_properly_cross(a0: Point2, a1: Point2, b0: Point2, b1: Point2) -> bool {
    fn cross(o: Point2, a: Point2, b: Point2) -> f64 {
        (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
    }
    let d1 = cross(a0, a1, b0);
    let d2 = cross(a0, a1, b1);
    let d3 = cross(b0, b1, a0);
    let d4 = cross(b0, b1, a1);
    // Strictly opposite signs on both segments ⇒ a proper interior crossing.
    ((d1 > 0.0) != (d2 > 0.0))
        && ((d3 > 0.0) != (d4 > 0.0))
        && d1 != 0.0
        && d2 != 0.0
        && d3 != 0.0
        && d4 != 0.0
}

/// Tiny helper so `check_well_formed` can be called from inside `proptest!`
/// closures (which need `prop_assert!`-style early return). We instead panic
/// with a message, which proptest captures and shrinks.
fn prop_assert_ok(cond: bool, msg: &str) {
    assert!(cond, "{msg}");
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 1200, ..ProptestConfig::default() })]

    /// area(A − B) + area(A ∩ B) = area(A).
    #[test]
    fn difference_plus_intersection_equals_a(a in shape_strategy(), b in shape_strategy()) {
        let diff = difference(&a, &b, &tol()).unwrap();
        let inter = intersection(&a, &b, &tol()).unwrap();
        let lhs = diff.area() + inter.area();
        let rhs = a.area();
        prop_assert!((lhs - rhs).abs() <= AREA_EPS * (1.0 + rhs.abs()),
            "area(A−B)+area(A∩B)={lhs} != area(A)={rhs}");
        check_well_formed(&diff);
        check_well_formed(&inter);
    }

    /// area(A ∪ B) = area(A) + area(B) − area(A ∩ B).
    #[test]
    fn union_inclusion_exclusion(a in shape_strategy(), b in shape_strategy()) {
        let uni = union(&a, &b, &tol()).unwrap();
        let inter = intersection(&a, &b, &tol()).unwrap();
        let lhs = uni.area();
        let rhs = a.area() + b.area() - inter.area();
        prop_assert!((lhs - rhs).abs() <= AREA_EPS * (1.0 + rhs.abs()),
            "area(A∪B)={lhs} != area(A)+area(B)−area(A∩B)={rhs}");
        check_well_formed(&uni);
    }

    /// Idempotence of union: (A ∪ B) ∪ B = A ∪ B.
    #[test]
    fn union_idempotent(a in rect_strategy(), b in rect_strategy()) {
        let ab = union(&a, &b, &tol()).unwrap();
        let abb = union(&ab, &b, &tol()).unwrap();
        prop_assert!((ab.area() - abb.area()).abs() <= AREA_EPS * (1.0 + ab.area()),
            "union not idempotent: {} vs {}", ab.area(), abb.area());
        check_well_formed(&abb);
    }

    /// Idempotence of intersection: (A ∩ B) ∩ B = A ∩ B.
    #[test]
    fn intersection_idempotent(a in rect_strategy(), b in rect_strategy()) {
        let ab = intersection(&a, &b, &tol()).unwrap();
        let abb = intersection(&ab, &b, &tol()).unwrap();
        prop_assert!((ab.area() - abb.area()).abs() <= AREA_EPS * (1.0 + ab.area()),
            "intersection not idempotent: {} vs {}", ab.area(), abb.area());
    }

    /// Difference removes nothing twice: (A − B) − B = A − B.
    #[test]
    fn difference_idempotent(a in rect_strategy(), b in rect_strategy()) {
        let ab = difference(&a, &b, &tol()).unwrap();
        let abb = difference(&ab, &b, &tol()).unwrap();
        prop_assert!((ab.area() - abb.area()).abs() <= AREA_EPS * (1.0 + ab.area()),
            "difference not idempotent: {} vs {}", ab.area(), abb.area());
    }

    /// A − A = ∅ for any shape.
    #[test]
    fn self_difference_empty(a in shape_strategy()) {
        let r = difference(&a, &a, &tol()).unwrap();
        prop_assert!(r.area() <= AREA_EPS, "A − A not empty: area {}", r.area());
    }

    /// A ∩ A = A and A ∪ A = A.
    #[test]
    fn self_intersection_union(a in shape_strategy()) {
        let i = intersection(&a, &a, &tol()).unwrap();
        let u = union(&a, &a, &tol()).unwrap();
        prop_assert!((i.area() - a.area()).abs() <= AREA_EPS * (1.0 + a.area()));
        prop_assert!((u.area() - a.area()).abs() <= AREA_EPS * (1.0 + a.area()));
        check_well_formed(&i);
        check_well_formed(&u);
    }

    /// Commutativity of union and intersection (by area).
    #[test]
    fn union_intersection_commute(a in shape_strategy(), b in shape_strategy()) {
        let uab = union(&a, &b, &tol()).unwrap();
        let uba = union(&b, &a, &tol()).unwrap();
        let iab = intersection(&a, &b, &tol()).unwrap();
        let iba = intersection(&b, &a, &tol()).unwrap();
        prop_assert!((uab.area() - uba.area()).abs() <= AREA_EPS * (1.0 + uab.area()));
        prop_assert!((iab.area() - iba.area()).abs() <= AREA_EPS * (1.0 + iab.area()));
    }
}

/// Off-grid (rotated) rectangles to exercise the non-axis-aligned arrangement
/// path with genuine transversal crossings. Fewer cases but still meaningful.
fn rotated_rect(cx: f64, cy: f64, hw: f64, hh: f64, theta: f64) -> Region {
    let (s, c) = theta.sin_cos();
    let corners = [(-hw, -hh), (hw, -hh), (hw, hh), (-hw, hh)];
    let pts: Vec<Point2> = corners
        .iter()
        .map(|&(x, y)| Point2::new(cx + x * c - y * s, cy + x * s + y * c))
        .collect();
    Region::from_points(&pts)
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 1200, ..ProptestConfig::default() })]

    #[test]
    fn rotated_area_conservation(
        ax in -3.0f64..3.0, ay in -3.0f64..3.0, at in 0.0f64..1.5,
        bx in -3.0f64..3.0, by in -3.0f64..3.0, bt in 0.0f64..1.5,
    ) {
        let a = rotated_rect(ax, ay, 2.0, 1.0, at);
        let b = rotated_rect(bx, by, 1.5, 1.5, bt);
        let diff = difference(&a, &b, &tol()).unwrap();
        let inter = intersection(&a, &b, &tol()).unwrap();
        let lhs = diff.area() + inter.area();
        let rhs = a.area();
        // Rotated cases accumulate more error; scale tolerance with perimeter.
        let scale = 1e-5 * (1.0 + rhs.abs());
        prop_assert!((lhs - rhs).abs() <= scale,
            "rotated area(A−B)+area(A∩B)={lhs} != area(A)={rhs}");
        check_well_formed(&diff);
        check_well_formed(&inter);
    }
}
