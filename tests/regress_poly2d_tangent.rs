//! Regression spec for the seg×arc tangent split-drop bug.
//!
//! ## Symptom
//! When a circular arc is tangent to one or more polygon edges and the contact
//! points are shared vertices, boolean ops silently drop a whole bounded face
//! and the result becomes argument-order dependent:
//!   * `intersection(circle, square)` = 0   but
//!   * `intersection(square, circle)` = π/4  (the correct quarter disc)
//!
//! and `V(C∩S)+V(C−S)` reads 2.356 instead of π. The architectural form
//! (a 4×2 wall minus a r=1 void tangent to the top/bottom edges) collapses its
//! difference to 0 instead of 8−π.
//!
//! ## Root cause (the DCEL tangent-pinch trace, NOT `intersect.rs`)
//! The tangent split is *already present* in the arrangement regardless of
//! `intersect.rs`, because (a) the circle's seam vertices already sit exactly on
//! the tangent points (`Region::circle` seams a circle at its top/bottom), and
//! (b) the segment-only vertex-on-edge snapping pass in `arrangement.rs`
//! (`build()`, the "Vertex-on-edge snapping for straight edges" loop) splits the
//! polygon edge at any existing vertex lying on its interior. The face was lost
//! later, in the **DCEL ring ordering** at the degree-4 tangent pinch: at the
//! shared vertex an arc departs with a tangent direction collinear with the
//! segment (and, in the wall/void case, with the opposite arc), so ordering the
//! outgoing half-edges by tangent angle alone left two half-edges tied. The tie
//! was resolved arbitrarily by the input arrival order, so `next`/`prev` linked
//! the wrong neighbour and `trace_faces` dropped the bounded face — which is why
//! the result was argument-order dependent.
//!
//! The fix (`arrangement.rs`, `leave_key` / `LeaveKey::cmp_ring`) adds a **signed
//! curvature** secondary key to the ring sort: when two half-edges leave a vertex
//! with the same tangent direction, the arc that bends toward its centre is
//! ordered to the correct side of the straight (or less-curved) edge, giving the
//! tangent pinch a deterministic, geometry-correct ring order.
//!
//! These tests encode the correct behaviour. The control case (which never had a
//! tangent contact) must stay green so the fix does not regress the non-tangent
//! path.

use archi_kernel::boolean::poly2d::{difference, intersection, union, Point2, Region};
use archi_kernel::tolerance::Tol;
use std::f64::consts::PI;

fn unit_square() -> Region {
    Region::from_points(&[
        Point2::new(0.0_f64, 0.0_f64),
        Point2::new(1.0_f64, 0.0_f64),
        Point2::new(1.0_f64, 1.0_f64),
        Point2::new(0.0_f64, 1.0_f64),
    ])
}

fn wall_4x2() -> Region {
    Region::from_points(&[
        Point2::new(0.0_f64, 0.0_f64),
        Point2::new(4.0_f64, 0.0_f64),
        Point2::new(4.0_f64, 2.0_f64),
        Point2::new(0.0_f64, 2.0_f64),
    ])
}

/// The intersection of the unit circle and the unit square is the quarter disc
/// in the first quadrant, area π/4, **independent of argument order**.
#[test]
fn circle_square_corner_tangent_intersection() {
    let tol = Tol::default();
    let circle = Region::circle(Point2::new(0.0_f64, 0.0_f64), 1.0_f64);
    let square = unit_square();

    let cs_area = intersection(&circle, &square, &tol).unwrap().area();
    let sc_area = intersection(&square, &circle, &tol).unwrap().area();
    let expected = PI / 4.0;
    assert!(
        (cs_area - expected).abs() <= 1e-6_f64,
        "intersection(circle, square) = {cs_area}, expected {expected}"
    );
    assert!(
        (sc_area - expected).abs() <= 1e-6_f64,
        "intersection(square, circle) = {sc_area}, expected {expected}"
    );
    assert!(
        (cs_area - sc_area).abs() <= 1e-9_f64,
        "intersection must be order-independent: {cs_area} vs {sc_area}"
    );
}

/// Volume identity V(C∩S) + V(C−S) = V(C) = π must hold.
#[test]
fn circle_square_volume_identity() {
    let tol = Tol::default();
    let circle = Region::circle(Point2::new(0.0_f64, 0.0_f64), 1.0_f64);
    let square = unit_square();

    let inter = intersection(&circle, &square, &tol).unwrap().area();
    let diff = difference(&circle, &square, &tol).unwrap().area();
    let sum = inter + diff;
    assert!(
        (sum - PI).abs() <= 1e-6_f64,
        "V(C∩S)+V(C−S) = {sum}, expected π = {PI}"
    );
}

/// Union sanity: V(C∪S) = V(C) + V(S) − V(C∩S) = π + 1 − π/4.
#[test]
fn circle_square_union_consistent() {
    let tol = Tol::default();
    let circle = Region::circle(Point2::new(0.0_f64, 0.0_f64), 1.0_f64);
    let square = unit_square();

    let u = union(&circle, &square, &tol).unwrap().area();
    let expected = PI + 1.0 - PI / 4.0;
    assert!(
        (u - expected).abs() <= 1e-6_f64,
        "V(C∪S) = {u}, expected {expected}"
    );
}

/// Architectural case: a 4×2 wall minus a r=1 circular void centred at (2,1)
/// that is tangent to the top/bottom edges at (2,2)/(2,0). The difference is
/// 8 − π; the bug collapses it to 0.
#[test]
fn wall_void_tangent_difference() {
    let tol = Tol::default();
    let wall = wall_4x2();
    let void = Region::circle(Point2::new(2.0_f64, 1.0_f64), 1.0_f64);
    let area = difference(&wall, &void, &tol).unwrap().area();
    let expected = 8.0 - PI;
    assert!(
        (area - expected).abs() <= 1e-6_f64,
        "wall − void = {area}, expected 8 − π = {expected}"
    );
}

/// Control: the same wall minus a r=0.9 void (pulled off the edges, no tangent
/// contact) is already correct and must stay correct. This active test guards
/// the non-tangent path so any future arrangement-side fix does not regress it.
#[test]
fn wall_void_pulled_off_edges_ok_control() {
    let tol = Tol::default();
    let wall = wall_4x2();
    let void = Region::circle(Point2::new(2.0_f64, 1.0_f64), 0.9_f64);
    let area = difference(&wall, &void, &tol).unwrap().area();
    let expected = 8.0 - PI * 0.9 * 0.9;
    assert!(
        (area - expected).abs() <= 1e-6_f64,
        "control = {area}, expected {expected}"
    );
}
