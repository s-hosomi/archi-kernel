//! Pinned regression for the sub-Tol jitter boolean identity defect.
//!
//! Found by `poly2d_tol_jitter_identity_proptest` (tests/adversarial_boolean.rs):
//! a 0.5×0.5 square `b` overlapping the corner `(0,0)` of the unit square `a`
//! with a *sub-Tol* offset on both axes. Snap-merge collapses `b`'s near corner
//! onto `a`'s corner, after which the classification of `b`'s edges relative to
//! `a`'s edges broke and the `0.25` intersection cell silently vanished, so
//! `area(A−B) + area(A∩B) = 0.7499995… ≠ 1.0`.

use archi_kernel::boolean::poly2d::{self, Point2, Region};
use archi_kernel::tolerance::Tol;

fn region_rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Region {
    Region::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

/// Minimal failing input from the proptest, pinned. Both the difference and the
/// intersection must be individually correct, and the identity must hold to
/// regularization precision.
#[test]
fn poly2d_jitter_corner_overlap_identity() {
    let tol = Tol::default();
    // sub-Tol jitter (Tol::length = 1e-6): b's corner is within eps of a's (0,0)
    let jx = 2.536902937174801e-7_f64;
    let jy = 4.773365684494552e-7_f64;
    let a = region_rect(0.0, 0.0, 1.0, 1.0);
    let b = region_rect(jx, jy, jx + 0.5, jy + 0.5);

    let diff = poly2d::difference(&a, &b, &tol).expect("a-b");
    let inter = poly2d::intersection(&a, &b, &tol).expect("a∩b");

    let da = diff.area();
    let ia = inter.area();
    assert!((ia - 0.25).abs() <= 1e-4, "A∩B area = {ia}, expected ≈0.25");
    assert!((da - 0.75).abs() <= 1e-4, "A−B area = {da}, expected ≈0.75");
    assert!(
        (da + ia - 1.0).abs() <= 1e-4,
        "identity {} vs 1.0 (A−B={da}, A∩B={ia})",
        da + ia
    );
}

/// Same input repeated 32× must give identical areas every call — guards against
/// `HashMap`-iteration-order nondeterminism in the snap / arrangement layers.
#[test]
fn poly2d_jitter_corner_overlap_deterministic() {
    let tol = Tol::default();
    let jx = 2.536902937174801e-7_f64;
    let jy = 4.773365684494552e-7_f64;
    let a = region_rect(0.0, 0.0, 1.0, 1.0);
    let b = region_rect(jx, jy, jx + 0.5, jy + 0.5);

    let mut diffs = Vec::with_capacity(32);
    let mut inters = Vec::with_capacity(32);
    for _ in 0..32 {
        let diff = poly2d::difference(&a, &b, &tol).expect("a-b");
        let inter = poly2d::intersection(&a, &b, &tol).expect("a∩b");
        diffs.push(diff.area());
        inters.push(inter.area());
    }
    let d0 = diffs[0];
    let i0 = inters[0];
    assert!(
        diffs.iter().all(|&x| (x - d0).abs() <= 1e-12_f64),
        "non-deterministic A−B areas: {diffs:?}"
    );
    assert!(
        inters.iter().all(|&x| (x - i0).abs() <= 1e-12_f64),
        "non-deterministic A∩B areas: {inters:?}"
    );
    assert!((d0 - 0.75).abs() <= 1e-4, "A−B = {d0}, expected 0.75");
    assert!((i0 - 0.25).abs() <= 1e-4, "A∩B = {i0}, expected 0.25");
}
