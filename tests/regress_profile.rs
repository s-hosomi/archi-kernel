//! Regression tests for two profile / curve bugs:
//!
//! 1. **H-section degenerate polygon** (`profile/mod.rs`): `h_section_outline`
//!    used `web > b` (strictly greater), so `web == b` slipped through and
//!    produced a 12-vertex polygon with 4 zero-length edges.  The fix changes
//!    the guard to `web >= b`.
//!
//! 2. **`Ellipse3::new` missing orthogonality check** (`primitives/curve.rs`):
//!    when `major_dir == normal` the cross product `v = normal × major_dir`
//!    degenerates to zero, making `point_at` return incorrect values.  The fix
//!    adds a Tol-guarded perpendicularity check and Gram-Schmidt correction.

use archi_kernel::csg::Profile2d;
use archi_kernel::error::KernelError;
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::Ellipse3;
use archi_kernel::profile::ProfileGeom;
use archi_kernel::tolerance::Tol;

// ── Bug 1: H-section web == b produces degenerate polygon ────────────────────

/// Before the fix: `web == b` slipped through `if web > b`, yielding a valid
/// `Ok(Polygon)` with 4 pairs of duplicate adjacent vertices (zero-length edges).
/// After the fix: this must return `Err(NonPositiveDimension)`.
#[test]
fn h_section_rejects_web_equal_to_b() {
    let half_w = 0.1_f64;
    let half_h = 0.2_f64;
    let web = 0.2_f64; // == b = 2 * half_w = 0.2  (boundary case)
    let flange = 0.02_f64;

    let p = Profile2d::h_section(half_w, half_h, web, flange).expect("ctor ok");
    assert!(
        matches!(
            p.outline(),
            Err(KernelError::NonPositiveDimension {
                name: "flange_width_minus_web",
                ..
            })
        ),
        "expected NonPositiveDimension error for web == b, got {:?}",
        p.outline()
    );
}

/// A valid H-section with `web < b` still succeeds after the fix.
#[test]
fn h_section_still_accepts_valid_web() {
    let half_w = 0.1_f64;
    let half_h = 0.2_f64;
    let web = 0.01_f64; // << b = 0.2
    let flange = 0.02_f64;

    let p = Profile2d::h_section(half_w, half_h, web, flange).expect("ctor ok");
    match p.outline().expect("valid outline") {
        ProfileGeom::Polygon(v) => {
            assert_eq!(v.len(), 12usize);
            // Verify no zero-length edges (all consecutive vertices distinct).
            let n = v.len();
            for i in 0..n {
                let a = v[i];
                let b = v[(i + 1) % n];
                let dx = (a[0] - b[0]).abs();
                let dy = (a[1] - b[1]).abs();
                assert!(
                    dx > 1e-12_f64 || dy > 1e-12_f64,
                    "zero-length edge at index {i}: {a:?} == {b:?}"
                );
            }
        }
        ProfileGeom::Circle { .. } => panic!("H must be a polygon"),
        _ => panic!("unexpected ProfileGeom variant"),
    }
}

/// `web` strictly greater than `b` must still be rejected (regression guard).
#[test]
fn h_section_rejects_web_greater_than_b() {
    let half_w = 0.1_f64;
    let half_h = 0.2_f64;
    let web = 0.3_f64; // > b = 0.2
    let flange = 0.02_f64;

    let p = Profile2d::h_section(half_w, half_h, web, flange).expect("ctor ok");
    assert!(
        matches!(p.outline(), Err(KernelError::NonPositiveDimension { .. })),
        "expected error for web > b"
    );
}

// ── Bug 2: Ellipse3::new missing major_dir ⊥ normal check ────────────────────

/// Before the fix: `Ellipse3::new` with `major_dir == normal` returned `Ok`,
/// producing an invalid instance where `point_at` degenerates.
/// After the fix: this must return `Err(MajorDirNotInPlane)`.
#[test]
fn ellipse3_rejects_major_dir_parallel_to_normal() {
    // normal = major_dir = Z → completely parallel, should be rejected.
    let result = Ellipse3::new(Point3::origin(), Vec3::Z, Vec3::Z, 2.0_f64, 1.0_f64);
    assert!(
        matches!(result, Err(KernelError::MajorDirNotInPlane { .. })),
        "expected MajorDirNotInPlane error, got {:?}",
        result
    );
}

/// Anti-parallel case: major_dir = -normal.  Also nearly parallel → error.
#[test]
fn ellipse3_rejects_major_dir_antiparallel_to_normal() {
    let result = Ellipse3::new(
        Point3::origin(),
        Vec3::Z,
        Vec3::new(0.0_f64, 0.0_f64, -1.0_f64),
        2.0_f64,
        1.0_f64,
    );
    assert!(
        matches!(result, Err(KernelError::MajorDirNotInPlane { .. })),
        "expected MajorDirNotInPlane error, got {:?}",
        result
    );
}

/// A perfectly orthogonal `major_dir` is still accepted and `point_at` is correct.
#[test]
fn ellipse3_accepts_orthogonal_major_dir_and_point_at_is_correct() {
    // normal = Z, major_dir = X  (perfectly orthogonal)
    let semi_major = 2.0_f64;
    let semi_minor = 1.0_f64;
    let e = Ellipse3::new(Point3::origin(), Vec3::Z, Vec3::X, semi_major, semi_minor)
        .expect("orthogonal major_dir must be accepted");

    // t = 0 → point along major axis at distance semi_major.
    let p0 = e.point_at(0.0_f64);
    assert!(
        (p0.x - semi_major).abs() < 1e-12_f64,
        "p0.x = {}, expected {}",
        p0.x,
        semi_major
    );
    assert!(p0.y.abs() < 1e-12_f64, "p0.y = {}, expected 0", p0.y);
    assert!(p0.z.abs() < 1e-12_f64, "p0.z = {}, expected 0", p0.z);

    // t = π/2 → point along minor axis at distance semi_minor.
    let p1 = e.point_at(std::f64::consts::FRAC_PI_2);
    assert!(p1.x.abs() < 1e-12_f64, "p1.x = {}, expected 0", p1.x);
    assert!(
        (p1.y - semi_minor).abs() < 1e-12_f64,
        "p1.y = {}, expected {}",
        p1.y,
        semi_minor
    );
    assert!(p1.z.abs() < 1e-12_f64, "p1.z = {}, expected 0", p1.z);
}

/// A slightly non-orthogonal `major_dir` (within angular tolerance) must be
/// silently Gram-Schmidt corrected and the resulting `point_at` must still be
/// correct to within geometric tolerance.
#[test]
fn ellipse3_corrects_slightly_non_orthogonal_major_dir() {
    let tol = Tol::default();
    // Inject a tiny out-of-plane component well within tol.angular.
    let eps = tol.angular * 0.1_f64; // 1e-10 << 1e-9
    let major_dir = Vec3::new(1.0_f64, 0.0_f64, eps); // nearly X, with tiny Z

    let semi_major = 3.0_f64;
    let semi_minor = 1.5_f64;
    let e = Ellipse3::new(Point3::origin(), Vec3::Z, major_dir, semi_major, semi_minor)
        .expect("slightly off-orthogonal must be accepted and corrected");

    // After Gram-Schmidt the major_dir is projected to lie in the XY plane.
    // t = π/2 must produce a point with |z| < geometric tolerance.
    let p = e.point_at(std::f64::consts::FRAC_PI_2);
    assert!(
        p.z.abs() < 1e-10_f64,
        "point_at(π/2).z = {} should be ~0 after Gram-Schmidt correction",
        p.z
    );
    // The point should be close to (0, semi_minor, 0).
    assert!(
        (p.y - semi_minor).abs() < 1e-10_f64,
        "point_at(π/2).y = {}, expected ~{}",
        p.y,
        semi_minor
    );
}
