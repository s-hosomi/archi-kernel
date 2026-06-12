//! Regression tests for two prismatic boolean bugs:
//!
//! 1. **void-void corner pinch** (`boolean/prismatic/build.rs`): two openings
//!    (voids) that touch at exactly one cross-section corner pinch the surrounding
//!    material into a non-manifold vertical edge (four walls sharing one curve),
//!    which `validate(Full)` rejected. This is the dual of the material-material
//!    checkerboard corner touch.
//!
//! 2. **vertex-on-edge (grazing) projection** (`boolean/prismatic/arrange.rs`): a
//!    sub-tolerance near-collinear overlapping pair of vertical edges emptied the
//!    result (covered by the now-active `prism_collinear_near_overlap_sliver_*`
//!    test in `adversarial_boolean.rs`; an extra ULP-window check is included here
//!    for defence in depth).

use archi_kernel::boolean::prismatic::{self, ExtrudeLeaf};
use archi_kernel::csg::Profile2d;
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::ValidateLevel;

/// A box with its cross-section's lower-left corner at `(x0, y0)`, width `wx`/`wy`,
/// extruded along `+z` from `z0` by height `wz`.
fn box_z(x0: f64, y0: f64, z0: f64, wx: f64, wy: f64, wz: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(wy / 2.0, wx / 2.0).expect("valid rect"),
        origin: Point3::new(x0 + wx / 2.0, y0 + wy / 2.0, z0),
        axis: Vec3::Z,
        length: wz,
    }
}

// ── Fix 1: void-void corner pinch ───────────────────────────────────────────

/// Two openings touching at the cross-section corner `(1,1)` must yield a valid
/// (manifold, `validate(Full)`-clean) result whose volume is the base minus the
/// two unit voids.
#[test]
fn void_void_corner_pinch_opening_subtraction() {
    let tol = Tol::default();
    // base x[0,3] y[0,3] height 1
    let base = box_z(0.0, 0.0, 0.0, 3.0, 3.0, 1.0);
    // void_A x[0,1] y[0,1]; void_B x[1,2] y[1,2]; they touch at corner (1,1).
    let void_a = box_z(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let void_b = box_z(1.0, 1.0, 0.0, 1.0, 1.0, 1.0);

    let brep = prismatic::opening_subtraction(&base, &[void_a, void_b], &tol)
        .expect("opening_subtraction must succeed (corner pinch)");
    brep.validate(&tol, ValidateLevel::Full)
        .expect("result must be manifold / watertight");
    let v: f64 = brep.signed_volume();
    // 9 (base) − 1 − 1 = 7.
    assert!((v - 7.0_f64).abs() <= 1e-9_f64, "volume = {v}, expected 7");
}

/// Control: shifting one void by `2e-6` (> `Tol::length`) so the corners do *not*
/// touch yields the same valid single solid — isolating the corner coincidence as
/// the trigger.
#[test]
fn void_void_corner_offset_passes() {
    let tol = Tol::default();
    let base = box_z(0.0, 0.0, 0.0, 3.0, 3.0, 1.0);
    let void_a = box_z(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let void_b = box_z(1.0 + 2e-6, 1.0 + 2e-6, 0.0, 1.0, 1.0, 1.0);

    let brep = prismatic::opening_subtraction(&base, &[void_a, void_b], &tol)
        .expect("offset case must succeed");
    brep.validate(&tol, ValidateLevel::Full)
        .expect("offset result must be valid");
    assert_eq!(brep.solids.len(), 1);
    assert!((brep.signed_volume() - 7.0_f64).abs() <= 1e-9_f64);
}

/// A slab pierced by two diagonal square columns touching at the slab centre
/// corner `(0,0)` — the everyday architectural form of the pinch.
#[test]
fn void_void_corner_pinch_diag_two_columns() {
    let tol = Tol::default();
    // slab x[-1,1] y[-1,1] thickness 0.2
    let slab = box_z(-1.0, -1.0, 0.0, 2.0, 2.0, 0.2);
    let col = |cx: f64, cy: f64| box_z(cx - 0.25, cy - 0.25, 0.0, 0.5, 0.5, 0.2);
    // Two diagonal columns touching at corner (0,0).
    let cols = [col(-0.25, -0.25), col(0.25, 0.25)];

    let brep = prismatic::opening_subtraction(&slab, &cols, &tol)
        .expect("diag two columns must succeed (corner pinch)");
    brep.validate(&tol, ValidateLevel::Full)
        .expect("result must be manifold / watertight");
    // slab vol 0.8 minus two 0.5×0.5×0.2 = 0.05 columns ⇒ 0.7.
    let v: f64 = brep.signed_volume();
    assert!(
        (v - 0.7_f64).abs() <= 1e-9_f64,
        "volume = {v}, expected 0.7"
    );
}

/// Control: four columns at all four quadrant corners (the slab centre becomes a
/// fully-removed cross point, not a pinch) must also succeed.
#[test]
fn void_four_columns_corner_passes() {
    let tol = Tol::default();
    let slab = box_z(-1.0, -1.0, 0.0, 2.0, 2.0, 0.2);
    let col = |cx: f64, cy: f64| box_z(cx - 0.25, cy - 0.25, 0.0, 0.5, 0.5, 0.2);
    let cols = [
        col(-0.25, -0.25),
        col(0.25, 0.25),
        col(-0.25, 0.25),
        col(0.25, -0.25),
    ];
    let brep =
        prismatic::opening_subtraction(&slab, &cols, &tol).expect("four columns must succeed");
    brep.validate(&tol, ValidateLevel::Full)
        .expect("result must be valid");
    let v: f64 = brep.signed_volume();
    // slab 0.8 minus four 0.05 columns = 0.6.
    assert!(
        (v - 0.6_f64).abs() <= 1e-9_f64,
        "volume = {v}, expected 0.6"
    );
}

// ── Fix 2: vertex-on-edge grazing projection (ULP window) ────────────────────

/// The conservation identity `V(A∩B) + V(A−B) == V(A)` must hold across the small
/// ULP window where B's right edge grazes A's right edge `x = 1.0` from just
/// inside — the configuration that emptied the result before the grazing
/// projection landed.
#[test]
fn prism_grazing_ulp_window_conserves_volume() {
    let tol = Tol::default();
    let a = box_z(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let x0 = -1.55_f64;
    for k in -2_i64..=2 {
        let right = 1.0_f64 + (k as f64) * f64::EPSILON;
        let wx = right - x0;
        let b = box_z(x0, 0.1, 0.0, wx, 2.9, 0.2);
        let vd = prismatic::difference(&a, &b, &tol)
            .expect("difference")
            .signed_volume();
        let vi = prismatic::intersection(&a, &b, &tol)
            .expect("intersection")
            .signed_volume();
        let sum = vd + vi;
        assert!(
            (sum - 1.0_f64).abs() <= 1e-9_f64,
            "k={k} right={right:.17}: V(A−B)+V(A∩B) = {sum} != 1.0",
        );
    }
}
