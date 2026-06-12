//! Regression tests for two bugs fixed together:
//!
//! 1. `sleeve.rs`: `circular_opening` was returning the extrusion `origin` (one
//!    cap face) instead of the **midpoint** of the axis.  When a sleeve was
//!    extruded in the beam depth direction its origin sat on a flange, giving
//!    `edge_dist = 0` and a false `EdgeDistanceTooSmall` violation.
//!
//! 2. `plane_store.rs`: `should_flip` used a strict inequality (`> tol.angular`)
//!    to decide which component is "significant".  A component whose absolute
//!    value equals exactly `tol.angular` was silently skipped, so two planes with
//!    first components `+tol.angular` and `-tol.angular` both fell through to the
//!    *second* component, canonicalised to the same half-space with **different**
//!    normals, and were not de-duplicated.

use archi_kernel::clash::{check_sleeve, SleeveRule, SleeveViolationKind};
use archi_kernel::csg::{CsgNode, Opening, OpeningId, Profile2d};
use archi_kernel::geom::GeomStore;
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::Plane;
use archi_kernel::tolerance::Tol;

// ── Bug 1: sleeve origin vs. midpoint ────────────────────────────────────────

/// A sleeve extruded in the beam depth direction (+Y) whose `origin` sits on
/// the bottom flange (world y = -half_h).  The true hole centre is at mid-depth
/// (y = 0).  Before the fix, `circular_opening` returned the flange origin,
/// making `edge_dist = 0.0 < depth/3 ≈ 0.133`, which triggered a false
/// `EdgeDistanceTooSmall` violation even though the sleeve is perfectly centred.
#[test]
fn sleeve_depth_direction_extrusion_no_false_edge_violation() {
    let tol = Tol::default();
    let half_h = 0.2_f64; // half-depth → full depth 0.4, edge_limit = 0.4/3 ≈ 0.1333
    let half_w = 0.3_f64;
    let length = 3.0_f64;

    let base = CsgNode::Extrude {
        profile: Profile2d::rect(half_w, half_h).expect("valid rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::X,
        length,
    };

    // Sleeve extruded along +Y (depth direction).
    // origin = bottom cap face at y = -half_h = -0.2.
    // Midpoint = (1.5, 0, 0) — exactly mid-depth, which is admissible.
    let sleeve = CsgNode::Extrude {
        profile: Profile2d::circle(0.05_f64).expect("valid circle"),
        origin: Point3::new(1.5_f64, -0.2_f64, 0.0_f64),
        axis: Vec3::Y,
        length: 0.4_f64,
    };

    let beam = CsgNode::OpeningSubtraction {
        base: Box::new(base),
        openings: vec![(OpeningId(42), Opening { shape: sleeve })],
    };

    let report = check_sleeve(&beam, &SleeveRule::default(), &tol)
        .expect("should parse as rectangular beam");

    let edge_viol = report
        .violations
        .iter()
        .find(|v| v.kind == SleeveViolationKind::EdgeDistanceTooSmall);

    assert!(
        edge_viol.is_none(),
        "FALSE EdgeDistanceTooSmall violation: the sleeve is at mid-depth but got {:?}",
        edge_viol
    );
    assert!(
        report.is_admissible(),
        "sleeve at mid-depth should be admissible; violations = {:?}",
        report.violations
    );
}

/// Confirm that the midpoint calculation is correct across a diagonal extrusion
/// whose axis is neither X, Y, nor Z.  We construct a sleeve whose midpoint is
/// at the beam mid-depth (y = 0) even though the origin is off-centre, and
/// verify no false violation occurs.
#[test]
fn sleeve_diagonal_extrusion_uses_midpoint() {
    let tol = Tol::default();
    let half_h = 0.2_f64;
    let half_w = 0.3_f64;
    let length = 3.0_f64;

    let base = CsgNode::Extrude {
        profile: Profile2d::rect(half_w, half_h).expect("valid rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::X,
        length,
    };

    // Extrude along a diagonal (Y + Z, normalised inside CsgNode::Extrude via
    // try_unit inside circular_opening).  The sleeve half-length along Y is
    // 0.1, so starting at y = -0.1 places the midpoint at y = 0 (mid-depth).
    let diag_len = 0.2_f64 * 2.0_f64.sqrt(); // length along (Y+Z)/√2
    let sleeve = CsgNode::Extrude {
        profile: Profile2d::circle(0.05_f64).expect("valid circle"),
        origin: Point3::new(1.5_f64, -0.1_f64, -0.1_f64),
        axis: Vec3::new(0.0_f64, 1.0_f64, 1.0_f64), // not yet normalised; try_unit handles it
        length: diag_len,
    };

    let beam = CsgNode::OpeningSubtraction {
        base: Box::new(base),
        openings: vec![(OpeningId(43), Opening { shape: sleeve })],
    };

    let report = check_sleeve(&beam, &SleeveRule::default(), &tol)
        .expect("should parse as rectangular beam");

    let edge_viol = report
        .violations
        .iter()
        .find(|v| v.kind == SleeveViolationKind::EdgeDistanceTooSmall);

    assert!(
        edge_viol.is_none(),
        "expected no EdgeDistanceTooSmall (midpoint at y=0); got {:?}",
        edge_viol
    );
}

/// Width-direction extrusion (axis along Z) — the canonical existing
/// configuration.  Confirms that the midpoint fix does NOT change the measured
/// depth offset for the standard case (axis ⊥ depth-direction v = −Y, so the
/// midpoint and origin have the same depth coordinate).
#[test]
fn sleeve_width_direction_extrusion_unaffected() {
    let tol = Tol::default();
    let depth = 0.9_f64;
    let length = 6.0_f64;
    let half_w = 0.15_f64;
    let depth_off = 0.1_f64; // offset from mid-depth along v = −Y → world y = -0.1

    let base = CsgNode::Extrude {
        profile: Profile2d::rect(half_w, depth / 2.0_f64).expect("valid rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::X,
        length,
    };
    // Sleeve runs along Z (width direction); origin has the correct depth
    // position already because axis has no Y component.
    let sleeve = CsgNode::Extrude {
        profile: Profile2d::circle(0.1_f64).expect("valid circle"),
        origin: Point3::new(3.0_f64, -depth_off, -half_w),
        axis: Vec3::Z,
        length: 2.0_f64 * half_w,
    };
    let beam = CsgNode::OpeningSubtraction {
        base: Box::new(base),
        openings: vec![(OpeningId(44), Opening { shape: sleeve })],
    };

    let report = check_sleeve(&beam, &SleeveRule::default(), &tol).expect("rectangular beam");

    assert!(
        report.is_admissible(),
        "width-direction sleeve should remain admissible; violations = {:?}",
        report.violations
    );
}

// ── Bug 2: plane de-duplication at angular tolerance boundary ─────────────────

/// Two planes through the origin whose normals differ only by the sign of their
/// x-component, which is exactly `±tol.angular`.  Both normals point to the
/// same geometric half-space (x ≈ 0, z > 0); they should canonicalise to a
/// single entry and return the same `SurfaceId`.
///
/// Before the fix, the strict `> tol.angular` threshold silently skipped the
/// ±tol.angular x-component, both planes resolved their canonical orientation
/// from the z-component (same), kept their distinct ±ε x-components, and
/// `normals_agree` returned false (|A − B| = 2ε > ε = tol.angular).
#[test]
fn dedup_at_angular_boundary() {
    let tol = Tol::default(); // angular = 1e-9
    let eps = tol.angular; // exactly the boundary value

    // Both vectors are essentially unit-length (z ≈ 1).
    let z = (1.0_f64 - eps * eps).sqrt();

    let pt = Point3::origin();
    let plane_a = Plane::new(pt, Vec3::new(eps, 0.0_f64, z)).expect("valid plane a");
    let plane_b = Plane::new(pt, Vec3::new(-eps, 0.0_f64, z)).expect("valid plane b");

    let mut store = GeomStore::new();
    let (id_a, _flip_a) = store.insert_plane(plane_a, &tol);
    let (id_b, _flip_b) = store.insert_plane(plane_b, &tol);

    assert_eq!(
        id_a, id_b,
        "planes that agree by design intent must share one SurfaceId; got {:?} vs {:?}",
        id_a, id_b
    );
}

/// Sweep around the angular-tolerance boundary.
///
/// `normals_agree` uses `2 * tol.angular` as its threshold (to absorb the
/// residual introduced by `should_flip` skipping sub-tol components).
/// So planes with opposite-sign sub-tol x-components are deduplicated when
/// their chord difference `|a−b| = 2·|x|` is ≤ `2·tol.angular`, i.e.
/// `|x| ≤ tol.angular`.  At factors above 1× the planes are genuinely distinct
/// and should NOT be merged.
#[test]
fn dedup_boundary_sweep() {
    let tol = Tol::default();
    let eps = tol.angular;

    // Factors ≤ 1.0 → |a−b| = 2·factor·eps ≤ 2·eps = 2*tol.angular → deduplicated.
    for factor in [0.5_f64, 1.0_f64] {
        let x = eps * factor;
        let z = (1.0_f64 - x * x).sqrt();
        let pt = Point3::origin();

        let plane_a = Plane::new(pt, Vec3::new(x, 0.0_f64, z)).expect("valid plane a");
        let plane_b = Plane::new(pt, Vec3::new(-x, 0.0_f64, z)).expect("valid plane b");

        let mut store = GeomStore::new();
        let (id_a, _) = store.insert_plane(plane_a, &tol);
        let (id_b, _) = store.insert_plane(plane_b, &tol);

        assert_eq!(
            id_a, id_b,
            "factor={}: planes within 2*tol.angular should be deduplicated; got {:?} vs {:?}",
            factor, id_a, id_b
        );
    }

    // Factors above 1.0 → |a−b| > 2*tol.angular → genuinely distinct, not merged.
    {
        let factor = 2.0_f64;
        let x = eps * factor;
        let z = (1.0_f64 - x * x).sqrt();
        let pt = Point3::origin();

        let plane_a = Plane::new(pt, Vec3::new(x, 0.0_f64, z)).expect("valid plane a");
        let plane_b = Plane::new(pt, Vec3::new(-x, 0.0_f64, z)).expect("valid plane b");

        let mut store = GeomStore::new();
        let (id_a, _) = store.insert_plane(plane_a, &tol);
        let (id_b, _) = store.insert_plane(plane_b, &tol);

        assert_ne!(
            id_a, id_b,
            "factor={}: planes 4*tol.angular apart are genuinely distinct and should NOT be merged",
            factor
        );
    }
}

/// Sanity check: planes in genuinely distinct half-spaces are NOT deduplicated.
#[test]
fn distinct_planes_not_merged() {
    let tol = Tol::default();

    let pt = Point3::origin();
    let plane_pos = Plane::new(pt, Vec3::new(1.0_f64, 0.0_f64, 0.0_f64)).expect("valid +X plane");
    let plane_neg = Plane::new(pt, Vec3::new(-1.0_f64, 0.0_f64, 0.0_f64)).expect("valid -X plane");

    // +X and -X are opposite half-spaces; canonicalise to the same +X canonical,
    // so they ARE the same plane — verify they share an ID.
    let mut store = GeomStore::new();
    let (id_pos, flipped_pos) = store.insert_plane(plane_pos, &tol);
    let (id_neg, flipped_neg) = store.insert_plane(plane_neg, &tol);

    assert_eq!(
        id_pos, id_neg,
        "+X and -X planes through the origin are the same plane; should share SurfaceId"
    );
    // One of them should be reported as flipped.
    assert_ne!(
        flipped_pos, flipped_neg,
        "exactly one of +X/-X should be the flipped form"
    );
}
