//! Interference-check tests (`DESIGN.md` §10 Phase 7, §6-6).
//!
//! Covers the hard-clash two-phase pipeline (coarse AABB → fine volume), the
//! honest `PotentialClash` degeneracy for pairs with no common prismatic
//! direction, the clip-pair exclusion, sleeve verification (admissible /
//! diameter / end-distance / edge-distance), and the `member_from_axis` adapter
//! entry. Every literal carries an `f64` annotation and an explicit tolerance
//! (`DESIGN.md` §12).

use archi_kernel::clash::{
    check_sleeve, clash_check, member_from_axis, ClashKind, ClashOptions, SleeveRule,
    SleeveViolationKind,
};
use archi_kernel::csg::{CsgNode, Member, Opening, OpeningId, Profile2d, StableId};
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::model::Model;
use archi_kernel::tolerance::Tol;

const VOL_EPS: f64 = 1e-9;

// ── geometry helpers ─────────────────────────────────────────────────────────

/// An axis-aligned box `[x0,x0+wx]×[y0,y0+wy]×[z0,z0+wz]` extruded along `+z`.
/// The `+z` profile axes are `(u, v) = (Y, −X)`, so the half-width is the `Y`
/// half and the half-height is the `X` half (matches `tests/prismatic.rs`).
fn box_z(x0: f64, y0: f64, z0: f64, wx: f64, wy: f64, wz: f64) -> CsgNode {
    CsgNode::Extrude {
        profile: Profile2d::rect(wy / 2.0, wx / 2.0).expect("valid rect"),
        origin: Point3::new(x0 + wx / 2.0, y0 + wy / 2.0, z0),
        axis: Vec3::Z,
        length: wz,
    }
}

/// An axis-aligned box extruded along `+x`. `(u, v) = (Z, −Y)`.
fn box_x(x0: f64, y0: f64, z0: f64, wx: f64, wy: f64, wz: f64) -> CsgNode {
    CsgNode::Extrude {
        profile: Profile2d::rect(wz / 2.0, wy / 2.0).expect("valid rect"),
        origin: Point3::new(x0, y0 + wy / 2.0, z0 + wz / 2.0),
        axis: Vec3::X,
        length: wx,
    }
}

fn id(n: u64) -> StableId {
    StableId(n)
}

// ── hard clash: rectangular beams crossing ───────────────────────────────────

#[test]
fn crossing_rectangular_beams_hard_clash_exact_volume() {
    let tol = Tol::default();
    let mut model = Model::new();

    // Beam A along x: section 0.3 (y) × 0.3 (z), x ∈ [0, 2].
    let beam_a = box_x(0.0_f64, -0.15_f64, -0.15_f64, 2.0_f64, 0.3_f64, 0.3_f64);
    // Beam B along z: section 0.2 (x) × 0.3 (y), passing through at x ∈ [0.9,1.1].
    let beam_b = box_z(0.9_f64, -0.15_f64, -0.5_f64, 0.2_f64, 0.3_f64, 1.0_f64);

    model.insert(id(1), Member::new(beam_a)).expect("insert a");
    model.insert(id(2), Member::new(beam_b)).expect("insert b");

    let result = clash_check(&mut model, &tol, &ClashOptions::from_tol(&tol));
    assert!(
        result.errors.is_empty(),
        "no eval errors: {:?}",
        result.errors
    );
    assert_eq!(result.clashes.len(), 1usize, "exactly one pair clashes");

    let c = result.clashes[0];
    assert_eq!((c.a, c.b), (id(1), id(2)));
    // Overlap box: x ∈ [0.9,1.1] (0.2), y ∈ [−0.15,0.15] (0.3), z ∈ [−0.15,0.15]
    // (0.3) = 0.2·0.3·0.3 = 0.018 m³.
    match c.kind {
        ClashKind::HardClash { volume } => {
            assert!(
                (volume - 0.018_f64).abs() <= VOL_EPS,
                "intersection volume = {volume}"
            );
        }
        other => panic!("expected HardClash, got {other:?}"),
    }
}

// ── hard clash: a pipe (cylinder) piercing a beam ─────────────────────────────

#[test]
fn pipe_piercing_beam_hard_clash() {
    let tol = Tol::default();
    let mut model = Model::new();

    // Beam along x: 0.4 (y) × 0.6 (z) section, x ∈ [0, 3].
    let beam = box_x(0.0_f64, -0.2_f64, -0.3_f64, 3.0_f64, 0.4_f64, 0.6_f64);
    // Round pipe of radius 0.1 along y through the beam at x = 1.5, z = 0,
    // y ∈ [−0.5, 0.5] (so it fully crosses the 0.4-deep beam in y).
    let pipe = CsgNode::Extrude {
        profile: Profile2d::circle(0.1_f64).expect("valid circle"),
        origin: Point3::new(1.5_f64, -0.5_f64, 0.0_f64),
        axis: Vec3::Y,
        length: 1.0_f64,
    };

    model.insert(id(1), Member::new(beam)).expect("insert beam");
    model.insert(id(2), Member::new(pipe)).expect("insert pipe");

    let result = clash_check(&mut model, &tol, &ClashOptions::from_tol(&tol));
    assert!(result.errors.is_empty(), "no eval errors");
    assert_eq!(result.clashes.len(), 1usize);

    // The pipe and beam share a common prismatic direction (the pipe axis y is a
    // box axis of the beam), so the exact volume is available: a cylinder of
    // radius 0.1 over the y-span where the beam is present (y ∈ [−0.2, 0.2],
    // length 0.4): π·0.1²·0.4 = 0.0125663706 m³.
    let c = result.clashes[0];
    match c.kind {
        ClashKind::HardClash { volume } => {
            let expected = std::f64::consts::PI * 0.1_f64 * 0.1_f64 * 0.4_f64;
            assert!(
                (volume - expected).abs() <= 1e-6_f64,
                "pipe∩beam volume = {volume}, expected {expected}"
            );
        }
        other => panic!("expected HardClash, got {other:?}"),
    }
}

// ── honest degeneracy: H × H crossed has no common direction ──────────────────

#[test]
fn crossed_h_sections_report_potential_clash() {
    let tol = Tol::default();
    let mut model = Model::new();

    // Two H-section members crossing at a right angle: an H is prismatic only
    // along its own axis, so the pair shares no common prismatic direction and
    // the exact volume is undecidable in this phase → PotentialClash.
    let h = || Profile2d::h_section(0.2_f64, 0.25_f64, 0.012_f64, 0.02_f64).expect("valid H");
    let beam_x_h = CsgNode::Extrude {
        profile: h(),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::X,
        length: 2.0_f64,
    };
    let beam_y_h = CsgNode::Extrude {
        profile: h(),
        origin: Point3::new(1.0_f64, -1.0_f64, 0.0_f64),
        axis: Vec3::Y,
        length: 2.0_f64,
    };

    model
        .insert(id(1), Member::new(beam_x_h))
        .expect("insert x");
    model
        .insert(id(2), Member::new(beam_y_h))
        .expect("insert y");

    let result = clash_check(&mut model, &tol, &ClashOptions::from_tol(&tol));
    assert!(result.errors.is_empty(), "no eval errors");
    assert_eq!(result.clashes.len(), 1usize);
    assert_eq!(
        result.clashes[0].kind,
        ClashKind::PotentialClash,
        "crossed H-sections must report PotentialClash, never a silent clear"
    );
}

// ── coarse phase: far-apart members produce no result ─────────────────────────

#[test]
fn far_apart_members_coarse_rejected() {
    let tol = Tol::default();
    let mut model = Model::new();

    // A grid of boxes well separated along x — every pair's AABBs are disjoint.
    for k in 0..8u64 {
        let x0 = 10.0_f64 * (k as f64);
        let b = box_z(x0, 0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64, 3.0_f64);
        model.insert(id(k), Member::new(b)).expect("insert");
    }

    let result = clash_check(&mut model, &tol, &ClashOptions::from_tol(&tol));
    assert!(result.errors.is_empty(), "no eval errors");
    assert!(
        result.clashes.is_empty(),
        "well-separated members must be coarse-rejected, got {:?}",
        result.clashes
    );
}

// ── clip pairs are intentional, not clashes ───────────────────────────────────

#[test]
fn clipped_girder_and_column_not_reported_by_default() {
    let tol = Tol::default();

    // A column (along z) and a girder (along x) crossing it. The girder is
    // priority-clipped by the column, so their overlap is intentional.
    let column = box_z(0.9_f64, -0.5_f64, 0.0_f64, 0.6_f64, 1.0_f64, 3.0_f64);
    let girder_gross = box_x(0.0_f64, -0.2_f64, 1.0_f64, 3.0_f64, 0.4_f64, 0.6_f64);
    let girder = CsgNode::Clip {
        base: Box::new(girder_gross),
        clippers: vec![id(1)],
        rule: archi_kernel::csg::ClipRule::Priority,
    };

    let mut model = Model::new();
    model
        .insert(id(1), Member::new(column))
        .expect("insert column");
    model
        .insert(id(2), Member::new(girder))
        .expect("insert girder");

    // Default: clip pairs excluded → no clash reported.
    let excluded = clash_check(&mut model, &tol, &ClashOptions::from_tol(&tol));
    assert!(
        excluded.errors.is_empty(),
        "no eval errors: {:?}",
        excluded.errors
    );
    assert!(
        excluded.clashes.is_empty(),
        "intentional clip overlap must be excluded by default, got {:?}",
        excluded.clashes
    );

    // Opt-in: include clip pairs → the overlap surfaces (the gross prisms do
    // overlap), as either a HardClash or PotentialClash, but it *is* reported.
    let mut opts = ClashOptions::from_tol(&tol);
    opts.include_clip_pairs = true;
    let included = clash_check(&mut model, &tol, &opts);
    assert_eq!(
        included.clashes.len(),
        1usize,
        "including clip pairs surfaces the overlap"
    );
}

// ── sleeve verification ───────────────────────────────────────────────────────

/// Build an RC beam (rectangular, depth = せい) along `+x` with one circular
/// sleeve. The beam profile is `Rect { half_w, half_h }`; for a `+x` extrusion
/// the depth direction `v = −Y` and the half-depth is `half_h`. The sleeve is a
/// circular extrusion whose centre line passes through `(cx, cy_offset_along_v,
/// _)`. We position the sleeve by its axial coordinate `cx` and its depth offset
/// `off` (signed along `v = −Y`, i.e. world `y = −off`).
fn rc_beam_with_sleeve(
    depth: f64,
    length: f64,
    radius: f64,
    cx: f64,
    depth_off: f64,
    opening_id: u64,
) -> CsgNode {
    let half_w = 0.15_f64; // beam width / 2 (along u = Z), irrelevant to the rules
    let base = CsgNode::Extrude {
        profile: Profile2d::rect(half_w, depth / 2.0).expect("valid rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::X,
        length,
    };
    // v = −Y for a +x extrusion, so a depth offset of `depth_off` along v sits at
    // world y = −depth_off. The sleeve runs along the width (Z) so it crosses the
    // web; its origin is on the web centre line at the chosen axial/depth point.
    let sleeve = CsgNode::Extrude {
        profile: Profile2d::circle(radius).expect("valid circle"),
        origin: Point3::new(cx, -depth_off, -half_w),
        axis: Vec3::Z,
        length: 2.0_f64 * half_w,
    };
    CsgNode::OpeningSubtraction {
        base: Box::new(base),
        openings: vec![(OpeningId(opening_id), Opening { shape: sleeve })],
    }
}

#[test]
fn sleeve_admissible_passes() {
    let tol = Tol::default();
    // Depth 0.9 m, span 6 m. Sleeve Ø 0.2 (≤ 0.9/3 = 0.3), centred axially at
    // x = 3 (end distance 3 ≥ 1.5·0.9 = 1.35), mid-depth (off 0 → edge distance
    // 0.45 ≥ 0.9/3 = 0.3). Admissible.
    let beam = rc_beam_with_sleeve(0.9_f64, 6.0_f64, 0.1_f64, 3.0_f64, 0.0_f64, 7);
    let report = check_sleeve(&beam, &SleeveRule::default(), &tol).expect("rectangular beam");
    assert_eq!(report.beam_depth, 0.9_f64);
    assert_eq!(report.sleeves_checked, 1usize);
    assert!(
        report.is_admissible(),
        "sleeve should be admissible: {:?}",
        report.violations
    );
}

#[test]
fn sleeve_diameter_too_large_flagged() {
    let tol = Tol::default();
    // Ø 0.4 (radius 0.2) > 0.9/3 = 0.3 → DiameterTooLarge. Positioned fine.
    let beam = rc_beam_with_sleeve(0.9_f64, 6.0_f64, 0.2_f64, 3.0_f64, 0.0_f64, 7);
    let report = check_sleeve(&beam, &SleeveRule::default(), &tol).expect("rectangular beam");
    let v = report
        .violations
        .iter()
        .find(|v| v.kind == SleeveViolationKind::DiameterTooLarge)
        .expect("diameter violation");
    assert_eq!(v.opening, OpeningId(7));
    assert!(
        (v.measured - 0.4_f64).abs() <= VOL_EPS,
        "measured Ø = {}",
        v.measured
    );
    assert!((v.limit - 0.3_f64).abs() <= VOL_EPS, "limit = {}", v.limit);
}

#[test]
fn sleeve_end_distance_too_small_flagged() {
    let tol = Tol::default();
    // Sleeve at x = 0.5; end distance 0.5 < 1.5·0.9 = 1.35 → EndDistanceTooSmall.
    let beam = rc_beam_with_sleeve(0.9_f64, 6.0_f64, 0.1_f64, 0.5_f64, 0.0_f64, 9);
    let report = check_sleeve(&beam, &SleeveRule::default(), &tol).expect("rectangular beam");
    let v = report
        .violations
        .iter()
        .find(|v| v.kind == SleeveViolationKind::EndDistanceTooSmall)
        .expect("end-distance violation");
    assert_eq!(v.opening, OpeningId(9));
    assert!(
        (v.measured - 0.5_f64).abs() <= VOL_EPS,
        "measured = {}",
        v.measured
    );
    assert!((v.limit - 1.35_f64).abs() <= VOL_EPS, "limit = {}", v.limit);
}

#[test]
fn sleeve_too_close_to_top_edge_flagged() {
    let tol = Tol::default();
    // Depth 0.9, half-depth 0.45. Offset the sleeve 0.35 toward an edge: edge
    // distance 0.45 − 0.35 = 0.10 < 0.9/3 = 0.30 → EdgeDistanceTooSmall.
    let beam = rc_beam_with_sleeve(0.9_f64, 6.0_f64, 0.1_f64, 3.0_f64, 0.35_f64, 11);
    let report = check_sleeve(&beam, &SleeveRule::default(), &tol).expect("rectangular beam");
    let v = report
        .violations
        .iter()
        .find(|v| v.kind == SleeveViolationKind::EdgeDistanceTooSmall)
        .expect("edge-distance violation");
    assert_eq!(v.opening, OpeningId(11));
    assert!(
        (v.measured - 0.10_f64).abs() <= 1e-9_f64,
        "measured = {}",
        v.measured
    );
    assert!(
        (v.limit - 0.30_f64).abs() <= 1e-9_f64,
        "limit = {}",
        v.limit
    );
}

// ── member_from_axis adapter entry ────────────────────────────────────────────

#[test]
fn member_from_axis_column_volume_matches() {
    let tol = Tol::default();
    // A 0.6 × 0.6 RC column, 4 m tall, built from its two axis nodes.
    let start = Point3::new(0.0_f64, 0.0_f64, 0.0_f64);
    let end = Point3::new(0.0_f64, 0.0_f64, 4.0_f64);
    let profile = Profile2d::rect(0.3_f64, 0.3_f64).expect("valid rect");
    let node = member_from_axis(profile, start, end).expect("distinct nodes");

    let mut member = Member::new(node);
    let vol = member.brep(&tol).expect("evaluate").signed_volume();
    // 0.6 · 0.6 · 4.0 = 1.44 m³.
    assert!((vol - 1.44_f64).abs() <= VOL_EPS, "column volume = {vol}");
}

#[test]
fn member_from_axis_beam_volume_matches() {
    let tol = Tol::default();
    // A 0.4 × 0.7 beam spanning 5 m horizontally from start to end.
    let start = Point3::new(1.0_f64, 2.0_f64, 3.0_f64);
    let end = Point3::new(6.0_f64, 2.0_f64, 3.0_f64);
    let profile = Profile2d::rect(0.2_f64, 0.35_f64).expect("valid rect");
    let node = member_from_axis(profile, start, end).expect("distinct nodes");

    let mut member = Member::new(node);
    let vol = member.brep(&tol).expect("evaluate").signed_volume();
    // 0.4 · 0.7 · 5.0 = 1.4 m³.
    assert!((vol - 1.4_f64).abs() <= VOL_EPS, "beam volume = {vol}");
}

#[test]
fn member_from_axis_rejects_zero_length() {
    let p = Point3::new(1.0_f64, 1.0_f64, 1.0_f64);
    let profile = Profile2d::rect(0.3_f64, 0.3_f64).expect("valid rect");
    assert!(member_from_axis(profile, p, p).is_err());
}
