//! Phase 5 acceptance: column-priority quantity take-off (`DESIGN.md` §6-2, §7).
//!
//! This is the commercial demo's basis: two RC columns and one RC girder, with
//! the girder clipped to its columns' inner faces. The take-off must reproduce
//! the hand-calculated concrete volume (inner-clear length) and formwork split,
//! the model's dirty propagation must make the girder follow a moved column, and
//! a clip cycle must isolate exactly the cyclic members.
//!
//! Every literal carries an `f64` annotation and an explicit tolerance
//! (`DESIGN.md` §12). The volume tolerance is `1e-9 m³` (axis-aligned boxes
//! accumulate only f64 round-off, far below this; the same bound the Phase 2/3
//! tests use).

use archi_kernel::csg::{
    ClipRule, CsgNode, EvalError, Member, Opening, OpeningId, Profile2d, StableId,
};
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::model::{takeoff, Model};
use archi_kernel::tolerance::Tol;

const EPS: f64 = 1e-9;

// ── Geometry helpers ─────────────────────────────────────────────────────────

/// A 0.5 m × 0.5 m RC column centred at `(cx, 0)` in plan, `z ∈ [0, 3]`,
/// extruded along `+Z`. The `+Z` profile axes are `(u, v) = (Y, −X)`, so the
/// profile half-width is the `Y` half (0.25) and the half-height the `X` half
/// (0.25).
fn column_node(cx: f64) -> CsgNode {
    CsgNode::Extrude {
        profile: Profile2d::rect(0.25_f64, 0.25_f64).expect("column rect"),
        origin: Point3::new(cx, 0.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 3.0_f64,
    }
}

/// The gross (un-clipped) RC girder extrusion: 0.4 m wide (`Y`) × 0.6 m deep
/// (`Z`), spanning `x ∈ [0, 6]` (column centre to column centre), sitting at the
/// column top (`z ∈ [2.4, 3.0]`), extruded along `+X`. The `+X` profile axes are
/// `(u, v) = (Z, −Y)`, so the profile half-width is the `Z` half (0.3) and the
/// half-height the `Y` half (0.2).
fn girder_extrude() -> CsgNode {
    CsgNode::Extrude {
        profile: Profile2d::rect(0.3_f64, 0.2_f64).expect("girder rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 2.7_f64),
        axis: Vec3::X,
        length: 6.0_f64,
    }
}

// Stable ids assigned by the caller.
const COL_LEFT: StableId = StableId(1);
const COL_RIGHT: StableId = StableId(2);
const GIRDER: StableId = StableId(3);

/// Build the standard model: two columns and a girder clipped by both columns
/// (column-priority: the girder lists the columns as clippers).
fn standard_model() -> Model {
    let mut model = Model::new();
    model
        .insert(COL_LEFT, Member::new(column_node(0.0_f64)))
        .expect("insert left column");
    model
        .insert(COL_RIGHT, Member::new(column_node(6.0_f64)))
        .expect("insert right column");
    let girder = CsgNode::Clip {
        base: Box::new(girder_extrude()),
        clippers: vec![COL_LEFT, COL_RIGHT],
        rule: ClipRule::Priority,
    };
    model
        .insert(GIRDER, Member::new(girder))
        .expect("insert girder");
    model
}

// ── Column-priority volume + dirty propagation ───────────────────────────────

#[test]
fn columns_full_volume_girder_inner_clear() {
    let tol = Tol::default();
    let mut model = standard_model();

    // Columns: full section × full height, no deduction (column priority means
    // nothing clips them). 0.5·0.5·3.0 = 0.75 each.
    let qty_left = takeoff(&mut model, COL_LEFT, &tol).expect("left column take-off");
    let qty_right = takeoff(&mut model, COL_RIGHT, &tol).expect("right column take-off");
    let col_vol = 0.5_f64 * 0.5_f64 * 3.0_f64;
    assert!(
        (qty_left.concrete_volume - col_vol).abs() <= EPS,
        "left column V = {}",
        qty_left.concrete_volume
    );
    assert!(
        (qty_right.concrete_volume - col_vol).abs() <= EPS,
        "right column V = {}",
        qty_right.concrete_volume
    );

    // Girder: inner-clear length = 6.0 − 0.5 = 5.5 (column inner faces at x =
    // 0.25 and x = 5.75); section 0.4 × 0.6. V = 5.5·0.4·0.6 = 1.32.
    let qty_girder = takeoff(&mut model, GIRDER, &tol).expect("girder take-off");
    let inner_len = 6.0_f64 - 0.5_f64;
    let girder_vol = inner_len * 0.4_f64 * 0.6_f64;
    assert!(
        (qty_girder.concrete_volume - girder_vol).abs() <= EPS,
        "girder V = {}, expected {girder_vol}",
        qty_girder.concrete_volume
    );
}

#[test]
fn girder_formwork_matches_hand_calc() {
    let tol = Tol::default();
    let mut model = standard_model();
    let qty = takeoff(&mut model, GIRDER, &tol).expect("girder take-off");

    let inner_len = 6.0_f64 - 0.5_f64; // 5.5
                                       // Side formwork: the two vertical web faces (normal ±Y), each inner_len ×
                                       // depth (0.6). The two column-contact end faces (normal ±X) are excluded.
    let side_expected = 2.0_f64 * inner_len * 0.6_f64; // 6.6
                                                       // Bottom formwork: the soffit (normal −Z), inner_len × width (0.4).
    let bottom_expected = inner_len * 0.4_f64; // 2.2
    assert!(
        (qty.formwork_side - side_expected).abs() <= EPS,
        "girder side formwork = {}, expected {side_expected}",
        qty.formwork_side
    );
    assert!(
        (qty.formwork_bottom - bottom_expected).abs() <= EPS,
        "girder bottom formwork = {}, expected {bottom_expected}",
        qty.formwork_bottom
    );
}

#[test]
fn moving_a_column_propagates_dirty_and_changes_girder_volume() {
    let tol = Tol::default();
    let mut model = standard_model();

    let before = takeoff(&mut model, GIRDER, &tol)
        .expect("girder before")
        .concrete_volume;
    let inner_len = 6.0_f64 - 0.5_f64;
    assert!((before - inner_len * 0.4_f64 * 0.6_f64).abs() <= EPS);

    // Widen the left column to a 0.7 m square (still centred at x = 0), so its
    // inner face moves out from x = 0.25 to x = 0.35. The girder follows: the
    // boolean re-clips against the larger column, deducting an extra
    // 0.1 m × 0.4 × 0.6 from the embedded end.
    {
        let col = model.get_mut(COL_LEFT).expect("left column present");
        *col.csg_mut() = CsgNode::Extrude {
            profile: Profile2d::rect(0.35_f64, 0.35_f64).expect("wide column rect"),
            origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
            axis: Vec3::Z,
            length: 3.0_f64,
        };
    }
    // The girder depends on the column: marking the column dirty must mark the
    // girder dirty too (column moves ⇒ beam follows).
    model.mark_dirty(COL_LEFT);
    assert!(
        model.get(GIRDER).expect("girder present").is_dirty(&tol),
        "girder must be dirty after the column it clips changed"
    );

    let after = takeoff(&mut model, GIRDER, &tol)
        .expect("girder after")
        .concrete_volume;
    // New inner-clear length: left inner face at x = 0.35, right at x = 5.75 ⇒
    // 5.75 − 0.35 = 5.4. V = 5.4·0.4·0.6 = 1.296.
    let new_inner = 5.75_f64 - 0.35_f64; // 5.4
    let expected = new_inner * 0.4_f64 * 0.6_f64;
    assert!(
        (after - expected).abs() <= EPS,
        "girder V after column widened = {after}, expected {expected}"
    );
    assert!(
        after < before,
        "the wider column deducts more, so the girder volume dropped"
    );
}

// ── Slab with an opening: formwork excludes the opening reveal ────────────────

#[test]
fn slab_with_opening_formwork_nets_the_hole() {
    let tol = Tol::default();
    // A 4 m × 4 m × 0.2 m slab (thin plate), extruded along +Z, base at z = 0.
    // +Z profile axes (u, v) = (Y, −X): half_w = Y half = 2.0, half_h = X half =
    // 2.0. A 1 m × 1 m vertical opening pierces it through the thickness.
    let slab_base = CsgNode::Extrude {
        profile: Profile2d::rect(2.0_f64, 2.0_f64).expect("slab rect"),
        origin: Point3::new(2.0_f64, 2.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 0.2_f64,
    };
    // Opening: 1 m × 1 m centred at (2, 2) in plan, through the 0.2 thickness
    // (z ∈ [−0.1, 0.3] covers it), extruded along +Z.
    let opening = CsgNode::Extrude {
        profile: Profile2d::rect(0.5_f64, 0.5_f64).expect("opening rect"),
        origin: Point3::new(2.0_f64, 2.0_f64, -0.1_f64),
        axis: Vec3::Z,
        length: 0.4_f64,
    };
    let slab = CsgNode::OpeningSubtraction {
        base: Box::new(slab_base),
        openings: vec![(OpeningId(1), Opening { shape: opening })],
    };

    let mut model = Model::new();
    model.insert(StableId(10), Member::new(slab)).expect("slab");
    let qty = takeoff(&mut model, StableId(10), &tol).expect("slab take-off");

    // Concrete volume: (4·4 − 1·1)·0.2 = 15·0.2 = 3.0.
    let vol = (4.0_f64 * 4.0_f64 - 1.0_f64 * 1.0_f64) * 0.2_f64;
    assert!(
        (qty.concrete_volume - vol).abs() <= EPS,
        "slab V = {}, expected {vol}",
        qty.concrete_volume
    );

    // Bottom formwork (the slab underside, normal −Z) is the net plan area, with
    // the opening removed because it appears as a hole loop in the bottom face:
    // 4·4 − 1·1 = 15.
    let bottom_expected = 4.0_f64 * 4.0_f64 - 1.0_f64 * 1.0_f64;
    assert!(
        (qty.formwork_bottom - bottom_expected).abs() <= EPS,
        "slab bottom formwork = {}, expected {bottom_expected}",
        qty.formwork_bottom
    );

    // Side formwork: the 4 outer edges (perimeter 16 × 0.2 = 3.2) plus the 4
    // opening-reveal faces (perimeter 4 × 0.2 = 0.8) = 4.0. The reveals are
    // vertical inner-loop walls, counted as side formwork.
    let side_expected = 16.0_f64 * 0.2_f64 + 4.0_f64 * 0.2_f64;
    assert!(
        (qty.formwork_side - side_expected).abs() <= EPS,
        "slab side formwork = {}, expected {side_expected}",
        qty.formwork_side
    );
}

// ── Cyclic clip: isolate the cycle, evaluate the rest ────────────────────────

#[test]
fn cyclic_clip_isolates_both_members_third_is_fine() {
    let tol = Tol::default();
    let mut model = Model::new();

    // A and B each clip the other — a dependency cycle. C is independent.
    let a = CsgNode::Clip {
        base: Box::new(column_node(0.0_f64)),
        clippers: vec![StableId(101)],
        rule: ClipRule::Priority,
    };
    let b = CsgNode::Clip {
        base: Box::new(column_node(3.0_f64)),
        clippers: vec![StableId(100)],
        rule: ClipRule::Priority,
    };
    let c = column_node(6.0_f64);
    model.insert(StableId(100), Member::new(a)).expect("A");
    model.insert(StableId(101), Member::new(b)).expect("B");
    model.insert(StableId(102), Member::new(c)).expect("C");

    let results = model.evaluate_all(&tol);

    // A and B are in the cycle: both isolated with CyclicDependency carrying both
    // ids; C evaluates normally.
    match results.get(&StableId(100)) {
        Some(Err(EvalError::CyclicDependency { members })) => {
            assert_eq!(members, &vec![100u64, 101u64], "cycle membership for A");
        }
        other => panic!("A should be CyclicDependency, got {other:?}"),
    }
    match results.get(&StableId(101)) {
        Some(Err(EvalError::CyclicDependency { members })) => {
            assert_eq!(members, &vec![100u64, 101u64], "cycle membership for B");
        }
        other => panic!("B should be CyclicDependency, got {other:?}"),
    }
    match results.get(&StableId(102)) {
        Some(Ok(brep)) => {
            let v = brep.signed_volume();
            assert!(
                (v - 0.5_f64 * 0.5_f64 * 3.0_f64).abs() <= EPS,
                "C (independent column) V = {v}"
            );
        }
        other => panic!("C should evaluate, got {other:?}"),
    }
}

// ── Duplicate id rejected ────────────────────────────────────────────────────

#[test]
fn duplicate_member_id_is_rejected() {
    let mut model = Model::new();
    model
        .insert(StableId(1), Member::new(column_node(0.0_f64)))
        .expect("first insert");
    assert!(
        model
            .insert(StableId(1), Member::new(column_node(1.0_f64)))
            .is_err(),
        "inserting a duplicate id must error"
    );
}

// ── Unknown clipper surfaces a machine-readable error ────────────────────────

#[test]
fn unknown_clipper_is_reported() {
    let tol = Tol::default();
    let mut model = Model::new();
    let girder = CsgNode::Clip {
        base: Box::new(girder_extrude()),
        clippers: vec![StableId(999)], // not in the model
        rule: ClipRule::Priority,
    };
    model.insert(GIRDER, Member::new(girder)).expect("girder");
    match takeoff(&mut model, GIRDER, &tol) {
        Err(EvalError::UnknownClipper { clipper }) => assert_eq!(clipper, 999u64),
        other => panic!("expected UnknownClipper, got {other:?}"),
    }
}
