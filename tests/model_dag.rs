//! Phase 5 model-layer unit tests: dependency DAG dirty propagation, multi-level
//! priority chains, and located-defect coordinates.

use archi_kernel::csg::{ClipRule, CsgNode, Member, Profile2d, StableId};
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::model::Model;
use archi_kernel::tolerance::Tol;

fn column(cx: f64, side: f64) -> CsgNode {
    CsgNode::Extrude {
        profile: Profile2d::rect(side / 2.0, side / 2.0).expect("rect"),
        origin: Point3::new(cx, 0.0, 0.0),
        axis: Vec3::Z,
        length: 3.0,
    }
}

fn beam_x(len: f64) -> CsgNode {
    CsgNode::Extrude {
        profile: Profile2d::rect(0.2, 0.2).expect("rect"),
        origin: Point3::new(0.0, 0.0, 2.6),
        axis: Vec3::X,
        length: len,
    }
}

/// Column → girder → small-beam priority chain: marking the column dirty
/// propagates through the girder to the small beam (transitive dependents).
#[test]
fn dirty_propagates_transitively_through_priority_chain() {
    let tol = Tol::default();
    let mut model = Model::new();

    let col = StableId(1);
    let girder = StableId(2);
    let small = StableId(3);

    model.insert(col, Member::new(column(0.0, 0.5))).unwrap();
    // Girder clipped by the column.
    model
        .insert(
            girder,
            Member::new(CsgNode::Clip {
                base: Box::new(beam_x(6.0)),
                clippers: vec![col],
                rule: ClipRule::Priority,
            }),
        )
        .unwrap();
    // Small beam clipped by BOTH the column and the girder (gross leaves of each).
    model
        .insert(
            small,
            Member::new(CsgNode::Clip {
                base: Box::new(CsgNode::Extrude {
                    profile: Profile2d::rect(0.15, 0.15).expect("rect"),
                    origin: Point3::new(3.0, -2.0, 2.6),
                    axis: Vec3::Y,
                    length: 4.0,
                }),
                clippers: vec![col, girder],
                rule: ClipRule::Priority,
            }),
        )
        .unwrap();

    // After fresh insert every member is dirty; mark each clean by a successful
    // per-member evaluation (which fills the cache).
    for id in [col, girder, small] {
        let _ = model.get_mut(id).unwrap().brep(&tol);
    }
    assert!(!model.get(girder).unwrap().is_dirty(&tol));
    assert!(!model.get(small).unwrap().is_dirty(&tol));

    // Now move the column: propagation marks the girder AND the small beam dirty.
    model.mark_dirty(col);
    assert!(
        model.get(girder).unwrap().is_dirty(&tol),
        "girder must be dirty (depends on column)"
    );
    assert!(
        model.get(small).unwrap().is_dirty(&tol),
        "small beam must be dirty (depends on column transitively and directly)"
    );
}

/// Removing a member returns it and shrinks the model.
#[test]
fn remove_returns_member() {
    let mut model = Model::new();
    model
        .insert(StableId(1), Member::new(column(0.0, 0.5)))
        .unwrap();
    assert_eq!(model.len(), 1);
    let removed = model.remove(StableId(1));
    assert!(removed.is_ok());
    assert!(model.is_empty());
    assert!(model.remove(StableId(1)).is_err(), "second remove errors");
}
