//! Regression tests for topo validation bugs fixed in this batch:
//!
//! 1. `pair_siblings` MultipleSiblings: candidates not marked matched, causing
//!    spurious sibling-map entries and an inflated edge count.
//! 2. `validate_topology` empty-shell / empty-solid structural degeneracy
//!    silently returning `Ok(())`.

use archi_kernel::geom::{CurveGeom, GeomStore, SurfaceGeom, VertexGeom};
use archi_kernel::topo::validate::{validate_topology, Defect};
use archi_kernel::topo::{Face, HalfEdge, Loop, Sense, Shell, Solid, TopoStore, Vertex};
use archi_kernel::{Line3, Plane, Point3, Tol, Vec3};

// ── Bug 1: pair_siblings MultipleSiblings candidate marking ──────────────────

/// A 4-way clash `[0,1],[1,0],[0,1],[1,0]` on one curve must produce
/// **zero** sibling pairs (empty sibling_map) and must not produce any
/// `LoopDiscontinuity` from a corrupt sibling_map.
///
/// Before the fix, `pair_siblings` left candidates[j] un-matched so a later
/// iteration found exactly one remaining candidate and formed a spurious pair,
/// giving `pairs == 1` and a corrupt `sibling_map`.  The fix marks every
/// candidate j as matched inside the MultipleSiblings branch.
#[test]
fn four_way_sibling_clash_produces_no_spurious_pair() {
    let mut geom = GeomStore::new();
    let curve = geom.insert_curve(CurveGeom::Line(
        Line3::new(Point3::origin(), Vec3::Z).expect("valid line"),
    ));
    let surf_id = geom.insert_surface(SurfaceGeom::Plane(
        Plane::new(Point3::origin(), Vec3::Z).expect("valid plane"),
    ));

    let mut store = TopoStore::new();
    let p = geom.insert_point(VertexGeom::Explicit(Point3::origin()));
    let v = store.add_vertex(Vertex { point: p });

    // Four half-edges: two forward [0,1] and two reverse [1,0].
    let boundaries: [[f64; 2]; 4] = [
        [0.0_f64, 1.0_f64],
        [1.0_f64, 0.0_f64],
        [0.0_f64, 1.0_f64],
        [1.0_f64, 0.0_f64],
    ];
    let mut he_ids = Vec::new();
    for b in boundaries {
        he_ids.push(store.add_half_edge(HalfEdge {
            start: v,
            curve,
            boundary: b,
        }));
    }

    // Wrap the half-edges in a minimal loop/face/shell/solid so that
    // sibling_pairs (the public helper) and validate_topology can reach them.
    let lp = store.add_loop(Loop {
        half_edges: he_ids.clone(),
    });
    let face = store.add_face(Face {
        surface: surf_id,
        sense: Sense::Same,
        outer: lp,
        inners: vec![],
    });
    let shell = store.add_shell(Shell { faces: vec![face] });
    let solid = store.add_solid(Solid {
        shells: vec![shell],
    });

    let tol = Tol::default();

    // ── Primary assertion: no spurious pair ──────────────────────────────────
    //
    // sibling_pairs returns Err because every half-edge is unpaired.
    // The critical invariant: the returned Err must NOT contain any
    // LoopDiscontinuity that would indicate a corrupt sibling_map.  Before the
    // fix a spurious (1↔2) entry caused check_loop_continuity to mis-derive
    // end vertices.
    let sib_result = archi_kernel::topo::sibling_pairs(&store, &[solid], &tol);
    assert!(
        sib_result.is_err(),
        "sibling_pairs should fail for a 4-way clash (every half-edge is ambiguous)"
    );
    let sib_defects = sib_result.unwrap_err();

    // No sibling pair should have been formed → no LoopDiscontinuity sourced
    // from a corrupt sibling_map.
    assert!(
        !sib_defects
            .iter()
            .any(|d| matches!(d, Defect::LoopDiscontinuity { .. })),
        "corrupt sibling_map produced spurious LoopDiscontinuity: {sib_defects:?}"
    );

    // No MissingSibling that "should" have been MultipleSiblings, caused by
    // the spurious pairing consuming one candidate:
    // After the fix the greedy scan marks all candidates j=1,3 as consumed
    // when processing i=0's MultipleSiblings clash.  i=2 ([0,1]) then has no
    // remaining sibling candidates and is properly reported as MissingSibling
    // (both [1,0] partners were already consumed).  This is the correct
    // greedy outcome — the important property is that i=1 and i=3 are NOT
    // matched together spuriously, which would have inflated pairs by 1.
    //
    // Assert: no MultipleSiblings defect that simultaneously has pairs > 0
    // is present.  We check this indirectly via validate_topology.

    // ── Secondary: validate_topology must not report GenusMismatch from ──────
    //    wrong E count, and must not report LoopDiscontinuity.
    let v_result = validate_topology(&store, &[solid], &tol, None);
    let v_errs =
        v_result.expect_err("4-way clash must produce at least one defect in validate_topology");

    assert!(
        !v_errs
            .iter()
            .any(|d| matches!(d, Defect::GenusMismatch { .. })),
        "GenusMismatch detected: spurious pair likely inflated E count: {v_errs:?}"
    );
    assert!(
        !v_errs
            .iter()
            .any(|d| matches!(d, Defect::LoopDiscontinuity { .. })),
        "LoopDiscontinuity detected: corrupt sibling_map likely produced it: {v_errs:?}"
    );

    // All defects reported must be either MultipleSiblings or MissingSibling.
    // (MissingSibling for the remainder after greedy candidate consumption is
    // correct and expected; the forbidden outcome is LoopDiscontinuity /
    // GenusMismatch from a corrupt/inflated pairing.)
    assert!(
        v_errs.iter().all(|d| matches!(
            d,
            Defect::MultipleSiblings { .. } | Defect::MissingSibling { .. }
        )),
        "unexpected defect kinds in 4-way clash: {v_errs:?}"
    );
}

// ── Bug 2: empty Shell / Solid passes validate_topology ───────────────────────

/// A `Shell` with `faces = []` must be rejected.  Before the fix the Euler
/// check saw `S=1, V=E=F=L=0 → two_g=2` (non-negative, even), so it silently
/// passed.
#[test]
fn empty_shell_is_rejected() {
    let mut store = TopoStore::new();
    let shell_id = store.add_shell(Shell { faces: vec![] });
    let solid_id = store.add_solid(Solid {
        shells: vec![shell_id],
    });

    let tol = Tol::default();
    let result = validate_topology(&store, &[solid_id], &tol, None);

    assert!(
        result.is_err(),
        "empty shell (faces=[]) must be rejected by validate_topology, got Ok(())"
    );
    let errs = result.unwrap_err();
    assert!(
        errs.iter().any(|d| matches!(d, Defect::EmptyShell { .. })),
        "expected EmptyShell defect, got {errs:?}"
    );
}

/// A `Solid` with `shells = []` must be rejected.  Before the fix the Euler
/// check saw `S=0, two_g=0` (genus 0), so it silently passed.
#[test]
fn empty_solid_is_rejected() {
    let mut store = TopoStore::new();
    let solid_id = store.add_solid(Solid { shells: vec![] });

    let tol = Tol::default();
    let result = validate_topology(&store, &[solid_id], &tol, None);

    assert!(
        result.is_err(),
        "empty solid (shells=[]) must be rejected by validate_topology, got Ok(())"
    );
    let errs = result.unwrap_err();
    assert!(
        errs.iter().any(|d| matches!(d, Defect::EmptySolid { .. })),
        "expected EmptySolid defect, got {errs:?}"
    );
}
