use archi_kernel::csg::{CsgNode, CurvedPanelNode, EvalError, Member, UnsupportedReason};
use archi_kernel::curved::{
    tessellate_cylinder_panel, tessellate_thick_cylinder_panel, CurvedError, CylinderPanel,
    CylinderPanelOptions, SurfaceMesh, ThickCylinderPanel, TrimLoop2d,
};
use archi_kernel::{Cylinder, Line3, Point3, Tol, Vec3};
use std::collections::HashMap;

fn cylinder() -> Cylinder {
    Cylinder::new(
        Line3::new(Point3::origin(), Vec3::Z).expect("axis"),
        2.0_f64,
    )
    .expect("cylinder")
}

#[test]
fn cylinder_panel_without_holes_tessellates() {
    let tol = Tol::default();
    let panel = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        std::f64::consts::FRAC_PI_2,
        0.0_f64,
        3.0_f64,
        Vec::new(),
        &tol,
    )
    .expect("panel");
    let mesh = tessellate_cylinder_panel(
        &panel,
        &CylinderPanelOptions::with_chord_tolerance(1e-4_f64),
        &tol,
    )
    .expect("mesh");

    assert!(mesh.vertex_count() > 0);
    assert!(mesh.triangle_count() > 0);
    let expected = 2.0_f64 * std::f64::consts::FRAC_PI_2 * 3.0_f64;
    assert!(
        (mesh.surface_area() - expected).abs() <= 2e-3_f64,
        "mesh area {} vs exact {expected}",
        mesh.surface_area()
    );
    assert!((panel.surface_area() - expected).abs() <= 1e-12_f64);
}

#[test]
fn rectangular_uv_hole_removes_panel_area() {
    let tol = Tol::default();
    let hole = TrimLoop2d::rectangle(0.4_f64, 0.8_f64, 1.0_f64, 2.0_f64, &tol)
        .expect("rectangular hole")
        .reversed();
    let panel = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.2_f64,
        0.0_f64,
        3.0_f64,
        vec![hole],
        &tol,
    )
    .expect("panel");
    let mesh = tessellate_cylinder_panel(
        &panel,
        &CylinderPanelOptions::with_chord_tolerance(5e-4_f64),
        &tol,
    )
    .expect("mesh");

    let exact = 2.0_f64 * ((1.2_f64 * 3.0_f64) - (0.4_f64 * 1.0_f64));
    assert!((panel.surface_area() - exact).abs() <= 1e-12_f64);
    assert!(
        (mesh.surface_area() - exact).abs() <= 5e-3_f64,
        "mesh area {} vs exact {exact}",
        mesh.surface_area()
    );
}

#[test]
fn hole_outside_panel_is_rejected() {
    let tol = Tol::default();
    let hole = TrimLoop2d::rectangle(0.8_f64, 1.4_f64, 1.0_f64, 2.0_f64, &tol).expect("hole");
    let err = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.2_f64,
        0.0_f64,
        3.0_f64,
        vec![hole],
        &tol,
    )
    .expect_err("outside hole");
    assert!(matches!(err, CurvedError::HoleOutsidePanel));
}

#[test]
fn circular_uv_hole_removes_panel_area() {
    let tol = Tol::default();
    let hole = TrimLoop2d::circle([0.7_f64, 1.5_f64], 0.2_f64, &tol).expect("circle");
    let panel = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.4_f64,
        0.0_f64,
        3.0_f64,
        vec![hole],
        &tol,
    )
    .expect("panel");
    let mesh = tessellate_cylinder_panel(
        &panel,
        &CylinderPanelOptions::with_chord_tolerance(1e-3_f64),
        &tol,
    )
    .expect("mesh");

    let exact = 2.0_f64 * (1.4_f64 * 3.0_f64 - std::f64::consts::PI * 0.2_f64 * 0.2_f64);
    assert!((panel.surface_area() - exact).abs() <= 1e-12_f64);
    assert!(
        (mesh.surface_area() - exact).abs() <= 4e-2_f64,
        "mesh area {} vs exact {exact}",
        mesh.surface_area()
    );
}

#[test]
fn overlapping_holes_are_rejected() {
    let tol = Tol::default();
    let a = TrimLoop2d::rectangle(0.2_f64, 0.7_f64, 0.5_f64, 1.5_f64, &tol).expect("a");
    let b = TrimLoop2d::rectangle(0.6_f64, 1.0_f64, 1.0_f64, 2.0_f64, &tol).expect("b");
    let err = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.2_f64,
        0.0_f64,
        3.0_f64,
        vec![a, b],
        &tol,
    )
    .expect_err("overlap");
    assert!(matches!(err, CurvedError::HoleOverlap));
}

#[test]
fn separate_circular_holes_with_overlapping_bounds_are_accepted() {
    let tol = Tol::default();
    let a = TrimLoop2d::circle([0.45_f64, 1.0_f64], 0.25_f64, &tol).expect("a");
    let b = TrimLoop2d::circle([0.82_f64, 1.37_f64], 0.25_f64, &tol).expect("b");
    // Bounding boxes overlap in both axes, but centre distance is greater than
    // the sum of radii, so the circular holes are genuinely disjoint.
    let panel = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.2_f64,
        0.0_f64,
        3.0_f64,
        vec![a, b],
        &tol,
    )
    .expect("disjoint circles must be accepted");
    assert_eq!(panel.holes.len(), 2);
}

#[test]
fn overlapping_circular_holes_are_rejected() {
    let tol = Tol::default();
    let a = TrimLoop2d::circle([0.45_f64, 1.0_f64], 0.25_f64, &tol).expect("a");
    let b = TrimLoop2d::circle([0.70_f64, 1.0_f64], 0.25_f64, &tol).expect("b");
    let err = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.2_f64,
        0.0_f64,
        3.0_f64,
        vec![a, b],
        &tol,
    )
    .expect_err("overlapping circles");
    assert!(matches!(err, CurvedError::HoleOverlap));
}

#[test]
fn thick_cylinder_panel_without_holes_is_closed_and_has_volume() {
    let tol = Tol::default();
    let mid = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.0_f64,
        0.0_f64,
        3.0_f64,
        Vec::new(),
        &tol,
    )
    .expect("mid");
    let panel = ThickCylinderPanel::new(mid, 0.2_f64).expect("thick");
    let mesh = tessellate_thick_cylinder_panel(
        &panel,
        &CylinderPanelOptions::with_chord_tolerance(5e-4_f64),
        &tol,
    )
    .expect("mesh");

    assert_mesh_edges_closed(&mesh);
    assert!(
        (mesh.signed_volume().abs() - panel.volume()).abs() <= 3e-3_f64,
        "mesh volume {} vs exact {}",
        mesh.signed_volume(),
        panel.volume()
    );
}

#[test]
fn thick_cylinder_panel_with_rectangular_hole_is_closed() {
    let tol = Tol::default();
    let hole = TrimLoop2d::rectangle(0.4_f64, 0.8_f64, 1.0_f64, 2.0_f64, &tol)
        .expect("hole")
        .reversed();
    let mid = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.2_f64,
        0.0_f64,
        3.0_f64,
        vec![hole],
        &tol,
    )
    .expect("mid");
    let panel = ThickCylinderPanel::new(mid, 0.2_f64).expect("thick");
    let mesh = tessellate_thick_cylinder_panel(
        &panel,
        &CylinderPanelOptions::with_chord_tolerance(5e-4_f64),
        &tol,
    )
    .expect("mesh");

    assert_mesh_edges_closed(&mesh);
    assert!(
        (mesh.signed_volume().abs() - panel.volume()).abs() <= 5e-3_f64,
        "mesh volume {} vs exact {}",
        mesh.signed_volume(),
        panel.volume()
    );
}

#[test]
fn thick_cylinder_panel_rejects_arc_holes_until_ruled_arc_sides_exist() {
    let tol = Tol::default();
    let hole = TrimLoop2d::circle([0.7_f64, 1.5_f64], 0.2_f64, &tol).expect("hole");
    let mid = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.4_f64,
        0.0_f64,
        3.0_f64,
        vec![hole],
        &tol,
    )
    .expect("mid");
    let panel = ThickCylinderPanel::new(mid, 0.2_f64).expect("thick");
    let err = tessellate_thick_cylinder_panel(
        &panel,
        &CylinderPanelOptions::with_chord_tolerance(5e-4_f64),
        &tol,
    )
    .expect_err("arc hole side walls not yet implemented");
    assert!(matches!(err, CurvedError::UnsupportedArcTrim));
}

#[test]
fn curved_panel_csg_node_is_explicitly_unsupported_by_brep_evaluator() {
    let tol = Tol::default();
    let mid = CylinderPanel::new(
        cylinder(),
        0.0_f64,
        1.0_f64,
        0.0_f64,
        3.0_f64,
        Vec::new(),
        &tol,
    )
    .expect("mid");
    let panel = ThickCylinderPanel::new(mid, 0.2_f64).expect("thick");
    let mut member = Member::new(CsgNode::CurvedPanel(CurvedPanelNode { panel }));

    let err = member.brep(&tol).expect_err("curved B-rep is unsupported");
    assert!(matches!(
        err,
        EvalError::Unsupported3dBoolean {
            reason: UnsupportedReason::CurvedPanel
        }
    ));
}

fn assert_mesh_edges_closed(mesh: &SurfaceMesh) {
    let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
    for tri in mesh.indices.chunks_exact(3) {
        for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if a <= b { (a, b) } else { (b, a) };
            *counts.entry(key).or_default() += 1;
        }
    }
    assert!(
        counts.values().all(|&n| n == 2),
        "non-closed edges: {:?}",
        counts
            .iter()
            .filter_map(|(k, &n)| (n != 2).then_some((*k, n)))
            .collect::<Vec<_>>()
    );
}
