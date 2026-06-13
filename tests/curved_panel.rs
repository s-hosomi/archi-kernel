use archi_kernel::curved::{
    tessellate_cylinder_panel, CurvedError, CylinderPanel, CylinderPanelOptions, TrimLoop2d,
};
use archi_kernel::{Cylinder, Line3, Point3, Tol, Vec3};

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
