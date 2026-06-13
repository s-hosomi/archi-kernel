use archi_kernel::curved::{
    tessellate_cone_panel, ConePanel, ConePanelOptions, CurvedError, TrimLoop2d,
};
use archi_kernel::{Point3, Tol, Unit3};

fn cone_panel(holes: Vec<TrimLoop2d>) -> ConePanel {
    let tol = Tol::default();
    ConePanel::new(
        Point3::origin(),
        Unit3::Z,
        0.35_f64,
        0.0_f64,
        1.0_f64,
        1.0_f64,
        2.2_f64,
        holes,
        &tol,
    )
    .expect("cone panel")
}

#[test]
fn cone_panel_without_holes_tessellates() {
    let tol = Tol::default();
    let panel = cone_panel(Vec::new());
    let mesh = tessellate_cone_panel(
        &panel,
        &ConePanelOptions::with_chord_tolerance(1e-3_f64),
        &tol,
    )
    .expect("mesh");

    assert!(mesh.vertex_count() > 0);
    assert!(mesh.triangle_count() > 0);
    let exact = panel.untrimmed_surface_area();
    assert!(
        (mesh.surface_area() - exact).abs() <= 8e-3_f64,
        "mesh area {} vs exact {exact}",
        mesh.surface_area()
    );
}

#[test]
fn rectangular_uv_hole_removes_conical_area() {
    let tol = Tol::default();
    let hole = TrimLoop2d::rectangle(0.25_f64, 0.55_f64, 1.3_f64, 1.7_f64, &tol)
        .expect("hole")
        .reversed();
    let panel = cone_panel(vec![hole]);
    let mesh = tessellate_cone_panel(
        &panel,
        &ConePanelOptions::with_chord_tolerance(1e-3_f64),
        &tol,
    )
    .expect("mesh");

    let half_angle = 0.35_f64;
    let hole_area = 0.5_f64 * half_angle.tan() / half_angle.cos()
        * (1.7_f64 * 1.7_f64 - 1.3_f64 * 1.3_f64)
        * (0.55_f64 - 0.25_f64);
    let exact = panel.untrimmed_surface_area() - hole_area;
    assert!(
        (mesh.surface_area() - exact).abs() <= 8e-3_f64,
        "mesh area {} vs exact {exact}",
        mesh.surface_area()
    );
}

#[test]
fn cone_panel_rejects_apex_crossing_domain() {
    let tol = Tol::default();
    let err = ConePanel::new(
        Point3::origin(),
        Unit3::Z,
        0.35_f64,
        0.0_f64,
        1.0_f64,
        0.0_f64,
        2.0_f64,
        Vec::new(),
        &tol,
    )
    .expect_err("apex must be rejected");
    assert!(matches!(err, CurvedError::ApexCrossing));
}

#[test]
fn cone_panel_rejects_invalid_half_angle() {
    let tol = Tol::default();
    let err = ConePanel::new(
        Point3::origin(),
        Unit3::Z,
        std::f64::consts::FRAC_PI_2,
        0.0_f64,
        1.0_f64,
        1.0_f64,
        2.0_f64,
        Vec::new(),
        &tol,
    )
    .expect_err("angle must be rejected");
    assert!(matches!(err, CurvedError::InvalidConeAngle { .. }));
}
