use archi_kernel::curved::{
    tessellate_sphere_panel, tessellate_thick_sphere_panel, CurvedError, SpherePanel,
    SpherePanelOptions, SurfaceMesh, ThickSpherePanel, TrimLoop2d,
};
use archi_kernel::{Point3, Tol, Unit3};
use std::collections::HashMap;

fn sphere_panel(holes: Vec<TrimLoop2d>) -> SpherePanel {
    let tol = Tol::default();
    SpherePanel::new(
        Point3::origin(),
        4.0_f64,
        Unit3::Z,
        0.0_f64,
        1.2_f64,
        -0.35_f64,
        0.45_f64,
        holes,
        &tol,
    )
    .expect("sphere panel")
}

#[test]
fn sphere_panel_without_holes_tessellates() {
    let tol = Tol::default();
    let panel = sphere_panel(Vec::new());
    let mesh = tessellate_sphere_panel(
        &panel,
        &SpherePanelOptions::with_chord_tolerance(5e-4_f64),
        &tol,
    )
    .expect("mesh");

    assert!(mesh.vertex_count() > 0);
    assert!(mesh.triangle_count() > 0);
    let exact = panel.untrimmed_surface_area();
    assert!(
        (mesh.surface_area() - exact).abs() <= 2e-2_f64,
        "mesh area {} vs exact {exact}",
        mesh.surface_area()
    );
}

#[test]
fn rectangular_uv_hole_removes_spherical_area() {
    let tol = Tol::default();
    let hole = TrimLoop2d::rectangle(0.35_f64, 0.65_f64, -0.1_f64, 0.2_f64, &tol)
        .expect("hole")
        .reversed();
    let panel = sphere_panel(vec![hole]);
    let mesh = tessellate_sphere_panel(
        &panel,
        &SpherePanelOptions::with_chord_tolerance(5e-4_f64),
        &tol,
    )
    .expect("mesh");

    let hole_area = 4.0_f64 * 4.0_f64 * (0.65_f64 - 0.35_f64) * (0.2_f64.sin() - (-0.1_f64).sin());
    let exact = panel.untrimmed_surface_area() - hole_area;
    assert!(
        (mesh.surface_area() - exact).abs() <= 2e-2_f64,
        "mesh area {} vs exact {exact}",
        mesh.surface_area()
    );
}

#[test]
fn sphere_panel_rejects_pole_crossing_domain() {
    let tol = Tol::default();
    let err = SpherePanel::new(
        Point3::origin(),
        4.0_f64,
        Unit3::Z,
        0.0_f64,
        1.0_f64,
        -std::f64::consts::FRAC_PI_2,
        0.2_f64,
        Vec::new(),
        &tol,
    )
    .expect_err("pole must be rejected");
    assert!(matches!(err, CurvedError::PoleCrossing));
}

#[test]
fn thick_sphere_panel_with_rectangular_hole_is_closed() {
    let tol = Tol::default();
    let hole = TrimLoop2d::rectangle(0.35_f64, 0.65_f64, -0.1_f64, 0.2_f64, &tol)
        .expect("hole")
        .reversed();
    let mid = sphere_panel(vec![hole]);
    let panel = ThickSpherePanel::new(mid, 0.2_f64).expect("thick");
    let mesh = tessellate_thick_sphere_panel(
        &panel,
        &SpherePanelOptions::with_chord_tolerance(5e-4_f64),
        &tol,
    )
    .expect("mesh");

    assert_mesh_edges_closed(&mesh);
}

#[test]
fn thick_sphere_panel_rejects_arc_holes_until_arc_sides_exist() {
    let tol = Tol::default();
    let hole = TrimLoop2d::circle([0.6_f64, 0.0_f64], 0.1_f64, &tol).expect("hole");
    let mid = sphere_panel(vec![hole]);
    let panel = ThickSpherePanel::new(mid, 0.2_f64).expect("thick");
    let err = tessellate_thick_sphere_panel(
        &panel,
        &SpherePanelOptions::with_chord_tolerance(5e-4_f64),
        &tol,
    )
    .expect_err("arc hole side walls not yet implemented");
    assert!(matches!(err, CurvedError::UnsupportedArcTrim));
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
