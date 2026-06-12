//! Phase 5 mass-property unit tests: centroid, surface area, formwork split,
//! and the elliptical-rim volume (`DESIGN.md` §6-4).

use archi_kernel::build::extrude;
use archi_kernel::csg::Profile2d;
use archi_kernel::mass::{centroid, formwork_area, surface_area};
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::Line3;
use archi_kernel::tolerance::Tol;

const EPS: f64 = 1e-9;

/// A box `w(x) × d(y) × h(z)` centred in plan at the origin, base at z = 0,
/// extruded along +z. The +z profile axes are (u, v) = (Y, −X): half_w = Y half,
/// half_h = X half.
fn box_brep(w: f64, d: f64, h: f64) -> archi_kernel::Brep {
    let tol = Tol::default();
    let profile = Profile2d::rect(d / 2.0, w / 2.0).expect("rect");
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    extrude(&profile, &line, h, &tol).expect("box")
}

#[test]
fn box_centroid_is_geometric_centre() {
    let tol = Tol::default();
    let brep = box_brep(2.0, 1.0, 4.0);
    let c = centroid(&brep, &tol).expect("centroid");
    // Plan-centred at origin, base z = 0, height 4 ⇒ centroid (0, 0, 2).
    assert!(c.x.abs() <= EPS, "cx = {}", c.x);
    assert!(c.y.abs() <= EPS, "cy = {}", c.y);
    assert!((c.z - 2.0_f64).abs() <= EPS, "cz = {}", c.z);
}

#[test]
fn box_surface_area_is_sum_of_faces() {
    let brep = box_brep(2.0, 3.0, 4.0);
    let a = surface_area(&brep).expect("surface area");
    // 2(wd + wh + dh) = 2(6 + 8 + 12) = 52.
    let expected = 2.0_f64 * (2.0 * 3.0 + 2.0 * 4.0 + 3.0 * 4.0);
    assert!(
        (a - expected).abs() <= EPS,
        "area = {a}, expected {expected}"
    );
}

#[test]
fn box_formwork_splits_side_and_bottom() {
    let tol = Tol::default();
    let brep = box_brep(2.0, 3.0, 4.0);
    let fw = formwork_area(&brep, &tol).expect("formwork");
    // Sides (the 4 vertical faces): 2·(w·h) + 2·(d·h) = 2·8 + 2·12 = 40.
    let side_expected = 2.0_f64 * (2.0 * 4.0) + 2.0_f64 * (3.0 * 4.0);
    // Bottom (−z face): w·d = 6. Top (+z) carries no formwork.
    let bottom_expected = 2.0_f64 * 3.0;
    assert!(
        (fw.side - side_expected).abs() <= EPS,
        "side = {}, expected {side_expected}",
        fw.side
    );
    assert!(
        (fw.bottom - bottom_expected).abs() <= EPS,
        "bottom = {}, expected {bottom_expected}",
        fw.bottom
    );
}

#[test]
fn round_column_volume_and_lateral_area() {
    use std::f64::consts::PI;
    let tol = Tol::default();
    let r = 0.3_f64;
    let h = 3.0_f64;
    let profile = Profile2d::circle(r).expect("circle");
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    let brep = extrude(&profile, &line, h, &tol).expect("column");

    // Volume π r² h.
    let v = brep.signed_volume();
    assert!((v - PI * r * r * h).abs() <= 1e-6, "V = {v}");

    // Surface area: lateral 2πrh + two caps 2·πr².
    let a = surface_area(&brep).expect("area");
    let expected = 2.0 * PI * r * h + 2.0 * PI * r * r;
    assert!(
        (a - expected).abs() <= 1e-6,
        "area = {a}, expected {expected}"
    );

    // Formwork: the round wall is vertical ⇒ all lateral area is side formwork;
    // the bottom cap (−z) is bottom formwork; the top cap carries none.
    let fw = formwork_area(&brep, &tol).expect("formwork");
    assert!(
        (fw.side - 2.0 * PI * r * h).abs() <= 1e-6,
        "side = {}, expected {}",
        fw.side,
        2.0 * PI * r * h
    );
    assert!(
        (fw.bottom - PI * r * r).abs() <= 1e-6,
        "bottom = {}, expected {}",
        fw.bottom,
        PI * r * r
    );
}
