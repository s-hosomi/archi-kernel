//! Phase 2 extrusion tests: counts, watertight `validate(Full)`, and analytic
//! volume checks (the orientation test). Every literal carries an `f64`
//! annotation and an explicit tolerance, per `DESIGN.md` §12.

use std::f64::consts::PI;

use archi_kernel::build::extrude;
use archi_kernel::csg::{CsgNode, EvalError, Member, Profile2d};
use archi_kernel::geom::SurfaceGeom;
use archi_kernel::topo::validate::validate_topology;
use archi_kernel::{Line3, Point3, Tol, ValidateLevel, Vec3};

/// Volume tolerance: well above the `f64` round-off accumulated over a few
/// dozen face integrals at building scale, yet far tighter than any real error.
const VOL_EPS: f64 = 1e-9;

fn z_axis() -> Line3 {
    Line3::new(Point3::origin(), Vec3::Z).expect("valid axis")
}

// ── Rectangle ────────────────────────────────────────────────────────────────

#[test]
fn rect_extrusion_counts_validate_and_volume() {
    let tol = Tol::default();
    // width w = 0.3, depth d = 0.6 ⇒ half_w = 0.15, half_h = 0.3.
    let half_w = 0.15_f64;
    let half_h = 0.3_f64;
    let length = 3.0_f64;
    let profile = Profile2d::rect(half_w, half_h).expect("valid rect");
    let brep = extrude(&profile, &z_axis(), length, &tol).expect("extrude rect");

    // V8 / E12 (24 half-edges) / F6.
    assert_eq!(brep.topo.vertices.len(), 8usize, "vertices");
    assert_eq!(brep.topo.half_edges.len(), 24usize, "half-edges");
    assert_eq!(brep.topo.faces.len(), 6usize, "faces");

    validate_topology(&brep.topo, &brep.solids, &tol, Some(0u32)).expect("genus 0");
    brep.validate(&tol, ValidateLevel::Full)
        .expect("rect extrusion full");

    // Volume = w · d · L.
    let w = 2.0_f64 * half_w;
    let d = 2.0_f64 * half_h;
    let expected = w * d * length;
    let v = brep.signed_volume();
    assert!(
        v > 0.0_f64,
        "outward orientation ⇒ positive volume, got {v}"
    );
    assert!(
        (v - expected).abs() < VOL_EPS,
        "v = {v}, expected {expected}"
    );

    // Canonicalisation must de-duplicate coincident planes: a box has 6 faces
    // but only 3 distinct plane orientations × ... actually opposite faces are
    // 1 m apart so are distinct planes; the point is no plane appears twice, so
    // the distinct-plane count equals the face count (6), never more.
    let surfaces = referenced_plane_count(&brep);
    assert_eq!(
        surfaces, 6usize,
        "each cap/side has its own canonical plane"
    );
}

/// Count how many distinct plane `SurfaceId`s are referenced by the faces.
fn referenced_plane_count(brep: &archi_kernel::Brep) -> usize {
    use std::collections::HashSet;
    let mut ids = HashSet::new();
    for solid_id in &brep.solids {
        let solid = brep.topo.solids.get(*solid_id).unwrap();
        for shell_id in &solid.shells {
            let shell = brep.topo.shells.get(*shell_id).unwrap();
            for face_id in &shell.faces {
                let face = brep.topo.faces.get(*face_id).unwrap();
                if let Some(SurfaceGeom::Plane(_)) = brep.geom.surface(face.surface) {
                    ids.insert(face.surface);
                }
            }
        }
    }
    ids.len()
}

// ── H-section ────────────────────────────────────────────────────────────────

#[test]
fn h_section_extrusion_counts_validate_and_volume() {
    let tol = Tol::default();
    // b = 0.2, h = 0.4, tw = 0.01, tf = 0.02.
    let half_w = 0.1_f64;
    let half_h = 0.2_f64;
    let web = 0.01_f64;
    let flange = 0.02_f64;
    let length = 5.0_f64;
    let profile = Profile2d::h_section(half_w, half_h, web, flange).expect("valid H");
    let brep = extrude(&profile, &z_axis(), length, &tol).expect("extrude H");

    // 12-gon ⇒ V24 / E36 (72 half-edges) / F14 (2 caps + 12 sides).
    assert_eq!(brep.topo.vertices.len(), 24usize, "vertices");
    assert_eq!(brep.topo.half_edges.len(), 72usize, "half-edges");
    assert_eq!(brep.topo.faces.len(), 14usize, "faces");

    validate_topology(&brep.topo, &brep.solids, &tol, Some(0u32)).expect("genus 0");
    brep.validate(&tol, ValidateLevel::Full)
        .expect("H extrusion full");

    // Volume = area · L, area = 2·b·tf + (h − 2·tf)·tw.
    let b = 2.0_f64 * half_w;
    let h = 2.0_f64 * half_h;
    let area = 2.0_f64 * b * flange + (h - 2.0_f64 * flange) * web;
    let expected = area * length;
    let v = brep.signed_volume();
    assert!(v > 0.0_f64, "positive volume, got {v}");
    assert!(
        (v - expected).abs() < VOL_EPS,
        "v = {v}, expected {expected}"
    );
}

// ── Circle ───────────────────────────────────────────────────────────────────

#[test]
fn circle_extrusion_validates_and_volume() {
    let tol = Tol::default();
    let radius = 0.3_f64;
    let length = 4.0_f64;
    let profile = Profile2d::circle(radius).expect("valid circle");
    let brep = extrude(&profile, &z_axis(), length, &tol).expect("extrude circle");

    // Seam split ⇒ V4 / E6 (12 half-edges) / F4 (2 caps + 2 half-cylinders).
    assert_eq!(brep.topo.vertices.len(), 4usize, "vertices");
    assert_eq!(brep.topo.half_edges.len(), 12usize, "half-edges");
    assert_eq!(brep.topo.faces.len(), 4usize, "faces");

    // One shell, genus 0.
    validate_topology(&brep.topo, &brep.solids, &tol, Some(0u32)).expect("genus 0, 1 shell");
    let solid = brep.topo.solids.get(brep.solids[0]).unwrap();
    assert_eq!(solid.shells.len(), 1usize, "single shell");

    brep.validate(&tol, ValidateLevel::Full)
        .expect("circle extrusion full");

    // Volume = π · r² · L.
    let expected = PI * radius * radius * length;
    let v = brep.signed_volume();
    assert!(v > 0.0_f64, "positive volume, got {v}");
    assert!(
        (v - expected).abs() < VOL_EPS,
        "v = {v}, expected {expected}"
    );
}

// ── Oblique axis ─────────────────────────────────────────────────────────────

#[test]
fn oblique_rect_extrusion_volume_and_validate() {
    let tol = Tol::default();
    let half_w = 0.2_f64;
    let half_h = 0.25_f64;
    let length = 2.5_f64;
    let axis = Line3::new(
        Point3::new(1.0_f64, -2.0_f64, 0.5_f64),
        Vec3::new(1.0_f64, 1.0_f64, 1.0_f64),
    )
    .expect("valid oblique axis");
    let profile = Profile2d::rect(half_w, half_h).expect("valid rect");
    let brep = extrude(&profile, &axis, length, &tol).expect("extrude oblique");

    brep.validate(&tol, ValidateLevel::Full)
        .expect("oblique extrusion full");

    let w = 2.0_f64 * half_w;
    let d = 2.0_f64 * half_h;
    let expected = w * d * length;
    let v = brep.signed_volume();
    assert!(v > 0.0_f64, "positive volume, got {v}");
    assert!(
        (v - expected).abs() < VOL_EPS,
        "v = {v}, expected {expected}"
    );
}

#[test]
fn oblique_circle_extrusion_volume_and_validate() {
    let tol = Tol::default();
    let radius = 0.25_f64;
    let length = 3.0_f64;
    let axis = Line3::new(
        Point3::new(-1.0_f64, 0.5_f64, 2.0_f64),
        Vec3::new(2.0_f64, -1.0_f64, 3.0_f64),
    )
    .expect("valid oblique axis");
    let profile = Profile2d::circle(radius).expect("valid circle");
    let brep = extrude(&profile, &axis, length, &tol).expect("extrude oblique circle");

    brep.validate(&tol, ValidateLevel::Full)
        .expect("oblique circle full");

    let expected = PI * radius * radius * length;
    let v = brep.signed_volume();
    assert!(
        (v - expected).abs() < VOL_EPS,
        "v = {v}, expected {expected}"
    );
}

// ── Invalid input ────────────────────────────────────────────────────────────

#[test]
fn non_positive_length_is_rejected() {
    let tol = Tol::default();
    let profile = Profile2d::rect(0.1_f64, 0.1_f64).expect("valid rect");
    assert!(extrude(&profile, &z_axis(), 0.0_f64, &tol).is_err());
    assert!(extrude(&profile, &z_axis(), -1.0_f64, &tol).is_err());
}

#[test]
fn degenerate_h_section_is_rejected() {
    let tol = Tol::default();
    // web (0.3) > flange width b (0.2): degenerate.
    let wide_web = Profile2d::h_section(0.1_f64, 0.2_f64, 0.3_f64, 0.02_f64).expect("ctor ok");
    assert!(extrude(&wide_web, &z_axis(), 1.0_f64, &tol).is_err());

    // 2·tf (0.4) = h (0.4): flanges meet, degenerate.
    let thick_flange = Profile2d::h_section(0.1_f64, 0.2_f64, 0.01_f64, 0.2_f64).expect("ctor ok");
    assert!(extrude(&thick_flange, &z_axis(), 1.0_f64, &tol).is_err());
}

// ── CSG path ─────────────────────────────────────────────────────────────────

#[test]
fn csg_extrude_evaluates_and_reevaluates_when_dirty() {
    let tol = Tol::default();
    let profile = Profile2d::rect(0.15_f64, 0.3_f64).expect("valid rect");
    let mut member = Member::new(CsgNode::Extrude {
        origin: Point3::origin(),
        profile,
        axis: Vec3::Z,
        length: 3.0_f64,
    });

    // First evaluation succeeds and caches.
    {
        let brep = member.brep(&tol).expect("first eval");
        let v = brep.signed_volume();
        let expected = 0.3_f64 * 0.6_f64 * 3.0_f64;
        assert!((v - expected).abs() < VOL_EPS, "v = {v}");
    }
    assert!(member.last_valid().is_some());
    assert!(!member.is_dirty(&tol), "clean after eval");

    // Mutate the tree and mark dirty: the next eval must rebuild with the new
    // length.
    member.csg = CsgNode::Extrude {
        origin: Point3::origin(),
        profile,
        axis: Vec3::Z,
        length: 6.0_f64,
    };
    member.mark_dirty();
    assert!(member.is_dirty(&tol));

    let brep = member.brep(&tol).expect("re-eval");
    let v = brep.signed_volume();
    let expected = 0.3_f64 * 0.6_f64 * 6.0_f64;
    assert!((v - expected).abs() < VOL_EPS, "re-eval v = {v}");
}

#[test]
fn csg_extrude_rejects_bad_length() {
    let tol = Tol::default();
    let profile = Profile2d::rect(0.1_f64, 0.1_f64).expect("valid rect");
    let mut member = Member::new(CsgNode::Extrude {
        origin: Point3::origin(),
        profile,
        axis: Vec3::Z,
        length: -1.0_f64,
    });
    assert!(matches!(member.brep(&tol), Err(EvalError::Construction(_))));
    assert!(member.last_valid().is_none());
}
