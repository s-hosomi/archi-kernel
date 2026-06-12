//! Phase 6 — watertight tessellation tests (`DESIGN.md` §6-5, §7, §10 Phase 6).
//!
//! Targets the member shapes the kernel produces: a cube and an H-prism
//! (planar, with concavity), a vertical and an obliquely-cut cylinder
//! (cylinder patches, circular and elliptical rims), a wall pierced by a
//! circular sleeve and a box with a rectangular through-hole (annulus caps),
//! and the mini-building members from the section-drawing fixture.
//!
//! Each case checks three things:
//!
//! 1. **Watertight** — every interior mesh edge `(vi, vj)` appears exactly twice,
//!    once in each direction (the `DESIGN.md` §7 invariant: 各エッジちょうど 2
//!    三角形). The checker is the test utility `assert_watertight`.
//! 2. **Orientation** — the mesh's signed volume is positive (outward normals).
//! 3. **Volume** — the mesh volume matches the analytic / `signed_volume_checked`
//!    volume to a tolerance derived from the chord error (exact for purely
//!    planar solids, a chord-order systematic deficit for cylinders).
//!
//! Every literal carries an `f64` annotation and an explicit tolerance
//! (`DESIGN.md` §12).

use std::collections::HashMap;
use std::f64::consts::PI;

use archi_kernel::boolean::prismatic::{self, ExtrudeLeaf};
use archi_kernel::brep::Brep;
use archi_kernel::build::extrude;
use archi_kernel::csg::Profile2d;
use archi_kernel::mass::signed_volume_checked;
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::{Line3, Plane};
use archi_kernel::tess::{tessellate, Mesh, TessOptions};
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::ValidateLevel;

// ── test utility: watertight check ────────────────────────────────────────────

/// Assert the mesh is watertight: every directed edge `(a, b)` of every triangle
/// has its reverse `(b, a)` present exactly once, and no edge is shared by more
/// than two triangles. Equivalent to "each undirected edge is used by exactly
/// two triangles with opposite orientation" (`DESIGN.md` §7).
fn assert_watertight(mesh: &Mesh) {
    assert!(mesh.triangle_count() > 0, "mesh has no triangles");
    // Count directed edges.
    let mut dir: HashMap<(u32, u32), i32> = HashMap::new();
    for k in 0..mesh.triangle_count() {
        let t = [
            mesh.indices[3 * k],
            mesh.indices[3 * k + 1],
            mesh.indices[3 * k + 2],
        ];
        assert!(
            t[0] != t[1] && t[1] != t[2] && t[2] != t[0],
            "degenerate triangle {t:?}"
        );
        for e in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *dir.entry(e).or_insert(0) += 1;
        }
    }
    // Each directed edge appears at most once (manifold), and its reverse must
    // appear exactly once (boundary-free).
    for (&(a, b), &count) in &dir {
        assert_eq!(
            count, 1,
            "directed edge ({a},{b}) used {count} times (non-manifold)"
        );
        let rev = dir.get(&(b, a)).copied().unwrap_or(0);
        assert_eq!(
            rev, 1,
            "edge ({a},{b}) has no matching reverse ({b},{a}) — open boundary"
        );
    }
}

/// The default display options.
fn opts() -> TessOptions {
    TessOptions::default()
}

/// Tessellate a single-solid brep's first solid, asserting success.
fn tess(brep: &Brep, o: &TessOptions, tol: &Tol) -> Mesh {
    tessellate(brep, brep.solids[0], o, tol).expect("tessellation")
}

// ── cube (planar, exact volume) ───────────────────────────────────────────────

#[test]
fn cube_is_watertight_and_exact() {
    let tol = Tol::default();
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    let cube = extrude(
        &Profile2d::rect(0.5_f64, 0.5_f64).expect("rect"),
        &line,
        1.0_f64,
        &tol,
    )
    .expect("cube");
    cube.validate(&tol, ValidateLevel::Full)
        .expect("valid cube");

    let mesh = tess(&cube, &opts(), &tol);
    assert_watertight(&mesh);
    assert!(mesh.signed_volume() > 0.0_f64, "outward orientation");
    // Planar solid ⇒ mesh volume is exact.
    let analytic = 1.0_f64 * 1.0_f64 * 1.0_f64;
    assert!(
        (mesh.signed_volume() - analytic).abs() < 1e-9_f64,
        "cube volume {} (expected {analytic})",
        mesh.signed_volume()
    );
    // face_of tags every triangle to a real face index.
    assert_eq!(mesh.face_of.len(), mesh.triangle_count());
}

// ── H-prism (planar, concave) ─────────────────────────────────────────────────

#[test]
fn h_prism_is_watertight_and_exact() {
    let tol = Tol::default();
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    // H-section: half_w 0.15, half_h 0.10, web 0.02, flange 0.02.
    let prof = Profile2d::h_section(0.15_f64, 0.10_f64, 0.02_f64, 0.02_f64).expect("h");
    let h = extrude(&prof, &line, 0.5_f64, &tol).expect("h prism");
    h.validate(&tol, ValidateLevel::Full)
        .expect("valid h prism");

    let mesh = tess(&h, &opts(), &tol);
    assert_watertight(&mesh);
    assert!(mesh.signed_volume() > 0.0_f64);
    // Compare against the kernel's own closed-form volume (exact for planar).
    let analytic = signed_volume_checked(&h).expect("closed-form volume");
    assert!(
        (mesh.signed_volume() - analytic).abs() < 1e-9_f64,
        "H-prism mesh volume {} vs analytic {analytic}",
        mesh.signed_volume()
    );
}

// ── vertical cylinder (circular rims) ─────────────────────────────────────────

#[test]
fn vertical_cylinder_is_watertight_and_converges() {
    let tol = Tol::default();
    let r = 0.3_f64;
    let len = 1.2_f64;
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    let cyl = extrude(&Profile2d::circle(r).expect("circle"), &line, len, &tol).expect("cylinder");
    cyl.validate(&tol, ValidateLevel::Full)
        .expect("valid cylinder");

    let analytic = PI * r * r * len;

    // Coarse and fine chord tolerances.
    let coarse = TessOptions::with_chord_tolerance(1e-2_f64);
    // 1/4 of coarse ⇒ ~1/2 the chord step.
    let fine = TessOptions::with_chord_tolerance(2.5e-3_f64);

    let m_coarse = tess(&cyl, &coarse, &tol);
    let m_fine = tess(&cyl, &fine, &tol);

    assert_watertight(&m_coarse);
    assert_watertight(&m_fine);
    assert!(m_coarse.signed_volume() > 0.0_f64, "outward orientation");

    // A polygonal approximation of a disk *under*-estimates its area, so the
    // mesh volume is below the analytic volume; the deficit is O(chord²).
    let err_coarse = (analytic - m_coarse.signed_volume()).abs();
    let err_fine = (analytic - m_fine.signed_volume()).abs();
    assert!(
        err_fine < err_coarse,
        "refining chord tolerance must reduce the volume error: coarse {err_coarse}, fine {err_fine}"
    );
    // Halving the chord step quarters the error (the inscribed-polygon area
    // deficit is quadratic in the segment angle). Allow generous slack.
    assert!(
        err_fine < err_coarse * 0.5_f64,
        "error should fall ~4× when the step halves: coarse {err_coarse}, fine {err_fine}"
    );

    // The fine mesh is within a chord-order tolerance of the analytic volume.
    // Systematic deficit bound: ≈ (1/2) L r · sinΔφ·… ; use a chord-derived
    // budget. With chord_tol 2.5e-3 and r 0.3, Δφ ≈ 2√(2·chord/r) ≈ 0.26 rad,
    // the relative area deficit ≈ Δφ²/24 ≈ 3e-3, so allow 1%.
    assert!(
        err_fine < analytic * 1e-2_f64,
        "fine cylinder volume {} within 1% of {analytic} (err {err_fine})",
        m_fine.signed_volume()
    );
}

// ── obliquely-cut cylinder (elliptical rim) ───────────────────────────────────

#[test]
fn oblique_cut_cylinder_is_watertight() {
    use archi_kernel::boolean::{cut, KeepSide};
    let tol = Tol::default();
    let r = 0.3_f64;
    let len = 3.0_f64;
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    let cyl = extrude(&Profile2d::circle(r).expect("circle"), &line, len, &tol).expect("cylinder");

    // Oblique plane through the axis at mid-height (0,0,1.5), normal (0,1,1):
    // tilts the section into an ellipse and keeps exactly half the cylinder
    // (the known-good oblique cut, mirroring `tests/half_space.rs`). Keeping
    // `Below` leaves the side patches with one circular bottom rim and one
    // elliptical cut rim.
    let cut_plane = Plane::new(
        Point3::new(0.0_f64, 0.0_f64, 1.5_f64),
        Vec3::new(0.0_f64, 1.0_f64, 1.0_f64),
    )
    .expect("oblique plane");
    let res = cut(&cyl, cyl.solids[0], &cut_plane, KeepSide::Below, &tol).expect("cut");
    let chopped = res.brep();
    chopped
        .validate(&tol, ValidateLevel::Full)
        .expect("valid oblique cut");

    let o = TessOptions::with_chord_tolerance(2e-3_f64);
    let mesh = tess(&chopped, &o, &tol);
    assert_watertight(&mesh);
    assert!(mesh.signed_volume() > 0.0_f64, "outward orientation");

    // Compare with the kernel's closed-form volume (handles the oblique patch +
    // ellipse cap). Tolerance is chord-order on the curved part.
    let analytic = signed_volume_checked(&chopped).expect("closed-form volume");
    assert!(
        (mesh.signed_volume() - analytic).abs() < analytic * 2e-2_f64,
        "oblique cylinder mesh volume {} vs analytic {analytic}",
        mesh.signed_volume()
    );
}

// ── wall − circular sleeve (annulus on the cylinder wall, prismatic result) ───

#[test]
fn wall_with_circular_sleeve_is_watertight() {
    let tol = Tol::default();
    let wall = ExtrudeLeaf {
        profile: Profile2d::rect(1.0_f64, 0.1_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 0.1_f64, 1.0_f64),
        axis: Vec3::X,
        length: 3.0_f64,
    };
    let r = 0.05_f64;
    let sleeve = ExtrudeLeaf {
        profile: Profile2d::circle(r).expect("circle"),
        origin: Point3::new(1.5_f64, -0.1_f64, 1.0_f64),
        axis: Vec3::Y,
        length: 0.4_f64,
    };
    let result = prismatic::difference(&wall, &sleeve, &tol).expect("sleeve difference");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight wall with circular hole");

    let o = TessOptions::with_chord_tolerance(2e-3_f64);
    let mesh = tess(&result, &o, &tol);
    assert_watertight(&mesh);
    assert!(mesh.signed_volume() > 0.0_f64);

    // V = V_wall − πr²·t, t = 0.2 (wall thickness). Curved hole ⇒ chord-order.
    let analytic = 3.0_f64 * 0.2_f64 * 2.0_f64 - PI * r * r * 0.2_f64;
    assert!(
        (mesh.signed_volume() - analytic).abs() < 1e-3_f64,
        "wall+sleeve mesh volume {} vs analytic {analytic}",
        mesh.signed_volume()
    );
}

// ── box with a rectangular through-hole (planar annulus caps) ─────────────────

#[test]
fn box_with_rectangular_through_hole_is_watertight_and_exact() {
    let tol = Tol::default();
    let base = ExtrudeLeaf {
        profile: Profile2d::rect(1.0_f64, 1.0_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 0.5_f64,
    };
    // A rectangular opening fully through the box along z.
    let hole = ExtrudeLeaf {
        profile: Profile2d::rect(0.3_f64, 0.3_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, -0.1_f64),
        axis: Vec3::Z,
        length: 0.7_f64,
    };
    let result = prismatic::difference(&base, &hole, &tol).expect("through-hole");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight through-hole box");

    let mesh = tess(&result, &opts(), &tol);
    assert_watertight(&mesh);
    assert!(mesh.signed_volume() > 0.0_f64);

    // Purely planar ⇒ exact. V = 2·2·0.5 − 0.6·0.6·0.5.
    let analytic = 2.0_f64 * 2.0_f64 * 0.5_f64 - 0.6_f64 * 0.6_f64 * 0.5_f64;
    assert!(
        (mesh.signed_volume() - analytic).abs() < 1e-9_f64,
        "through-hole mesh volume {} vs analytic {analytic}",
        mesh.signed_volume()
    );
}

// ── mini-building members (section-drawing fixture analogues) ──────────────────

#[test]
fn mini_building_members_each_tessellate_watertight() {
    let tol = Tol::default();

    // Rectangular column along +z.
    let rect_col = ExtrudeLeaf {
        profile: Profile2d::rect(0.2_f64, 0.3_f64).expect("rect"),
        origin: Point3::new(1.0_f64, 1.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 3.0_f64,
    };
    // Circular column along +z.
    let circ_col = ExtrudeLeaf {
        profile: Profile2d::circle(0.25_f64).expect("circle"),
        origin: Point3::new(2.0_f64, 3.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 3.0_f64,
    };
    // Girder along +x.
    let girder = ExtrudeLeaf {
        profile: Profile2d::rect(0.15_f64, 0.25_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 1.0_f64, 2.75_f64),
        axis: Vec3::X,
        length: 4.0_f64,
    };

    for leaf in [rect_col, girder] {
        let line = Line3::new(leaf.origin, leaf.axis).expect("axis");
        let brep = extrude(&leaf.profile, &line, leaf.length, &tol).expect("extrude");
        brep.validate(&tol, ValidateLevel::Full).expect("valid");
        let mesh = tess(&brep, &opts(), &tol);
        assert_watertight(&mesh);
        assert!(mesh.signed_volume() > 0.0_f64);
        // Planar members are exact.
        let analytic = signed_volume_checked(&brep).expect("volume");
        assert!(
            (mesh.signed_volume() - analytic).abs() < 1e-9_f64,
            "member mesh volume {} vs analytic {analytic}",
            mesh.signed_volume()
        );
    }

    // Circular column: watertight, outward, chord-order volume.
    let line = Line3::new(circ_col.origin, circ_col.axis).expect("axis");
    let brep = extrude(&circ_col.profile, &line, circ_col.length, &tol).expect("extrude");
    brep.validate(&tol, ValidateLevel::Full).expect("valid");
    let o = TessOptions::with_chord_tolerance(2e-3_f64);
    let mesh = tess(&brep, &o, &tol);
    assert_watertight(&mesh);
    assert!(mesh.signed_volume() > 0.0_f64);
    let analytic = signed_volume_checked(&brep).expect("volume");
    assert!(
        (mesh.signed_volume() - analytic).abs() < analytic * 1e-2_f64,
        "circular column mesh volume {} vs analytic {analytic}",
        mesh.signed_volume()
    );
}

// ── slab with a rectangular opening and a circular sleeve ─────────────────────

#[test]
fn slab_with_two_openings_is_watertight() {
    let tol = Tol::default();
    let base = ExtrudeLeaf {
        profile: Profile2d::rect(2.0_f64, 2.0_f64).expect("rect"),
        origin: Point3::new(2.0_f64, 2.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 0.2_f64,
    };
    let rect = ExtrudeLeaf {
        profile: Profile2d::rect(0.3_f64, 0.4_f64).expect("rect"),
        origin: Point3::new(1.0_f64, 1.0_f64, -0.1_f64),
        axis: Vec3::Z,
        length: 0.4_f64,
    };
    let circ = ExtrudeLeaf {
        profile: Profile2d::circle(0.3_f64).expect("circle"),
        origin: Point3::new(3.0_f64, 3.0_f64, -0.1_f64),
        axis: Vec3::Z,
        length: 0.4_f64,
    };
    let slab =
        prismatic::opening_subtraction(&base, &[rect, circ], &tol).expect("slab with openings");
    slab.validate(&tol, ValidateLevel::Full)
        .expect("valid slab");

    let o = TessOptions::with_chord_tolerance(2e-3_f64);
    let mesh = tess(&slab, &o, &tol);
    assert_watertight(&mesh);
    assert!(mesh.signed_volume() > 0.0_f64);

    // V = 4·4·0.2 − (0.6·0.8·0.2) − (π·0.3²·0.2). Curved sleeve ⇒ chord-order.
    let analytic =
        16.0_f64 * 0.2_f64 - 0.6_f64 * 0.8_f64 * 0.2_f64 - PI * 0.3_f64 * 0.3_f64 * 0.2_f64;
    assert!(
        (mesh.signed_volume() - analytic).abs() < 1e-3_f64,
        "slab mesh volume {} vs analytic {analytic}",
        mesh.signed_volume()
    );
}

// ── error surface: non-positive chord tolerance ───────────────────────────────

#[test]
fn non_positive_chord_tolerance_is_an_error() {
    use archi_kernel::tess::TessError;
    let tol = Tol::default();
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    let cube = extrude(
        &Profile2d::rect(0.5_f64, 0.5_f64).expect("rect"),
        &line,
        1.0_f64,
        &tol,
    )
    .expect("cube");
    let bad = TessOptions::with_chord_tolerance(0.0_f64);
    assert!(matches!(
        tessellate(&cube, cube.solids[0], &bad, &tol),
        Err(TessError::NonPositiveChordTolerance { .. })
    ));
}
