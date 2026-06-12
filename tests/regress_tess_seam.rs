//! Regression — oblique cylinder side-patch seam winding (`tess/cylinder.rs`).
//!
//! An obliquely cut cylinder whose bottom rim is a circle and whose top rim is a
//! tall ellipse used to produce a sliver triangle near the seam whose winding was
//! decided by its own (round-off-dominated) face normal. That single triangle
//! flipped against its neighbours, so a directed mesh edge appeared twice in the
//! same direction and the module's watertight invariant ("each interior edge is
//! shared by exactly two triangles with opposite orientation", `DESIGN.md` §7)
//! was broken even though the mesh stayed manifold.
//!
//! These tests pin the fix: the side patch is now oriented as a whole by an
//! area-weighted outward test, so no sliver can flip independently. We assert the
//! mesh is fully watertight (every directed edge appears exactly once, with its
//! reverse exactly once) and outward-oriented.
//!
//! Every literal carries an `f64` annotation and an explicit tolerance
//! (`DESIGN.md` §12).

use std::collections::HashMap;

use archi_kernel::boolean::{cut, KeepSide};
use archi_kernel::brep::Brep;
use archi_kernel::build::extrude;
use archi_kernel::csg::Profile2d;
use archi_kernel::mass::signed_volume_checked;
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::{Line3, Plane};
use archi_kernel::tess::{tessellate, Mesh, TessOptions};
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::ValidateLevel;

/// Assert the mesh is watertight: every directed edge `(a, b)` appears exactly
/// once and its reverse `(b, a)` exactly once — i.e. each undirected edge is
/// shared by exactly two triangles with opposite orientation (`DESIGN.md` §7).
fn assert_watertight(mesh: &Mesh) {
    assert!(mesh.triangle_count() > 0, "mesh has no triangles");
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

/// Tessellate a single-solid brep's first solid, asserting success.
fn tess(brep: &Brep, o: &TessOptions, tol: &Tol) -> Mesh {
    tessellate(brep, brep.solids[0], o, tol).expect("tessellation")
}

/// The exact reported scenario: r=0.3 / L=3 cylinder, cut by point (0,0,1.5)
/// normal (0,5,1) Below. The (5,1) tilt stretches the elliptical top rim far in
/// z relative to the circular bottom rim, which used to spawn a seam sliver whose
/// winding flipped. The mesh must now be fully watertight and outward-oriented.
#[test]
fn oblique_cut_0_5_1_seam_winding_is_watertight() {
    let tol = Tol::default();
    let r = 0.3_f64;
    let len = 3.0_f64;
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    let cyl = extrude(&Profile2d::circle(r).expect("circle"), &line, len, &tol).expect("cylinder");

    let cut_plane = Plane::new(
        Point3::new(0.0_f64, 0.0_f64, 1.5_f64),
        Vec3::new(0.0_f64, 5.0_f64, 1.0_f64),
    )
    .expect("oblique plane");
    let res = cut(&cyl, cyl.solids[0], &cut_plane, KeepSide::Below, &tol).expect("cut");
    let chopped = res.brep();
    chopped
        .validate(&tol, ValidateLevel::Full)
        .expect("valid oblique cut");

    // Sanity: the closed-form volume matches the report (~0.42412), confirming we
    // are tessellating the same Full-valid solid the bug was found on.
    let analytic = signed_volume_checked(&chopped).expect("closed-form volume");
    assert!(
        (analytic - 0.42412_f64).abs() < 1e-3_f64,
        "unexpected solid: closed-form volume {analytic}"
    );

    // Default display options, exactly as the report specifies.
    let mesh = tess(&chopped, &TessOptions::default(), &tol);

    // The invariant the bug broke: no directed edge is reused, every edge has its
    // opposite-direction twin.
    assert_watertight(&mesh);
    assert!(mesh.signed_volume() > 0.0_f64, "outward orientation");

    // The mesh volume tracks the closed form to chord order (curved patch).
    assert!(
        (mesh.signed_volume() - analytic).abs() < analytic * 2e-2_f64,
        "mesh volume {} vs analytic {analytic}",
        mesh.signed_volume()
    );
}

/// Guard against orientation sensitivity to the seam location and tilt sign: a
/// family of strongly oblique cuts (large in-plane normal, small axial) that all
/// produce a tall elliptical rim must each be watertight and outward-oriented at
/// the default chord tolerance — the regime where seam slivers appear.
#[test]
fn strongly_oblique_cuts_are_watertight() {
    let tol = Tol::default();
    let r = 0.3_f64;
    let len = 3.0_f64;
    let line = Line3::new(Point3::origin(), Vec3::Z).expect("axis");
    let cyl = extrude(&Profile2d::circle(r).expect("circle"), &line, len, &tol).expect("cylinder");

    // (lateral, axial) normal components — steep tilts that elongate the ellipse
    // while the cut plane still meets only the side wall (a clean elliptical rim,
    // not the caps). r=0.3 over a 3 m column gives ample room around mid-height.
    let tilts: [(f64, f64); 4] = [
        (5.0_f64, 1.0_f64),
        (-5.0_f64, 1.0_f64),
        (4.0_f64, 1.0_f64),
        (3.0_f64, 1.0_f64),
    ];
    for (lat, ax) in tilts {
        let cut_plane = Plane::new(
            Point3::new(0.0_f64, 0.0_f64, 1.5_f64),
            Vec3::new(0.0_f64, lat, ax),
        )
        .expect("oblique plane");
        let res = cut(&cyl, cyl.solids[0], &cut_plane, KeepSide::Below, &tol).expect("cut");
        let chopped = res.brep();
        chopped
            .validate(&tol, ValidateLevel::Full)
            .unwrap_or_else(|e| panic!("valid cut for tilt ({lat},{ax}): {e:?}"));

        let mesh = tess(&chopped, &TessOptions::default(), &tol);
        assert_watertight(&mesh);
        assert!(
            mesh.signed_volume() > 0.0_f64,
            "outward orientation for tilt ({lat},{ax})"
        );
    }
}
