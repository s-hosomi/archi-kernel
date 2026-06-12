//! Coplanar degeneracy suite for the prismatic boolean (`DESIGN.md` §4.3).
//!
//! The first two tests are the two reference counter-examples the design fixes
//! the coplanar residency rule against (a difference `A − B` where a face of `A`
//! lies on `∂B`):
//!
//! 1. **Normals agree → drop.** A through difference whose shared bottom face at
//!    `z = 0` must *not* appear in the result.
//! 2. **Normals oppose → keep.** A pure-contact difference that removes nothing,
//!    so `A`'s bottom face survives and the result is `A` whole.
//!
//! The remaining tests pin the rest of the truth table (identical-solid
//! difference, shared-face union, contact-only intersection, overlap
//! intersection). Every literal carries an `f64` annotation and an explicit
//! tolerance, per `DESIGN.md` §12.

use archi_kernel::boolean::prismatic::{self, ExtrudeLeaf};
use archi_kernel::csg::Profile2d;
use archi_kernel::geom::SurfaceGeom;
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::ValidateLevel;
use archi_kernel::Brep;

/// Volume tolerance: well above accumulated `f64` round-off at building scale,
/// far below any real error.
const VOL_EPS: f64 = 1e-9;

/// A box `[x0, x0+wx] × [y0, y0+wy] × [z0, z0+wz]` as a rectangle extruded along
/// `+z`. Members are placed in world coordinates (`DESIGN.md` §5.1).
///
/// For a `+z` extrusion the profile's local axes are `(u, v) = (Y, −X)`
/// ([`plane_basis`](crate::primitives) seed rule), so the profile half-width
/// runs along `Y` and the half-height along `X`. The arguments are nonetheless
/// world `wx, wy` and are mapped accordingly.
fn box_leaf(x0: f64, y0: f64, z0: f64, wx: f64, wy: f64, wz: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(wy / 2.0, wx / 2.0).expect("valid rect"),
        origin: Point3::new(x0 + wx / 2.0, y0 + wy / 2.0, z0),
        axis: Vec3::Z,
        length: wz,
    }
}

/// The unit cube `[0,1]³`.
fn unit_cube() -> ExtrudeLeaf {
    box_leaf(0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64, 1.0_f64)
}

/// `true` if the B-rep has a planar face whose plane passes through `z = z0`
/// (within tolerance) and is horizontal (normal ≈ ±z).
fn has_horizontal_face_at_z(brep: &Brep, z0: f64, tol: &Tol) -> bool {
    for solid_id in &brep.solids {
        let Some(solid) = brep.topo.solids.get(*solid_id) else {
            continue;
        };
        for shell_id in &solid.shells {
            let Some(shell) = brep.topo.shells.get(*shell_id) else {
                continue;
            };
            for face_id in &shell.faces {
                let Some(face) = brep.topo.faces.get(*face_id) else {
                    continue;
                };
                if let Some(SurfaceGeom::Plane(plane)) = brep.geom.surface(face.surface) {
                    let n = plane.normal().as_vec();
                    let horizontal = n.x.abs() <= 1e-9_f64 && n.y.abs() <= 1e-9_f64;
                    let at_z = plane
                        .signed_distance(Point3::new(0.0_f64, 0.0_f64, z0))
                        .abs()
                        <= tol.length;
                    if horizontal && at_z {
                        return true;
                    }
                }
            }
        }
    }
    false
}

#[test]
fn coplanar_normals_agree_face_dropped() {
    // A = [0,1]³, B = [0,1]² × [0,0.5] (a through difference from the bottom).
    // Their bottom faces at z = 0 coincide with the *same* outward normal (−z),
    // so the correct result [0,1]² × [0.5,1] has NO face at z = 0.
    let tol = Tol::default();
    let a = unit_cube();
    let b = box_leaf(0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64, 0.5_f64);

    let result = prismatic::difference(&a, &b, &tol).expect("difference");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight result");

    assert!(
        (result.signed_volume() - 0.5_f64).abs() <= VOL_EPS,
        "volume = {}",
        result.signed_volume()
    );
    assert!(
        !has_horizontal_face_at_z(&result, 0.0_f64, &tol),
        "the shared bottom face at z = 0 must be dropped (normals agree)"
    );
    // The new top of the remaining material is at z = 0.5.
    assert!(has_horizontal_face_at_z(&result, 0.5_f64, &tol));
}

#[test]
fn coplanar_normals_oppose_face_kept() {
    // A = [0,1]³, B = [0,1]² × [−0.5,0] (B sits below A, sharing only the z = 0
    // plane in pure contact). A − B = A: nothing is removed, and A's bottom face
    // at z = 0 must remain.
    let tol = Tol::default();
    let a = unit_cube();
    let b = box_leaf(0.0_f64, 0.0_f64, -0.5_f64, 1.0_f64, 1.0_f64, 0.5_f64);

    let result = prismatic::difference(&a, &b, &tol).expect("difference");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight result");

    assert!(
        (result.signed_volume() - 1.0_f64).abs() <= VOL_EPS,
        "volume = {}",
        result.signed_volume()
    );
    assert!(
        has_horizontal_face_at_z(&result, 0.0_f64, &tol),
        "A's bottom face at z = 0 must survive (normals oppose)"
    );
}

#[test]
fn identical_solid_difference_is_empty() {
    let tol = Tol::default();
    let a = unit_cube();
    let result = prismatic::difference(&a, &a, &tol).expect("difference");
    assert!(result.solids.is_empty(), "A − A must be empty");
    assert!(result.signed_volume().abs() <= VOL_EPS);
}

#[test]
fn shared_face_union_merges_to_one_box() {
    // Two unit-square columns stacked with a shared face at z = 0.5; their union
    // is the single box [0,1]² × [0,1] with the internal face gone.
    let tol = Tol::default();
    let lower = box_leaf(0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64, 0.5_f64);
    let upper = box_leaf(0.0_f64, 0.0_f64, 0.5_f64, 1.0_f64, 1.0_f64, 0.5_f64);

    let result = prismatic::union_pair(&lower, &upper, &tol).expect("union");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight result");

    assert_eq!(result.solids.len(), 1usize, "one merged solid");
    assert!(
        (result.signed_volume() - 1.0_f64).abs() <= VOL_EPS,
        "volume = {}",
        result.signed_volume()
    );
    assert!(
        !has_horizontal_face_at_z(&result, 0.5_f64, &tol),
        "the internal shared face at z = 0.5 must vanish"
    );
}

#[test]
fn contact_only_intersection_is_empty() {
    // Two boxes touching only on the z = 0.5 plane share zero volume.
    let tol = Tol::default();
    let lower = box_leaf(0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64, 0.5_f64);
    let upper = box_leaf(0.0_f64, 0.0_f64, 0.5_f64, 1.0_f64, 1.0_f64, 0.5_f64);

    let result = prismatic::intersection(&lower, &upper, &tol).expect("intersection");
    assert!(result.solids.is_empty(), "contact-only ∩ is empty");
    assert!(result.signed_volume().abs() <= VOL_EPS);
}

#[test]
fn overlap_intersection_is_the_shared_box() {
    // Two unit cubes offset by (0.5, 0.5, 0.5): their intersection is the
    // [0.5,1]³ box of volume 0.125.
    let tol = Tol::default();
    let a = unit_cube();
    let b = box_leaf(0.5_f64, 0.5_f64, 0.5_f64, 1.0_f64, 1.0_f64, 1.0_f64);

    let result = prismatic::intersection(&a, &b, &tol).expect("intersection");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight result");

    assert_eq!(result.solids.len(), 1usize);
    assert!(
        (result.signed_volume() - 0.125_f64).abs() <= VOL_EPS,
        "volume = {}",
        result.signed_volume()
    );
}
