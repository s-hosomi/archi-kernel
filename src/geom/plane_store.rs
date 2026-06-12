//! Plane canonicalisation and de-duplication.
//!
//! Planes derived independently from different members differ by ULPs even
//! when they are "the same" plane by design intent. Before storing one we:
//!
//! 1. **Canonicalise the orientation.** The normal is flipped so that its first
//!    component that is non-zero (under [`Tol::angular`]) is positive. The
//!    returned `flipped` flag records whether this reversed the caller's
//!    normal.
//! 2. **De-duplicate.** The store is scanned for an existing plane whose normal
//!    agrees within [`Tol::angular`] and whose signed distance from the origin
//!    agrees within [`Tol::length`]. If one is found its identifier is
//!    returned; otherwise the canonicalised plane is inserted.
//!
//! The scan is linear today. The number of distinct planes within a single
//! member is only a few dozen, so this is adequate.
//
// TODO(perf): replace the linear scan with a spatial hash keyed on the
// quantised (normal, origin-distance) tuple once member sizes grow or
// inter-member stores appear (Phase 3+).

use crate::geom::ids::SurfaceId;
use crate::geom::surface::SurfaceGeom;
use crate::math::{Unit3, Vec3};
use crate::primitives::Plane;
use crate::tolerance::Tol;
use crate::topo::arena::Arena;

/// Canonicalise `plane`'s orientation and de-duplicate it against `surfaces`.
///
/// Returns the shared [`SurfaceId`] and whether the caller's normal was
/// flipped to reach the canonical orientation.
pub(crate) fn insert_plane(
    surfaces: &mut Arena<SurfaceGeom>,
    plane: Plane,
    tol: &Tol,
) -> (SurfaceId, bool) {
    let (canonical, flipped) = canonicalise(plane, tol);

    // Linear search for an existing coincident plane.
    let canon_normal = canonical.normal();
    let canon_dist = canonical.signed_distance(crate::math::Point3::origin());
    for (id, surf) in surfaces.iter() {
        if let SurfaceGeom::Plane(existing) = surf {
            if normals_agree(existing.normal(), canon_normal, tol)
                && tol.eq_length(
                    existing.signed_distance(crate::math::Point3::origin()),
                    canon_dist,
                )
            {
                return (SurfaceId(id), flipped);
            }
        }
    }

    let id = surfaces.insert(SurfaceGeom::Plane(canonical));
    (SurfaceId(id), flipped)
}

/// Flip `plane` so its normal points to the canonical half-space, returning the
/// canonicalised plane and whether the input normal was reversed.
///
/// The canonical rule: the first component (x, then y, then z) that is non-zero
/// under [`Tol::angular`] must be positive.
fn canonicalise(plane: Plane, tol: &Tol) -> (Plane, bool) {
    let n = plane.normal().as_vec();
    let flip = should_flip(n, tol);
    if flip {
        let flipped_normal = Unit3::new_unchecked(-n);
        (Plane::new_unchecked(plane.point(), flipped_normal), true)
    } else {
        (plane, false)
    }
}

/// Decide whether a normal points to the negative canonical half-space.
///
/// The first component (x, then y, then z) whose absolute value is **strictly
/// greater** than [`Tol::angular`] determines the canonical orientation: if
/// that component is negative the normal is flipped.  Components whose absolute
/// value is at most `tol.angular` are treated as "zero-like" and skipped.
///
/// Keeping the threshold strict (`>`) rather than inclusive (`>=`) is
/// deliberate: a component at exactly `tol.angular` is ambiguous (either sign
/// may be noise), so skipping it lets the next, more dominant component decide.
/// The consequence is that two normals whose first non-skipped component have
/// opposite signs at exactly `tol.angular` end up with the same canonical
/// orientation (both skip to the dominant component).  Any residual difference
/// of up to `2·tol.angular` between canonical normals is absorbed by
/// [`normals_agree`], which uses `2·tol.angular` as its threshold for exactly
/// this reason.
fn should_flip(n: Vec3, tol: &Tol) -> bool {
    for component in [n.x, n.y, n.z] {
        if component.abs() > tol.angular {
            return component < 0.0;
        }
    }
    // All components within angular tolerance of zero: degenerate normal.
    // Plane construction already guarantees a unit normal, so this is
    // unreachable in practice; leave the orientation unchanged.
    false
}

/// `true` if two canonicalised unit normals agree within `2·angular` tolerance.
///
/// Both inputs are already canonicalised to the same half-space by
/// [`should_flip`].  Because `should_flip` treats components up to
/// `tol.angular` as "zero-like" and skips them, two normals that are truly the
/// same plane (within `tol.angular`) may end up with canonical forms whose
/// sub-dominant component differs by up to `2·tol.angular` (the residual from
/// both sides of the skip boundary).  Using `2·tol.angular` as the comparison
/// threshold ensures those cases are correctly identified as duplicates.
///
/// The chord-length bound `|a − b|` is used rather than the exact angle; for
/// the small differences in question the two are numerically equivalent.
fn normals_agree(a: Unit3, b: Unit3, tol: &Tol) -> bool {
    let d = a.as_vec() - b.as_vec();
    d.norm() <= 2.0 * tol.angular
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point3;

    fn plane(px: f64, py: f64, pz: f64, nx: f64, ny: f64, nz: f64) -> Plane {
        Plane::new(Point3::new(px, py, pz), Vec3::new(nx, ny, nz)).expect("valid plane")
    }

    #[test]
    fn flip_makes_first_component_positive() {
        let tol = Tol::default();
        let (canon, flipped) = canonicalise(
            plane(0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, -1.0_f64),
            &tol,
        );
        assert!(flipped);
        assert!(canon.normal().z > 0.0_f64);
    }

    #[test]
    fn no_flip_when_already_canonical() {
        let tol = Tol::default();
        let (_canon, flipped) = canonicalise(
            plane(0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64),
            &tol,
        );
        assert!(!flipped);
    }
}
