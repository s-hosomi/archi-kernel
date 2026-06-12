//! Geometric predicates — the single funnel for sign decisions.
//!
//! Every geometric predicate in the kernel is routed through this module
//! rather than being open-coded inside boolean or validation logic
//! (`DESIGN.md` §3.5). The funnel deliberately distinguishes **two layers of
//! meaning**:
//!
//! 1. **Tolerant 3-value classification** ([`Sign3`]). Answers "is this *on* the
//!    plane, as a matter of design intent?" The `On` band is a first-class
//!    answer, because coincident geometry (a slab face on a wall face) is the
//!    common case in buildings. This is the layer implemented now.
//! 2. **Exact sign** (Phase 3a). Guards against true round-off mistakes
//!    *outside* the `On` band, by introducing Shewchuk's `orient3d` via the
//!    `robust` crate. The non-transitivity of `ε` comparison is inherent to the
//!    tolerant layer and does not disappear, but the exact layer keeps the
//!    ground under the classification from shifting. Not yet implemented.
//!
//! The argument type is [`VertexGeom`] (not a bare `[f64; 3]`) so that the
//! signature is already compatible with the future symbolic / implicit-point
//! representations; only `Explicit` is handled today (`synthesis.md` §2-5).

use crate::geom::{GeomStore, VertexGeom};
use crate::primitives::Plane;
use crate::tolerance::{Sign3, Tol};

/// Tolerant classification of which side of `plane` the point `p` lies on.
///
/// Returns [`Sign3::Above`] on the plane's normal side, [`Sign3::Below`] on the
/// far side, and [`Sign3::On`] within [`Tol::length`] of the plane.
///
/// `geom` is accepted so that symbolic vertex representations can be resolved
/// without the caller pre-evaluating coordinates; today every [`VertexGeom`] is
/// explicit, so it is currently unused but kept in the signature for the exact
/// path.
pub fn side_of_plane(plane: &Plane, p: &VertexGeom, geom: &GeomStore, tol: &Tol) -> Sign3 {
    let _ = geom; // reserved for symbolic point resolution (Phase 3a / Phase 8)
    match p {
        VertexGeom::Explicit(point) => tol.classify_length(plane.signed_distance(*point)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Point3, Vec3};

    #[test]
    fn classifies_three_ways() {
        let geom = GeomStore::new();
        let tol = Tol::default();
        let plane = Plane::new(Point3::origin(), Vec3::Z).expect("valid plane");
        let above = VertexGeom::Explicit(Point3::new(0.0_f64, 0.0_f64, 1.0_f64));
        let below = VertexGeom::Explicit(Point3::new(0.0_f64, 0.0_f64, -1.0_f64));
        let on = VertexGeom::Explicit(Point3::new(5.0_f64, 5.0_f64, 0.0_f64));
        assert_eq!(side_of_plane(&plane, &above, &geom, &tol), Sign3::Above);
        assert_eq!(side_of_plane(&plane, &below, &geom, &tol), Sign3::Below);
        assert_eq!(side_of_plane(&plane, &on, &geom, &tol), Sign3::On);
    }
}
