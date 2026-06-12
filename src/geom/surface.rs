//! Surface geometry.

use crate::math::Point3;
use crate::primitives::{Cylinder, Plane};

/// Geometric description of the surface a face lies on.
///
/// The surface carries no orientation of its own; the face's
/// [`Sense`](crate::topo::Sense) decides whether the face normal agrees with
/// the surface normal (`DESIGN.md` §3.3, after Fornjot).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum SurfaceGeom {
    /// A planar surface.
    Plane(Plane),
    /// A circular cylinder surface.
    Cylinder(Cylinder),
}

impl SurfaceGeom {
    /// Signed distance from `p` to this surface.
    ///
    /// For a plane this is positive on the normal side; for a cylinder it is
    /// negative inside. Used by the full geometric validation to check that a
    /// loop's vertices lie on the face surface.
    pub fn signed_distance(&self, p: Point3) -> f64 {
        match self {
            SurfaceGeom::Plane(plane) => plane.signed_distance(p),
            SurfaceGeom::Cylinder(cyl) => cyl.signed_distance(p),
        }
    }
}
