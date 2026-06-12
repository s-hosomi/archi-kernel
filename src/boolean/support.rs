//! Shared helpers for the boolean cut: a 2-D frame on the cutting plane,
//! coordinate quantisation for vertex de-duplication, and the exact-predicate
//! 2-D geometry (signed area, point-in-polygon) used to nest cap loops.
//!
//! Every orientation decision here goes through
//! [`orient2d_exact`](crate::predicate::orient2d_exact) so cap-loop winding and
//! containment stay combinatorially consistent (`synthesis.md` §1); naive cross
//! products are never used for sign decisions.

use crate::math::{Point3, Vec3};
use crate::predicate::orient2d_exact;
use crate::primitives::plane_basis;
use crate::primitives::Plane;
use crate::tolerance::Sign3;

/// Integer-quantised coordinate so identical points de-duplicate despite `f64`
/// not being `Hash`/`Eq`. The scale (1e9) resolves building dimensions
/// (1e-3..1e2 m) far below the 1e-6 m length tolerance without overflowing
/// `i64`; it matches the extrusion builder so a cut interoperates with
/// extruded input.
pub(crate) type CoordKey = (i64, i64, i64);

/// The single quantisation scale (1e9) shared by every coordinate-keying site
/// (the extruder, the cut, and the prismatic builder), so identical points
/// produce identical keys regardless of which subsystem created them. It is the
/// one authoritative definition; all other modules import [`key`] /
/// [`quantize`] from here rather than re-deriving the scale.
pub(crate) const QUANT_SCALE: f64 = 1.0e9_f64;

/// Quantise a single scalar coordinate (or curve parameter) to an integer key
/// at the shared [`QUANT_SCALE`].
#[inline]
pub(crate) fn quantize(x: f64) -> i64 {
    (x * QUANT_SCALE).round() as i64
}

/// Quantise a 3-D point to a [`CoordKey`].
pub(crate) fn key(p: Point3) -> CoordKey {
    (quantize(p.x), quantize(p.y), quantize(p.z))
}

/// An orthonormal 2-D coordinate frame on the cutting plane.
///
/// `u`, `v` are an orthonormal basis of the plane (with `u × v = normal`), and
/// `origin` is the plane's reference point. A 3-D point on the plane projects to
/// `(s, t) = ((p − origin)·u, (p − origin)·v)`. The basis comes from
/// [`plane_basis`], the same seed rule the rest of the kernel uses, so the frame
/// is deterministic.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PlaneFrame {
    origin: Point3,
    u: Vec3,
    v: Vec3,
}

impl PlaneFrame {
    /// Build the frame from a plane.
    pub(crate) fn new(plane: &Plane) -> Self {
        let (u, v) = plane_basis(plane.normal());
        Self {
            origin: plane.point(),
            u,
            v,
        }
    }

    /// Project a 3-D point to plane-local 2-D coordinates.
    #[inline]
    pub(crate) fn project(&self, p: Point3) -> [f64; 2] {
        let d = p - self.origin;
        [d.dot(self.u), d.dot(self.v)]
    }
}

/// Exact signed area of a simple 2-D polygon given by its vertices.
///
/// Positive for counter-clockwise winding, negative for clockwise. The
/// accumulation reuses the same fan decomposition as the volume integral, but
/// each triangle's contribution is summed from `orient2d`'s exact determinant so
/// the winding sign is never lost to round-off.
pub(crate) fn signed_area_2d(poly: &[[f64; 2]]) -> f64 {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let mut acc = 0.0_f64;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        acc += a[0] * b[1] - b[0] * a[1];
    }
    acc / 2.0_f64
}

/// `true` if the 2-D point `p` lies strictly inside the simple polygon `poly`.
///
/// Uses the winding/crossing count with [`orient2d_exact`] to classify each edge
/// crossing, so a point that is collinear with an edge is handled without a
/// round-off flip. Points exactly on the boundary return `false` (we only use
/// this with interior test points, so boundary cases do not arise in practice).
pub(crate) fn point_in_polygon(p: [f64; 2], poly: &[[f64; 2]]) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    // Ray-crossing count along +x, with the exact orientation deciding which
    // side of each edge the point is on.
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let a = poly[i];
        let b = poly[j];
        let crosses = (a[1] > p[1]) != (b[1] > p[1]);
        if crosses {
            // x-coordinate of the edge at height p.y, decided exactly: the point
            // is to the left of edge (a→b) iff orient2d(a, b, p) has the sign of
            // the edge's upward direction.
            let orient = orient2d_exact(a, b, p);
            let upward = b[1] > a[1];
            let left = match orient {
                Sign3::Above => true,
                Sign3::Below => false,
                Sign3::On => false,
            };
            if left == upward {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn signed_area_ccw_positive() {
        let sq = [
            [0.0_f64, 0.0_f64],
            [1.0_f64, 0.0_f64],
            [1.0_f64, 1.0_f64],
            [0.0_f64, 1.0_f64],
        ];
        assert!((signed_area_2d(&sq) - 1.0_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn signed_area_cw_negative() {
        let sq = [
            [0.0_f64, 0.0_f64],
            [0.0_f64, 1.0_f64],
            [1.0_f64, 1.0_f64],
            [1.0_f64, 0.0_f64],
        ];
        assert!(signed_area_2d(&sq) < 0.0_f64);
    }

    #[test]
    fn point_in_polygon_inside_and_outside() {
        let sq = [
            [0.0_f64, 0.0_f64],
            [2.0_f64, 0.0_f64],
            [2.0_f64, 2.0_f64],
            [0.0_f64, 2.0_f64],
        ];
        assert!(point_in_polygon([1.0_f64, 1.0_f64], &sq));
        assert!(!point_in_polygon([3.0_f64, 1.0_f64], &sq));
        assert!(!point_in_polygon([-1.0_f64, 1.0_f64], &sq));
    }

    #[test]
    fn plane_frame_projects_onto_plane() {
        let plane = Plane::new(Point3::origin(), Vec3::Z).expect("plane");
        let frame = PlaneFrame::new(&plane);
        let p = Point3::new(3.0_f64, 4.0_f64, 0.0_f64);
        let xy = frame.project(p);
        assert!((xy[0] * xy[0] + xy[1] * xy[1] - 25.0_f64).abs() < 1e-9_f64);
    }
}
