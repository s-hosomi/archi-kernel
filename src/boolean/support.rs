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

/// Signed area of a simple 2-D polygon given by its vertices (shoelace formula).
///
/// Positive for counter-clockwise winding, negative for clockwise. This is a
/// naive floating-point shoelace implementation — **not** an exact-arithmetic
/// computation. Round-off can affect the result for nearly-degenerate polygons,
/// but the sign is reliable for well-conditioned inputs encountered in practice
/// (building-scale geometry in SI metres).
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

/// Union-find: find with path halving.
#[inline]
pub(crate) fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// Union-find: union by root.
#[inline]
pub(crate) fn uf_union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (uf_find(parent, a), uf_find(parent, b));
    if ra != rb {
        parent[ra] = rb;
    }
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

/// `true` if the 2-D point `p` lies strictly inside the closed loop given by its
/// ordered vertex ring `ring`. A re-export of [`point_in_polygon`] under a name
/// the section module uses for loops; arc bulges are ignored (the ring's chord
/// polygon is tested), which is exact for the containment of an interior
/// representative point well away from the boundary.
pub(crate) fn point_in_loop_2d(p: [f64; 2], ring: &[[f64; 2]]) -> bool {
    point_in_polygon(p, ring)
}

/// Signed area of a loop in a 2-D frame: the shoelace area of the vertex `ring`
/// **plus a circular-segment correction for every arc edge**.
///
/// This mirrors [`crate::boolean::poly2d`] `Contour::signed_area`: the shoelace
/// term integrates the chord polygon, and each arc adds the signed lens area
/// `½r²(Δθ − sinΔθ)` between the arc and its chord. `arcs` carries one
/// `(radius, signed_sweep)` per arc edge, the sweep already in `(−π, π]` (a
/// section arc spans at most a semicircle). So a loop mixing straight edges and
/// arcs reports its exact enclosed area, and its sign is the loop's winding.
pub(crate) fn loop_signed_area_2d(ring: &[[f64; 2]], arcs: &[(f64, f64)]) -> f64 {
    let mut area = signed_area_2d(ring);
    for &(radius, dtheta) in arcs {
        area += 0.5 * radius * radius * (dtheta - dtheta.sin());
    }
    area
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
