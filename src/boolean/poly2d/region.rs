//! Polygonal regions: the input and output domain of the boolean engine.
//!
//! A [`Region`] is a set of [`Contour`]s. Each contour is a closed ring of
//! [`Edge2`] edges. By convention an **outer** boundary runs counter-clockwise
//! (CCW, positive signed area) and a **hole** runs clockwise (CW, negative
//! signed area). A region may contain several connected components and any
//! nesting of holes-in-islands.
//!
//! The engine does *not* require the caller to label which contour is an outer
//! and which is a hole, nor to supply consistent orientation: [`Region::new`]
//! stores contours as-given, and the arrangement determines containment from
//! geometry. Orientation is only *normalized on output* (see
//! [`crate::boolean::poly2d::reconstruct`]). [`Contour::signed_area`] is the tool
//! used for that.

use crate::boolean::poly2d::geom::{Arc, Edge2, Point2};
use crate::tolerance::Tol;

/// A closed ring of edges.
///
/// The ring is implicitly closed: the last edge's `end` connects back to the
/// first edge's `start`. The engine validates and, where helpful, repairs this
/// when ingesting input (consecutive endpoints must coincide within `eps`).
#[derive(Debug, Clone, PartialEq)]
pub struct Contour {
    /// The directed edges of the ring, in traversal order.
    pub edges: Vec<Edge2>,
}

impl Contour {
    /// Construct a contour from a list of directed edges.
    #[inline]
    pub fn new(edges: Vec<Edge2>) -> Self {
        Self { edges }
    }

    /// Construct a contour from an ordered list of vertices, connecting them
    /// with straight segments and closing the ring (last → first).
    ///
    /// Convenience for the common all-segment case. Vertices must not repeat the
    /// closing point (the closure is implicit).
    pub fn from_points(points: &[Point2]) -> Self {
        let n = points.len();
        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            let a = points[i];
            let b = points[(i + 1) % n];
            edges.push(Edge2::seg(a, b));
        }
        Self { edges }
    }

    /// Construct a CCW circular contour from a centre and radius, as two
    /// semicircular arcs over a two-point seam, matching the extruder's seam
    /// convention (`DESIGN.md` §6-1). The result is a closed ring of two
    /// [`Edge2::Arc`] edges traversed counter-clockwise.
    ///
    /// The seam is placed **vertically** (φ = π/2 and φ = 3π/2, i.e. the top and
    /// bottom of the circle) rather than horizontally, so the everyday
    /// horizontally-adjacent tangent-sleeve case (two voids side by side touching
    /// on the x axis) does not put the tangent contact on a seam vertex — which
    /// would pinch the arrangement. A tangency landing exactly on a seam is still
    /// reported as an explicit degeneracy rather than mis-answered.
    pub fn circle(center: Point2, radius: f64) -> Self {
        use std::f64::consts::PI;
        let half = PI;
        Self {
            edges: vec![
                Edge2::Arc(Arc::new(center, radius, PI / 2.0, half)),
                Edge2::Arc(Arc::new(center, radius, 3.0 * PI / 2.0, half)),
            ],
        }
    }

    /// The vertices of the ring (each edge's start point, in order).
    ///
    /// Only meaningful for all-segment contours; arc edges contribute only their
    /// start point here.
    pub fn vertices(&self) -> Vec<Point2> {
        self.edges.iter().map(|e| e.start()).collect()
    }

    /// Signed area of the ring via the shoelace formula over edge endpoints,
    /// **plus a circular-segment correction for every arc edge**.
    ///
    /// Positive for a CCW ring, negative for a CW ring. The shoelace term
    /// integrates the chord polygon; each arc additionally contributes the signed
    /// lens area between the arc and its chord, `½r²(Δθ − sinΔθ)` with `Δθ` the
    /// arc's signed sweep — so a contour mixing straight edges and arcs reports
    /// its exact enclosed area (the same formula `mass::volume` uses).
    pub fn signed_area(&self) -> f64 {
        let mut acc = 0.0_f64;
        for e in &self.edges {
            let a = e.start();
            let b = e.end();
            acc += a.x * b.y - b.x * a.y;
        }
        let mut area = 0.5 * acc;
        for e in &self.edges {
            if let Edge2::Arc(arc) = e {
                let dtheta = arc.sweep;
                area += 0.5 * arc.radius * arc.radius * (dtheta - dtheta.sin());
            }
        }
        area
    }

    /// Reverse the ring's traversal direction (flips orientation sign).
    pub fn reverse(&mut self) {
        self.edges.reverse();
        for e in &mut self.edges {
            *e = match *e {
                Edge2::Seg { start, end } => Edge2::seg(end, start),
                // Arcs are not yet supported; a placeholder reversal keeps the
                // match exhaustive without claiming correctness.
                Edge2::Arc(mut arc) => {
                    arc.start_angle += arc.sweep;
                    arc.sweep = -arc.sweep;
                    Edge2::Arc(arc)
                }
            };
        }
    }

    /// `true` if any edge of this contour is an arc.
    pub fn has_arc(&self) -> bool {
        self.edges.iter().any(|e| e.is_arc())
    }

    /// Number of distinct vertices, merging endpoints within `tol.length`.
    ///
    /// Used to reject degenerate contours (fewer than three distinct vertices
    /// bound no area).
    pub fn distinct_vertex_count(&self, tol: &Tol) -> usize {
        let verts = self.vertices();
        let mut count = 0_usize;
        let n = verts.len();
        for i in 0..n {
            let mut is_new = true;
            for j in 0..i {
                if verts[i].coincident(verts[j], tol) {
                    is_new = false;
                    break;
                }
            }
            if is_new {
                count += 1;
            }
        }
        count
    }
}

/// A polygonal region: any number of contours (outers and holes, possibly
/// several connected components and nested holes-in-islands).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Region {
    /// The contours making up the region.
    pub contours: Vec<Contour>,
}

impl Region {
    /// An empty region (the empty set).
    #[inline]
    pub fn empty() -> Self {
        Self {
            contours: Vec::new(),
        }
    }

    /// Construct a region from its contours.
    #[inline]
    pub fn new(contours: Vec<Contour>) -> Self {
        Self { contours }
    }

    /// `true` if the region has no contours.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.contours.is_empty()
    }

    /// Construct a single-contour region from an ordered vertex list (segments).
    pub fn from_points(points: &[Point2]) -> Self {
        Self {
            contours: vec![Contour::from_points(points)],
        }
    }

    /// Construct a single-contour circular region (CCW), as two semicircle arcs.
    pub fn circle(center: Point2, radius: f64) -> Self {
        Self {
            contours: vec![Contour::circle(center, radius)],
        }
    }

    /// Total signed area: sum of all contour signed areas.
    ///
    /// For a well-formed region (CCW outers, CW holes) this equals the true
    /// enclosed area (outer areas minus hole areas).
    pub fn signed_area(&self) -> f64 {
        self.contours.iter().map(|c| c.signed_area()).sum()
    }

    /// Absolute enclosed area.
    pub fn area(&self) -> f64 {
        self.signed_area().abs()
    }

    /// `true` if any contour of this region contains an arc edge.
    pub fn has_arc(&self) -> bool {
        self.contours.iter().any(|c| c.has_arc())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_square_ccw() -> Contour {
        Contour::from_points(&[
            Point2::new(0.0_f64, 0.0_f64),
            Point2::new(1.0_f64, 0.0_f64),
            Point2::new(1.0_f64, 1.0_f64),
            Point2::new(0.0_f64, 1.0_f64),
        ])
    }

    #[test]
    fn signed_area_ccw_positive() {
        let c = unit_square_ccw();
        assert!((c.signed_area() - 1.0_f64).abs() <= 1e-12_f64);
    }

    #[test]
    fn reverse_flips_sign() {
        let mut c = unit_square_ccw();
        c.reverse();
        assert!((c.signed_area() + 1.0_f64).abs() <= 1e-12_f64);
    }

    #[test]
    fn distinct_vertices_merges_near_duplicates() {
        let tol = Tol::default();
        let c = Contour::from_points(&[
            Point2::new(0.0_f64, 0.0_f64),
            Point2::new(1.0_f64, 0.0_f64),
            Point2::new(1.0_f64, 1e-9_f64), // ~coincident with previous
            Point2::new(0.0_f64, 1.0_f64),
        ]);
        assert_eq!(c.distinct_vertex_count(&tol), 3);
    }
}
