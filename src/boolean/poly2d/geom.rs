//! 2-D geometric primitives: points, vectors, the [`Edge2`] curve enum, and the
//! robust [`orient2d`] predicate.
//!
//! # Edge model — designed for arcs from day one
//!
//! [`Edge2`] is an `enum { Seg, Arc }` from the start, even though only [`Seg`]
//! is implemented. This is a **design acceptance condition**: adding arc support
//! must be *adding code*, never *rewriting* it. The structural seams that make
//! this true are:
//!
//! * **Intersection** dispatches on the edge-kind pair (`seg×seg`, `seg×arc`,
//!   `arc×arc`) in [`crate::boolean::poly2d::arrangement`]. Today only `seg×seg` is
//!   wired; the other arms return [`crate::boolean::poly2d::Poly2Error::ArcNotYetSupported`].
//! * **Point-on-edge classification** and **edge sampling** (midpoint, tangent)
//!   are methods on [`Edge2`] that already `match` on the variant, so arc
//!   geometry plugs in without touching call sites.
//! * **Face classification** uses a representative interior point + a ray cast
//!   that asks each edge "how many times does a ray cross you?" — again an
//!   [`Edge2`] method, arc-ready by construction.
//!
//! [`Seg`]: Edge2::Seg

use crate::tolerance::Tol;

/// Squared length tolerance derived from a `Tol`.
///
/// Helper used inside this module because the parent `Tol` only exposes
/// `length`, not a pre-computed `eps_sq`. Computed inline to avoid storing
/// redundant state.
#[inline]
pub(super) fn eps_sq(tol: &Tol) -> f64 {
    tol.length * tol.length
}

/// A point in the 2-D plane (coordinates in metres, SI).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    /// x coordinate.
    pub x: f64,
    /// y coordinate.
    pub y: f64,
}

impl Point2 {
    /// Construct a point from its two coordinates.
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Vector from `self` to `other` (`other - self`).
    #[inline]
    pub fn to(self, other: Point2) -> Vec2 {
        Vec2::new(other.x - self.x, other.y - self.y)
    }

    /// Squared Euclidean distance to `other`.
    #[inline]
    pub fn dist_sq(self, other: Point2) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    /// Euclidean distance to `other`.
    #[inline]
    pub fn dist(self, other: Point2) -> f64 {
        self.dist_sq(other).sqrt()
    }

    /// `true` if `self` and `other` are within `tol.length` of each other.
    #[inline]
    pub fn coincident(self, other: Point2, tol: &Tol) -> bool {
        self.dist_sq(other) <= eps_sq(tol)
    }

    /// Linear interpolation: `self` at `t = 0`, `other` at `t = 1`.
    #[inline]
    pub fn lerp(self, other: Point2, t: f64) -> Point2 {
        Point2::new(
            self.x + (other.x - self.x) * t,
            self.y + (other.y - self.y) * t,
        )
    }
}

/// A free vector in the 2-D plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec2 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
}

impl Vec2 {
    /// Construct a vector from its two components.
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, other: Vec2) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// 2-D cross product (z component of the 3-D cross), aka the perp-dot.
    #[inline]
    pub fn cross(self, other: Vec2) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Squared length.
    #[inline]
    pub fn len_sq(self) -> f64 {
        self.x * self.x + self.y * self.y
    }

    /// Length.
    #[inline]
    pub fn len(self) -> f64 {
        self.len_sq().sqrt()
    }
}

/// Winding / turn orientation of an ordered point triple, classified with the
/// robust adaptive-precision predicate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orient {
    /// `a → b → c` turns left (counter-clockwise); signed area positive.
    Left,
    /// The three points are collinear (to full floating-point precision).
    Collinear,
    /// `a → b → c` turns right (clockwise); signed area negative.
    Right,
}

/// Robust orientation of the triple `(a, b, c)` using Shewchuk's exact
/// adaptive-precision `orient2d` (via `crate::predicate::orient2d_exact`).
///
/// This is the **foundation of the arrangement's combinatorial robustness**:
/// every left/right/collinear decision that determines topology goes through
/// here, so the decisions are sign-exact and mutually consistent even when the
/// floating-point determinant would round the wrong way. Tolerant (`eps`-based)
/// classification is layered *on top* for "is this a deliberate coincidence"
/// questions; raw geometric sign questions use this.
///
/// Routed through `crate::predicate::orient2d_exact` rather than calling
/// `robust` directly, in compliance with the kernel's isolation rule
/// (`DESIGN.md` §3.5): the `robust` crate is confined to `src/predicate/`.
#[inline]
pub fn orient2d(a: Point2, b: Point2, c: Point2) -> Orient {
    use crate::predicate::orient2d_exact;
    use crate::tolerance::Sign3;
    match orient2d_exact([a.x, a.y], [b.x, b.y], [c.x, c.y]) {
        Sign3::Above => Orient::Left,
        Sign3::Below => Orient::Right,
        Sign3::On => Orient::Collinear,
    }
}

/// A circular arc: centre, radius, signed angular sweep, and orientation.
///
/// **Not yet implemented** — present so the [`Edge2`] enum and every algorithm
/// that dispatches on edge kind is arc-ready. Angles are in radians; the arc
/// runs from `start_angle` to `start_angle + sweep` (so `sweep`'s sign encodes
/// CCW / CW traversal). The fields are kept even though unused to lock the
/// representation; arc geometry is filled in by *adding* method arms.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Arc {
    /// Centre of the circle the arc lies on.
    pub center: Point2,
    /// Radius (positive).
    pub radius: f64,
    /// Angle of the arc's start point, radians, measured from +x.
    pub start_angle: f64,
    /// Signed angular sweep, radians. Positive = CCW, negative = CW.
    pub sweep: f64,
}

impl Arc {
    /// Construct an arc.
    #[inline]
    pub fn new(center: Point2, radius: f64, start_angle: f64, sweep: f64) -> Self {
        Self {
            center,
            radius,
            start_angle,
            sweep,
        }
    }

    /// The arc's start point.
    #[inline]
    pub fn start(&self) -> Point2 {
        self.point_at_angle(self.start_angle)
    }

    /// The arc's end point.
    #[inline]
    pub fn end(&self) -> Point2 {
        self.point_at_angle(self.start_angle + self.sweep)
    }

    #[inline]
    fn point_at_angle(&self, ang: f64) -> Point2 {
        Point2::new(
            self.center.x + self.radius * ang.cos(),
            self.center.y + self.radius * ang.sin(),
        )
    }
}

/// A directed boundary edge: a line segment or a circular arc.
///
/// The direction matters: `start → end` is the traversal direction, and the
/// region's interior is conventionally on the left of a CCW outer loop.
///
/// Only [`Edge2::Seg`] is implemented. Any operation that encounters an
/// [`Edge2::Arc`] returns [`crate::boolean::poly2d::Poly2Error::ArcNotYetSupported`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Edge2 {
    /// A straight segment from `start` to `end`.
    Seg {
        /// Segment start point.
        start: Point2,
        /// Segment end point.
        end: Point2,
    },
    /// A circular arc (not yet implemented; see [`Arc`]).
    Arc(Arc),
}

impl Edge2 {
    /// Construct a segment edge.
    #[inline]
    pub fn seg(start: Point2, end: Point2) -> Self {
        Edge2::Seg { start, end }
    }

    /// The edge's start point (traversal origin).
    #[inline]
    pub fn start(&self) -> Point2 {
        match self {
            Edge2::Seg { start, .. } => *start,
            Edge2::Arc(a) => a.start(),
        }
    }

    /// The edge's end point (traversal terminus).
    #[inline]
    pub fn end(&self) -> Point2 {
        match self {
            Edge2::Seg { end, .. } => *end,
            Edge2::Arc(a) => a.end(),
        }
    }

    /// `true` if this edge is an arc (used by callers to fail fast).
    #[inline]
    pub fn is_arc(&self) -> bool {
        matches!(self, Edge2::Arc(_))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orient_left_right_collinear() {
        let a = Point2::new(0.0_f64, 0.0_f64);
        let b = Point2::new(1.0_f64, 0.0_f64);
        assert_eq!(orient2d(a, b, Point2::new(0.0_f64, 1.0_f64)), Orient::Left);
        assert_eq!(
            orient2d(a, b, Point2::new(0.0_f64, -1.0_f64)),
            Orient::Right
        );
        assert_eq!(
            orient2d(a, b, Point2::new(2.0_f64, 0.0_f64)),
            Orient::Collinear
        );
    }

    #[test]
    fn cross_and_dot() {
        let u = Vec2::new(1.0_f64, 0.0_f64);
        let v = Vec2::new(0.0_f64, 1.0_f64);
        assert!((u.cross(v) - 1.0_f64).abs() <= 1e-15_f64);
        assert!(u.dot(v).abs() <= 1e-15_f64);
    }

    #[test]
    fn edge_endpoints() {
        let e = Edge2::seg(Point2::new(0.0_f64, 0.0_f64), Point2::new(3.0_f64, 4.0_f64));
        assert_eq!(e.start(), Point2::new(0.0_f64, 0.0_f64));
        assert_eq!(e.end(), Point2::new(3.0_f64, 4.0_f64));
        assert!(!e.is_arc());
    }
}
