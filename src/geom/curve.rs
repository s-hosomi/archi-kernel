//! Curve geometry.

use crate::math::Point3;
use crate::primitives::{Circle3, Ellipse3, Line3};

/// Geometric description of the curve a half-edge runs along.
///
/// A half-edge stores a parameter interval (its `boundary`) on one of these
/// curves; the curve itself is shared with the sibling half-edge. The
/// parameterisation matches the underlying primitive: arc length for
/// [`Line`](CurveGeom::Line), and angle in radians for the conic curves.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum CurveGeom {
    /// A straight line, parameterised by signed arc length from its origin.
    Line(Line3),
    /// A circle, parameterised by angle in radians.
    Circle(Circle3),
    /// An ellipse, parameterised by angle in radians.
    Ellipse(Ellipse3),
}

impl CurveGeom {
    /// Evaluate the curve at parameter `t`.
    ///
    /// `t` is arc length for [`Line`](Self::Line) and an angle in radians for
    /// the conic curves — the same parameterisation a half-edge `boundary`
    /// uses.
    pub fn point_at(&self, t: f64) -> Point3 {
        match self {
            CurveGeom::Line(l) => l.point_at(t),
            CurveGeom::Circle(c) => c.point_at(t),
            CurveGeom::Ellipse(e) => e.point_at(t),
        }
    }

    /// Reverse a boundary interval `[a, b]` into `[b, a]`.
    ///
    /// A sibling half-edge shares the same curve but runs the opposite way, so
    /// its boundary is this curve's boundary reversed. This is a pure
    /// interval-endpoint swap — the parameterisation is shared, so no
    /// per-curve transformation is needed — but it is exposed here so callers
    /// reverse boundaries through the curve abstraction rather than open-coding
    /// the swap.
    #[inline]
    pub fn reverse_param(boundary: [f64; 2]) -> [f64; 2] {
        [boundary[1], boundary[0]]
    }
}
