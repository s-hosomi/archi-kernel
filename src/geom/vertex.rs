//! Vertex geometry.

use crate::math::Point3;

/// Geometric description of a topological vertex.
///
/// Only [`Explicit`](VertexGeom::Explicit) is implemented today. The enum is
/// `#[non_exhaustive]` so that the symbolic representations reserved for the
/// future exact predicate path — a line/plane intersection or a three-plane
/// intersection carried without ever evaluating coordinates — can be added in
/// a semver-compatible way (`DESIGN.md` §3.4, `synthesis.md` §2-5).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum VertexGeom {
    /// An explicit design coordinate (e.g. a profile corner).
    Explicit(Point3),
}

impl VertexGeom {
    /// The explicit coordinate of this vertex.
    ///
    /// Returns `None` for symbolic variants that have not been evaluated to a
    /// point. Today every variant is explicit, so this always returns `Some`.
    pub fn as_point(&self) -> Option<Point3> {
        match self {
            VertexGeom::Explicit(p) => Some(*p),
        }
    }
}
