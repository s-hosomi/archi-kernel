//! Topology layer — the combinatorial B-rep structure.
//!
//! # Invariant: topology carries no geometry
//!
//! This module (and everything under `topo/`) must **never** import the math
//! or primitives layers. Topology entities hold no
//! coordinates; they reference geometry only through the opaque identifiers
//! [`PointId`](crate::geom::PointId), [`CurveId`](crate::geom::CurveId) and
//! [`SurfaceId`](crate::geom::SurfaceId), which carry no geometric value. This
//! is enforced mechanically by a CI grep check, because the type system cannot
//! express it (`DESIGN.md` §3.2; Fornjot #2116 is the cautionary tale).
//!
//! The boundary parameters on [`HalfEdge`] are `f64`, but they are *curve
//! parameters*, not coordinates — they index a point on the shared curve
//! exactly as in Fornjot. Carrying them here does not violate the invariant.
//!
//! # Hierarchy
//!
//! `Solid → Shell → Face → Loop → HalfEdge → Vertex`, a half-edge structure in
//! which the twin of a half-edge is **not** a stored field: siblings (same
//! curve, reversed boundary) are paired up at validation time. There are no
//! upward back-references (`HalfEdge → Loop`, etc.); they would cost integrity
//! maintenance on every face split and are added only when proven necessary
//! (`DESIGN.md` §3.3, `synthesis.md` §3).

pub mod arena;
pub mod store;
pub mod validate;

pub use store::TopoStore;
pub use validate::{Defect, ValidateLevel};

use crate::geom::{CurveId, PointId, SurfaceId};
use arena::Id;

/// A topological vertex. Holds only a reference to its geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Vertex {
    /// The geometry (a point) this vertex stands for.
    pub point: PointId,
}

/// One side of an edge: a directed walk along a curve between two vertices.
///
/// The end vertex is not stored — it is the start vertex of the next half-edge
/// in the loop (Fornjot's redundancy-elimination design). The sibling
/// half-edge shares the same `curve` and has the reversed `boundary`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HalfEdge {
    /// The vertex this half-edge starts at.
    pub start: Id<Vertex>,
    /// The curve this half-edge runs along, shared with its sibling.
    pub curve: CurveId,
    /// The curve-parameter interval `[a, b]` traversed. Direction is encoded by
    /// the ordering of the interval (`a < b` vs `a > b`); the sibling has the
    /// reversed interval. These are curve parameters, not coordinates.
    pub boundary: [f64; 2],
}

/// An ordered, cyclic chain of half-edges bounding a region of a face.
///
/// The successor of `half_edges[i]` is `half_edges[(i + 1) % n]`; there is no
/// stored `next` pointer. A face's outer loop and each of its inner loops are
/// `Loop`s.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Loop {
    /// The half-edges of this loop, in cyclic order.
    pub half_edges: Vec<Id<HalfEdge>>,
}

/// A bounded region of a surface.
///
/// The face carries the orientation (via [`sense`](Self::sense)); the surface
/// itself is unoriented. `outer` is the single outer boundary loop and
/// `inners` are the hole loops (interior loops), which appear when an opening
/// is subtracted.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Face {
    /// The surface this face lies on.
    pub surface: SurfaceId,
    /// Whether the face normal agrees with or opposes the surface normal.
    pub sense: Sense,
    /// The outer boundary loop.
    pub outer: Id<Loop>,
    /// The interior (hole) loops, if any.
    pub inners: Vec<Id<Loop>>,
}

/// An oriented, closed collection of faces.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Shell {
    /// The faces of this shell.
    pub faces: Vec<Id<Face>>,
}

/// A solid: an outer shell plus one shell per internal cavity.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Solid {
    /// The shells of this solid (outer shell first, then cavity shells).
    pub shells: Vec<Id<Shell>>,
}

/// Whether a face's normal agrees with the underlying surface normal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum Sense {
    /// The face normal points the same way as the surface normal.
    Same,
    /// The face normal points opposite to the surface normal.
    Reversed,
}
