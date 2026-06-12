//! Opaque geometry identifiers.
//!
//! These newtypes wrap a generational arena handle but expose no geometric
//! value. The topology layer references geometry exclusively through them, so
//! it can name a point / curve / surface without ever holding coordinates.

use crate::geom::{CurveGeom, SurfaceGeom, VertexGeom};
use crate::topo::arena::Id;

/// Opaque handle to vertex geometry in a [`GeomStore`](crate::geom::GeomStore).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PointId(pub(crate) Id<VertexGeom>);

/// Opaque handle to curve geometry in a [`GeomStore`](crate::geom::GeomStore).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CurveId(pub(crate) Id<CurveGeom>);

/// Opaque handle to surface geometry in a [`GeomStore`](crate::geom::GeomStore).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SurfaceId(pub(crate) Id<SurfaceGeom>);
