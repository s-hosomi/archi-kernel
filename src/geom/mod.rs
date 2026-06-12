//! Geometry store: the analytic geometry that topology entities refer to.
//!
//! The topology layer ([`crate::topo`]) holds no coordinates — it references
//! geometry only through the *opaque* identifiers defined here ([`PointId`],
//! [`CurveId`], [`SurfaceId`]). Those identifiers carry an arena index and
//! generation but no geometric value, so importing them into `topo` does not
//! break the "topology carries no geometry" invariant (`DESIGN.md` §3.2).
//!
//! Geometry itself (points, curves, surfaces) lives in a [`GeomStore`]. Two
//! points worth highlighting:
//!
//! * **Plane canonicalisation** ([`GeomStore::insert_plane`]). The same
//!   conceptual plane derived independently from two members differs by ULPs.
//!   On insertion the store normalises the orientation to a canonical sense and
//!   merges planes that agree within [`Tol`], returning a single shared
//!   [`SurfaceId`]. This is the prerequisite for coplanar detection, sibling
//!   pairing and any future exact path (`DESIGN.md` §2.2).
//! * **Vertex geometry is an open enum** ([`VertexGeom`]). Only `Explicit`
//!   coordinates are implemented today; the symbolic variants reserved for the
//!   future exact predicate path are added later (`synthesis.md` §2-5).

mod curve;
mod ids;
mod plane_store;
mod surface;
mod vertex;

pub use curve::CurveGeom;
pub use ids::{CurveId, PointId, SurfaceId};
pub use surface::SurfaceGeom;
pub use vertex::VertexGeom;

use crate::primitives::Plane;
use crate::tolerance::Tol;
use crate::topo::arena::Arena;

/// Storage for the analytic geometry referenced by topology.
///
/// Points and curves are inserted verbatim. Planes go through
/// [`insert_plane`](Self::insert_plane), which canonicalises and de-duplicates
/// them; cylinders are stored verbatim (no canonicalisation — see
/// `synthesis.md`).
#[derive(Debug, Clone, Default)]
pub struct GeomStore {
    points: Arena<VertexGeom>,
    curves: Arena<CurveGeom>,
    surfaces: Arena<SurfaceGeom>,
}

impl GeomStore {
    /// Create an empty geometry store.
    pub fn new() -> Self {
        Self::default()
    }

    // ── points ───────────────────────────────────────────────────────────

    /// Insert vertex geometry and return its identifier.
    pub fn insert_point(&mut self, geom: VertexGeom) -> PointId {
        PointId(self.points.insert(geom))
    }

    /// Resolve a point identifier.
    pub fn point(&self, id: PointId) -> Option<&VertexGeom> {
        self.points.get(id.0)
    }

    /// Number of stored points.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    // ── curves ───────────────────────────────────────────────────────────

    /// Insert curve geometry and return its identifier.
    ///
    /// Curves are stored verbatim — there is no canonicalisation, because
    /// sibling half-edges already share the same `CurveId` by construction.
    pub fn insert_curve(&mut self, geom: CurveGeom) -> CurveId {
        CurveId(self.curves.insert(geom))
    }

    /// Resolve a curve identifier.
    pub fn curve(&self, id: CurveId) -> Option<&CurveGeom> {
        self.curves.get(id.0)
    }

    /// Number of stored curves.
    pub fn curve_count(&self) -> usize {
        self.curves.len()
    }

    // ── surfaces ─────────────────────────────────────────────────────────

    /// Insert a canonicalised plane, de-duplicating against existing surfaces.
    ///
    /// The returned `flipped` flag is `true` when the supplied plane's normal
    /// pointed opposite to the canonical orientation (so a face built on it
    /// must account for the flip via its [`Sense`](crate::topo::Sense)).
    ///
    /// See [`plane_store`](self) for the canonical-orientation rule and the
    /// tolerance used for de-duplication.
    pub fn insert_plane(&mut self, plane: Plane, tol: &Tol) -> (SurfaceId, bool) {
        plane_store::insert_plane(&mut self.surfaces, plane, tol)
    }

    /// Insert surface geometry verbatim (no canonicalisation).
    ///
    /// Use [`insert_plane`](Self::insert_plane) for planes; this is the escape
    /// hatch for cylinders and for tests that need a raw surface.
    pub fn insert_surface(&mut self, geom: SurfaceGeom) -> SurfaceId {
        SurfaceId(self.surfaces.insert(geom))
    }

    /// Resolve a surface identifier.
    pub fn surface(&self, id: SurfaceId) -> Option<&SurfaceGeom> {
        self.surfaces.get(id.0)
    }

    /// Number of stored surfaces.
    pub fn surface_count(&self) -> usize {
        self.surfaces.len()
    }
}
