//! archi-kernel — domain-specific B-rep geometry kernel for building simulation.
//!
//! Scope is deliberately restricted to the surfaces that occur in building
//! structures — planes and circular cylinders — so that every surface-surface
//! intersection has a closed-form solution. This removes the numerical
//! marching that makes general-purpose B-rep kernels fragile.
//!
//! All lengths are in metres (SI), all angles in radians. See `DESIGN.md`
//! for the architecture roadmap and verification strategy.

pub mod brep;
pub mod build;
pub mod csg;
pub mod error;
pub mod geom;
pub mod intersect;
pub mod mass;
pub mod math;
pub mod predicate;
pub mod primitives;
pub mod profile;
pub mod tolerance;
pub mod topo;

pub use brep::Brep;
pub use build::extrude;
pub use error::KernelError;
pub use geom::{CurveGeom, CurveId, GeomStore, PointId, SurfaceGeom, SurfaceId, VertexGeom};
pub use math::{Point3, Unit3, Vec3};
pub use primitives::{Circle3, Cylinder, Ellipse3, Line3, Plane};
pub use profile::ProfileGeom;
pub use tolerance::{Sign3, Tol};
pub use topo::{
    Defect, Face, HalfEdge, Loop, Sense, Shell, Solid, TopoStore, ValidateLevel, Vertex,
};
