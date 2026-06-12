//! archi-kernel — domain-specific B-rep geometry kernel for building simulation.
//!
//! Scope is deliberately restricted to the surfaces that occur in building
//! structures — planes and circular cylinders — so that every surface-surface
//! intersection has a closed-form solution. This removes the numerical
//! marching that makes general-purpose B-rep kernels fragile.
//!
//! All lengths are in metres (SI), all angles in radians. See `DESIGN.md`
//! for the architecture roadmap and verification strategy.

pub mod error;
pub mod intersect;
pub mod math;
pub mod primitives;
pub mod tolerance;

pub use error::KernelError;
pub use math::{Point3, Unit3, Vec3};
pub use primitives::{Circle3, Cylinder, Ellipse3, Line3, Plane};
pub use tolerance::{Sign3, Tol};
