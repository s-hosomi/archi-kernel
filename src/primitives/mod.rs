//! Analytic geometric primitives.
//!
//! The kernel restricts surfaces to planes and circular cylinders, and
//! curves to lines, circles and ellipses — the closure of plane/cylinder
//! intersections. Everything is exact analytic geometry; there are no
//! NURBS approximations anywhere in the kernel.

mod curve;
mod cylinder;
mod line;
mod plane;

pub use curve::{Circle3, Ellipse3};
pub use cylinder::Cylinder;
pub use line::Line3;
pub use plane::Plane;
