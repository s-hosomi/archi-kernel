//! Closed-form intersections between analytic primitives.
//!
//! Because the kernel only admits planes and circular cylinders, every
//! surface-surface intersection here is solved in closed form — there is no
//! numerical curve marching. Cylinder × cylinder (a degree-4 space curve) is
//! intentionally out of scope for v0; round members piercing each other are
//! rare in building structures and will be added later if needed.

mod line_plane;
mod plane_cylinder;
mod plane_plane;

pub use line_plane::line_plane;
pub use plane_cylinder::{plane_cylinder, PlaneCylinder};
pub use plane_plane::{plane_plane, PlanePlane};
