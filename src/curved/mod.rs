//! Curved panel primitives with UV-space trim loops.
//!
//! This module is intentionally separate from the main B-rep/CSG evaluator. It
//! starts the curved-opening path as a trimmed surface primitive: holes are
//! supplied as loops in the surface parameter domain and tessellated for display.
//! General surface-surface boolean evaluation remains out of scope.

mod cylinder_panel;
mod domain;
mod error;
mod mesh;
mod sphere_panel;
mod trim;

pub use cylinder_panel::{
    tessellate_cylinder_panel, tessellate_thick_cylinder_panel, CylinderPanel,
    CylinderPanelOptions, ThickCylinderPanel,
};
pub use error::CurvedError;
pub use mesh::SurfaceMesh;
pub use sphere_panel::{
    tessellate_sphere_panel, tessellate_thick_sphere_panel, SpherePanel, SpherePanelOptions,
    ThickSpherePanel,
};
pub use trim::{TrimEdge2d, TrimLoop2d};
