//! Curved panel primitives with UV-space trim loops.
//!
//! This module is intentionally separate from the main B-rep/CSG evaluator. It
//! starts the curved-opening path as a trimmed surface primitive: holes are
//! supplied as loops in the surface parameter domain and tessellated for display.
//! General surface-surface boolean evaluation remains out of scope.

mod cylinder_panel;
mod error;
mod mesh;
mod trim;

pub use cylinder_panel::{tessellate_cylinder_panel, CylinderPanel, CylinderPanelOptions};
pub use error::CurvedError;
pub use mesh::SurfaceMesh;
pub use trim::{TrimEdge2d, TrimLoop2d};
