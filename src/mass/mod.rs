//! Mass properties (`DESIGN.md` §6-4).
//!
//! The signed [`volume`](signed_volume) is the orientation check Phase 2 relies
//! on (a correctly built outward-oriented solid has `V > 0`) and the concrete
//! quantity the take-off needs. Phase 5 completes the suite with [`centroid`],
//! total [`surface_area`], and the [`formwork_area`] split (side vs bottom) that
//! is the heart of the quantity take-off. Every entry point is closed-form and
//! returns a `Result` so an unsupported face is surfaced, never silently zero.

mod face_area;
mod properties;
mod volume;

pub use face_area::AreaError;
pub use properties::{centroid, formwork_area, surface_area, CentroidError, FormworkArea};
pub use volume::{signed_volume, signed_volume_checked, VolumeError};

/// The per-face net area (outer minus holes), shared by the take-off's
/// contact-aware formwork. Crate-internal; the public entry is
/// [`formwork_area`].
pub(crate) fn face_area_of(
    brep: &crate::brep::Brep,
    face: &crate::topo::Face,
) -> Result<f64, AreaError> {
    face_area::face_area(brep, face)
}
