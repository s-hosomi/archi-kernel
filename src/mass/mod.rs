//! Mass properties (Phase 5, minimal early implementation).
//!
//! Phase 2 needs a signed volume purely as an *orientation* check: a correctly
//! built solid with outward-facing normals has strictly positive volume, so the
//! extrusion tests assert `signed_volume > 0`. The full mass-property suite
//! (centroid, section properties, the `V(A−B) = V(A) − V(A∩B)` identity under
//! proptest) arrives in Phase 5; only the divergence-theorem volume lives here
//! for now.

mod volume;

pub use volume::signed_volume;
