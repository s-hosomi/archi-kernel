//! # poly2d — a robust 2-D polygon boolean engine
//!
//! This module is the "second kernel" of archi-kernel: a robust 2-D polygon
//! boolean engine over [`Region`]s (outer boundaries + holes, multiple
//! connected components). The overwhelming majority of 3-D booleans in the
//! building domain (prismatic differences) reduce to 2-D and land here.
//!
//! Integrated from the self-contained `poly2d-dev/` scratch package. The
//! local `tol.rs` mirror has been deleted; the parent `crate::tolerance::Tol`
//! and `crate::tolerance::Sign3` are used directly.  The `orient2d` predicate
//! is routed through `crate::predicate::orient2d_exact` in compliance with the
//! kernel's isolation rule (DESIGN.md §3.5): the `robust` crate is confined to
//! `src/predicate/` and must not be imported elsewhere.
//!
//! ## Algorithm choice — snapped planar arrangement + winding classification
//!
//! Following the kernel's boolean research (docs/research/05-boolean.md), the
//! engine uses a **snapped planar arrangement with face winding classification**
//! rather than a vertex-pair-walking method (Greiner–Hormann / Weiler–Atherton):
//!
//! 1. **Snap / merge** all vertices and intersection points with one absolute
//!    `eps`, so the pervasive "exact coincidence" of building geometry (shared
//!    edges, vertex-on-edge, opening touching a wall) collapses to identical
//!    vertices *before* any topology is built.
//! 2. **Split** every edge at its intersection vertices and **dedup** collinear
//!    coincident fragments into single arrangement edges, tracking each
//!    operand's directed winding contribution. Shared edges whose windings
//!    cancel carry no boundary and vanish — the heart of degeneracy handling.
//! 3. **Build a half-edge arrangement (DCEL)** and extract its bounded faces.
//! 4. **Classify** each face by its winding number in each operand (taken at an
//!    interior representative point, away from every edge, so the classification
//!    is unambiguous), and select with the operation's keep table.
//! 5. **Reconstruct** the result [`Region`], cancelling shared edges between
//!    selected faces and tracing outer/hole loops with nesting.
//!
//! **Why this family.** Greiner–Hormann and Weiler–Atherton trace intersection
//! pairs around the polygons; they are structurally fragile under degeneracy
//! (collinear overlaps, vertex-on-edge, identical inputs), which is the *common*
//! case in buildings, not the exception. The arrangement-plus-winding approach
//! never traces "the other polygon" — it builds one combined subdivision and
//! asks a local, exact question per face. Combined with snap-merge for
//! coincidence collapse and Shewchuk's exact `orient2d` for every left/right
//! decision, it matches the kernel's "Manifold-leaning" posture: single absolute
//! tolerance, `On` as a first-class state, degeneracies *detected and merged*
//! rather than perturbed away (no SoS).
//!
//! ## Edge model
//!
//! Edges are [`Edge2`] `{ Seg, Arc }` from the start. Only `Seg` is implemented;
//! any `Arc` yields [`Poly2Error::ArcNotYetSupported`]. The intersection
//! dispatch, point-on-edge classification, and face classification are all
//! structured so arc support is *added*, never retrofitted.
//!
//! ## Scope and conventions
//!
//! * Operations: [`difference`], [`union`], [`intersection`] — one shared
//!   classifier, three keep tables ([`classify::Op`]).
//! * Input/output: [`Region`] (outer CCW, holes CW; orientation is normalized on
//!   output via signed area).
//! * Single absolute length tolerance ([`crate::tolerance::Tol::length`], default
//!   `1e-6` m). No scale adaptation.
//! * Public API is panic-free: every entry point returns `Result`.

#![forbid(unsafe_code)]

mod arrangement;
mod classify;
mod error;
pub(crate) mod geom;
pub(crate) mod intersect;
mod reconstruct;
mod region;
pub(crate) mod snap;

pub use classify::Op;
pub use error::Poly2Error;
pub use geom::{Arc, Edge2, Orient, Point2, Vec2};
pub use region::{Contour, Region};

use crate::tolerance::Tol;
use arrangement::{Arrangement, Operand};
use classify::inside;

/// Compute `A − B` (the part of `a` not covered by `b`).
///
/// # Errors
/// Returns [`Poly2Error::ArcNotYetSupported`] if any edge is an arc, and the
/// validation errors for malformed input (degenerate / self-intersecting).
pub fn difference(a: &Region, b: &Region, tol: &Tol) -> Result<Region, Poly2Error> {
    boolean(a, b, Op::Difference, tol)
}

/// Compute `A ∪ B`.
///
/// # Errors
/// See [`difference`].
pub fn union(a: &Region, b: &Region, tol: &Tol) -> Result<Region, Poly2Error> {
    boolean(a, b, Op::Union, tol)
}

/// Compute `A ∩ B`.
///
/// # Errors
/// See [`difference`].
pub fn intersection(a: &Region, b: &Region, tol: &Tol) -> Result<Region, Poly2Error> {
    boolean(a, b, Op::Intersection, tol)
}

/// Shared driver: build the arrangement, classify faces with the op's keep
/// table, and reconstruct the result region.
fn boolean(a: &Region, b: &Region, op: Op, tol: &Tol) -> Result<Region, Poly2Error> {
    // Fast paths for empty operands keep the engine cheap and well-defined on
    // the identities the callers rely on.
    if a.is_empty() && b.is_empty() {
        return Ok(Region::empty());
    }
    match op {
        Op::Intersection if a.is_empty() || b.is_empty() => return Ok(Region::empty()),
        Op::Difference if a.is_empty() => return Ok(Region::empty()),
        _ => {}
    }

    let arr = Arrangement::build(a, b, tol)?;

    let mut selected: Vec<Vec<snap::VertexId>> = Vec::new();
    for face in &arr.faces {
        // Classify the face that lies to the *left* of this directed loop. A
        // kept CCW loop is an outer boundary; a kept CW loop is a hole boundary
        // of the same selected face. The unbounded outer wrap classifies as
        // outside both operands and is never kept.
        let p = face.face_sample_point();
        let in_a = inside(arr.winding(p, Operand::A));
        let in_b = inside(arr.winding(p, Operand::B));
        if op.keep(in_a, in_b) {
            selected.push(face.vertex_ids.clone());
        }
    }

    Ok(reconstruct::reconstruct(arr.store(), &selected, tol))
}
