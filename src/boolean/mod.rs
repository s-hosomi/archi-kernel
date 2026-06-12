//! Boolean machinery — the kernel's first real boolean mechanism.
//!
//! Phase 3a implements the easiest of the classical four boolean stages applied
//! to the simplest operand: **solid × half-space cut** ([`half_space::cut`]).
//! Cutting is the shared foundation for the difference boolean (Phase 3b) and
//! the section drawing (Phase 4): both reuse the inside/outside classifier and
//! the face splitter that the cut produces (`DESIGN.md` §4.1,
//! `docs/research/05-boolean.md` §5).
//!
//! A cut runs the classical pipeline restricted to *Imprint → Classify → Cap*:
//!
//! 1. classify every vertex against the cutting plane (`predicate::side_of_plane`);
//! 2. split the half-edges that straddle the plane at their closed-form
//!    intersection point (the new point is `On`);
//! 3. split each face by the plane, keeping the wanted side and recording the
//!    `On` boundary segments;
//! 4. chain those `On` segments into closed loops and cap the opening, with hole
//!    loops handled so a through-hole solid yields an annulus cap;
//! 5. validate the result at [`ValidateLevel::Full`](crate::topo::ValidateLevel).
//!
//! The output is always a fresh [`Brep`](crate::brep::Brep) with its own
//! geometry store (the input is read-only — `DESIGN.md` §5.2 ownership model).

pub mod half_space;
pub mod poly2d;
pub mod prismatic;

mod support;

pub use half_space::{cut, CutResult, KeepSide};
