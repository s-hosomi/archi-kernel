//! Solid construction from cross-section profiles (`DESIGN.md` §6-1, Phase 2).
//!
//! The only constructor today is [`extrude`], which sweeps a 2-D
//! [`Profile2d`](crate::csg::Profile2d) along an axis to produce a closed,
//! validated [`Brep`](crate::brep::Brep). This is the direct geometric
//! counterpart of an ST-Bridge member (an axis plus a section reference).

mod extrude;

pub use extrude::extrude;
