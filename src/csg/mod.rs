//! CSG history and lazy B-rep evaluation — the container only.
//!
//! Building members are represented as a shallow CSG tree ("sum of extrusions
//! minus openings"), which is the true source of truth; the B-rep is evaluated
//! lazily and cached (`DESIGN.md` §2.3, §5). This module defines the vocabulary
//! ([`CsgNode`], [`Opening`], [`ClipRule`], [`Profile2d`]), the per-member cache
//! and dirty machinery ([`Member`]), and the evaluation error type
//! ([`EvalError`]). Actual evaluation is a later phase; [`Member::brep`] returns
//! [`EvalError::NotYetImplemented`] for now.

mod ids;
mod member;
mod node;
mod profile;

pub use crate::boolean::prismatic::Operand;
pub use ids::{OpeningId, StableId};
pub use member::{EvalError, Member, UnsupportedReason};
pub use node::{ClipRule, CsgNode, Opening};
pub use profile::Profile2d;
