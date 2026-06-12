//! Stable identifiers for members and openings.
//!
//! These are opaque `u64` newtypes assigned by the *caller* (the adapter or
//! application layer); the kernel never invents them. Using stable ids — rather
//! than B-rep face indices — for external references structurally avoids the
//! topological-naming problem (`DESIGN.md` §5.1).

/// A stable identifier for a member (column, beam, wall, …).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StableId(pub u64);

/// A stable identifier for an opening (an `IfcRelVoidsElement`-equivalent void).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OpeningId(pub u64);
