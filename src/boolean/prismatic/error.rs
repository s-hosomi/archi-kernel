//! Error type for the prismatic (2.5-D) boolean reduction.
//!
//! Hand-written `Display` / `Error` impls — no `thiserror`, matching the parent
//! kernel's minimal-dependency policy. Every variant is machine-readable: it
//! says *which* operand and *why* the reduction could not proceed
//! (`synthesis.md` §2-15), so a failing member can be isolated and reported.

use std::fmt;

use crate::boolean::poly2d::Poly2Error;
use crate::topo::validate::Defect;

/// Which operand of a binary prismatic boolean a failure refers to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum Operand {
    /// The positive / first operand (`A`).
    A,
    /// The negative / second operand (`B`).
    B,
}

impl fmt::Display for Operand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Operand::A => write!(f, "operand A"),
            Operand::B => write!(f, "operand B"),
        }
    }
}

/// Why a prismatic (2.5-D) reduction could not be carried out, or failed.
///
/// `#[non_exhaustive]` so new failure modes can be added in a semver-compatible
/// way (the parent kernel's enum policy). The variants are intentionally
/// fine-grained so that [`EvalError`](crate::csg::EvalError) can carry the exact
/// reason a member fell back to the unsupported / limit path.
///
/// This is an *internal* diagnostic: it is converted to the public
/// [`EvalError`](crate::csg::EvalError) (which carries plain strings / defects)
/// before it can reach a serde boundary, so it does not itself derive serde —
/// it nests [`Poly2Error`], whose `&'static str` payload cannot round-trip.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PrismError {
    /// No prismatic direction is common to both operands, so the difference does
    /// not reduce to 2.5-D. The caller surfaces this as
    /// [`EvalError::Unsupported3dBoolean`](crate::csg::EvalError).
    NoCommonDirection,
    /// An operand's extrusion axis is zero (or below `MIN_POSITIVE`), i.e. it has
    /// no direction. This is malformed input, not an unsupported configuration;
    /// the caller surfaces it as a construction error.
    DegenerateAxis,
    /// The operand involves a circular (cylindrical) cross-section whose curved
    /// boundary would land on the 2-D side as an arc; arc support in the 2-D
    /// engine arrives in Phase 3c. Reported with the offending operand.
    CircularInvolved {
        /// Which operand is circular.
        operand: Operand,
    },
    /// The 2-D engine reported an arc it cannot yet handle (Phase 3c). Equivalent
    /// in effect to [`CircularInvolved`](PrismError::CircularInvolved) but raised
    /// from inside the 2-D pipeline rather than at reprofiling.
    ArcNotYetSupported,
    /// The reprofiled region or the band structure exceeded the configured
    /// complexity budget; the member is isolated rather than ground through an
    /// unbounded computation (`DESIGN.md` §4.5).
    ComplexityLimit {
        /// The complexity measure that was hit (total boundary vertex count).
        measure: usize,
        /// The budget it exceeded.
        budget: usize,
    },
    /// The 2-D engine rejected an input or failed internally. The underlying
    /// [`Poly2Error`] is carried for diagnosis.
    Poly2(Poly2Error),
    /// The assembled 3-D result failed structural validation (`DESIGN.md` §4.2
    /// item 7, §7). The defects are carried so the offending member can be
    /// reported; the caller surfaces this as
    /// [`EvalError::InvalidResult`](crate::csg::EvalError).
    InvalidResult(Vec<Defect>),
}

impl fmt::Display for PrismError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PrismError::NoCommonDirection => {
                write!(f, "no common prismatic direction (not reducible to 2.5-D)")
            }
            PrismError::DegenerateAxis => {
                write!(f, "extrusion axis is degenerate (zero direction)")
            }
            PrismError::CircularInvolved { operand } => {
                write!(f, "{operand} has a circular cross-section (arc, Phase 3c)")
            }
            PrismError::ArcNotYetSupported => {
                write!(f, "the 2-D engine encountered an arc (Phase 3c)")
            }
            PrismError::ComplexityLimit { measure, budget } => {
                write!(
                    f,
                    "complexity {measure} exceeds the budget {budget} (member isolated)"
                )
            }
            PrismError::Poly2(e) => write!(f, "2-D boolean failed: {e}"),
            PrismError::InvalidResult(defects) => {
                write!(f, "assembled result is invalid ({} defects)", defects.len())
            }
        }
    }
}

impl std::error::Error for PrismError {}

impl From<Poly2Error> for PrismError {
    fn from(e: Poly2Error) -> Self {
        match e {
            Poly2Error::ArcNotYetSupported => PrismError::ArcNotYetSupported,
            other => PrismError::Poly2(other),
        }
    }
}
