//! Error type for the 2-D polygon boolean engine.
//!
//! Hand-written `Display` / `Error` impls — no `thiserror`, matching the parent
//! kernel's "minimal dependency" policy.

use std::fmt;

/// Errors returned by the public boolean operations.
///
/// `#[non_exhaustive]` so new failure modes (e.g. complexity limits) can be
/// added without a breaking change, matching the parent kernel's enum policy.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Poly2Error {
    /// An input edge is a circular arc. Arc support is structurally present in
    /// the edge model but not yet implemented; the engine fails fast rather than
    /// silently mis-handling the arc.
    ArcNotYetSupported,
    /// An input contour has fewer than three distinct vertices, so it bounds no
    /// area. The offending contour index (within its region) is reported.
    DegenerateContour {
        /// Index of the contour in the input region.
        contour_index: usize,
    },
    /// An input contour is non-simple: two of its edges cross at a point that is
    /// not a shared endpoint. Self-intersecting input is outside the engine's
    /// contract (the output of a boolean op is always simple; its inputs must be
    /// too). Reports an approximate location of the crossing.
    SelfIntersectingInput {
        /// Approximate x of the crossing.
        x: f64,
        /// Approximate y of the crossing.
        y: f64,
    },
    /// An internal invariant of the arrangement was violated (e.g. a half-edge
    /// failed to find its successor). This indicates a bug, not bad input; it is
    /// surfaced as a `Result` rather than a panic so the engine stays
    /// panic-free at the public boundary.
    Internal {
        /// Human-readable description of the violated invariant.
        what: &'static str,
    },
    /// A circular-arc degeneracy outside the supported set was encountered
    /// (`DESIGN.md` §4.2, §13-3). Phase 3c supports the rectangle × circular-void
    /// family and the everyday circle/circle cases (separate, overlapping,
    /// tangent, concentric); a genuinely ambiguous tangent / nested configuration
    /// that the closed-form path cannot resolve is reported here rather than
    /// silently mis-answered. The string names the specific degeneracy.
    UnsupportedArcDegeneracy {
        /// Human-readable description of the unsupported degeneracy.
        what: &'static str,
    },
}

impl fmt::Display for Poly2Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Poly2Error::ArcNotYetSupported => {
                write!(f, "arc edges are not yet supported (only line segments)")
            }
            Poly2Error::DegenerateContour { contour_index } => {
                write!(
                    f,
                    "contour {contour_index} is degenerate (fewer than 3 distinct vertices)"
                )
            }
            Poly2Error::SelfIntersectingInput { x, y } => {
                write!(f, "self-intersecting input near ({x}, {y})")
            }
            Poly2Error::Internal { what } => {
                write!(f, "internal arrangement invariant violated: {what}")
            }
            Poly2Error::UnsupportedArcDegeneracy { what } => {
                write!(f, "unsupported circular-arc degeneracy: {what}")
            }
        }
    }
}

impl std::error::Error for Poly2Error {}
