//! Errors for curved panel construction and tessellation.

use std::fmt;

/// A curved-panel input or operation could not be accepted.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum CurvedError {
    /// A parameter range was not strictly increasing and finite.
    InvalidRange {
        /// Name of the offending parameter range.
        name: &'static str,
        /// Lower endpoint.
        min: f64,
        /// Upper endpoint.
        max: f64,
    },
    /// A chord tolerance was not strictly positive and finite.
    NonPositiveChordTolerance {
        /// The offending value.
        value: f64,
    },
    /// A panel thickness was not strictly positive and finite.
    NonPositiveThickness {
        /// The offending value.
        value: f64,
    },
    /// A radius was not strictly positive and finite.
    NonPositiveRadius {
        /// The offending value.
        value: f64,
    },
    /// The inner offset radius of a thick cylindrical panel is not positive.
    NonPositiveInnerRadius {
        /// The computed inner radius.
        radius: f64,
    },
    /// A trim loop has no edges.
    EmptyLoop,
    /// A trim edge has invalid parameters.
    InvalidTrimEdge,
    /// Consecutive trim edges do not connect within tolerance.
    OpenLoop,
    /// The trim loop's signed area is too small to define a stable region.
    DegenerateLoop {
        /// Signed UV-space area.
        area: f64,
    },
    /// Arc trim edges are represented but not supported by the requested
    /// operation.
    UnsupportedArcTrim,
    /// A hole is not fully contained in the panel's outer UV rectangle.
    HoleOutsidePanel,
    /// Two hole loops overlap or cross.
    HoleOverlap,
    /// A trim loop crosses the cylinder parameter seam; this phase requires
    /// loops to live inside one unwrapped `theta` interval.
    SeamCrossing,
    /// A spherical panel includes a pole where the longitude parameter
    /// collapses.
    PoleCrossing,
}

impl fmt::Display for CurvedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CurvedError::InvalidRange { name, min, max } => {
                write!(
                    f,
                    "{name} range must be finite and increasing, got [{min}, {max}]"
                )
            }
            CurvedError::NonPositiveChordTolerance { value } => {
                write!(f, "chord tolerance must be strictly positive, got {value}")
            }
            CurvedError::NonPositiveThickness { value } => {
                write!(f, "panel thickness must be strictly positive, got {value}")
            }
            CurvedError::NonPositiveRadius { value } => {
                write!(f, "radius must be strictly positive, got {value}")
            }
            CurvedError::NonPositiveInnerRadius { radius } => {
                write!(f, "inner cylinder radius must stay positive, got {radius}")
            }
            CurvedError::EmptyLoop => write!(f, "trim loop must contain at least one edge"),
            CurvedError::InvalidTrimEdge => write!(f, "trim edge has invalid parameters"),
            CurvedError::OpenLoop => write!(f, "trim loop edges must form a closed chain"),
            CurvedError::DegenerateLoop { area } => {
                write!(f, "trim loop area is degenerate: {area}")
            }
            CurvedError::UnsupportedArcTrim => {
                write!(f, "arc trim edges are not tessellated in this phase")
            }
            CurvedError::HoleOutsidePanel => write!(f, "hole loop must lie inside the panel"),
            CurvedError::HoleOverlap => write!(f, "hole loops must not overlap or cross"),
            CurvedError::SeamCrossing => {
                write!(f, "trim loop must not cross the cylinder theta seam")
            }
            CurvedError::PoleCrossing => {
                write!(f, "spherical panel must not include a parameter pole")
            }
        }
    }
}

impl std::error::Error for CurvedError {}
