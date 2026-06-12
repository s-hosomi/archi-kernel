//! Error types for the archi-kernel public API.
//!
//! Errors are plain enums with hand-written `Display` and `std::error::Error`
//! implementations. The `thiserror` crate is intentionally not used (see
//! DESIGN.md §8 dependency policy).

use std::fmt;

/// Errors that can be produced by kernel constructors and operations.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum KernelError {
    /// A plane or circle was constructed with a zero-length normal vector.
    ZeroNormal,
    /// A line or cylinder axis was constructed with a zero-length direction.
    ZeroDirection,
    /// A radius or semi-axis value was not strictly positive.
    NonPositiveRadius {
        /// The offending radius value.
        radius: f64,
    },
    /// A semi-axis value was not strictly positive (used for ellipses).
    NonPositiveSemiAxis {
        /// The offending semi-axis value.
        semi_axis: f64,
    },
    /// The semi-major axis is smaller than the semi-minor axis.
    SemiMajorLessThanSemiMinor {
        /// Semi-major axis value.
        semi_major: f64,
        /// Semi-minor axis value.
        semi_minor: f64,
    },
    /// A profile dimension (width, height, …) was not strictly positive.
    NonPositiveDimension {
        /// The name of the offending dimension.
        name: &'static str,
        /// The offending value.
        value: f64,
    },
    /// The `major_dir` supplied to `Ellipse3::new` is (nearly) parallel to
    /// `normal`, so the two directions cannot span the ellipse plane.
    ///
    /// The dot product `|major_dir · normal|` must be strictly less than
    /// `Tol::angular` for a valid ellipse (orthogonality invariant).
    MajorDirNotInPlane {
        /// Absolute value of `dot(major_dir, normal)` at the time of the check.
        dot: f64,
    },
}

impl fmt::Display for KernelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            KernelError::ZeroNormal => write!(f, "plane normal must be a non-zero vector"),
            KernelError::ZeroDirection => {
                write!(f, "line or axis direction must be a non-zero vector")
            }
            KernelError::NonPositiveRadius { radius } => {
                write!(f, "radius must be strictly positive, got {radius}")
            }
            KernelError::NonPositiveSemiAxis { semi_axis } => {
                write!(f, "semi-axis must be strictly positive, got {semi_axis}")
            }
            KernelError::SemiMajorLessThanSemiMinor {
                semi_major,
                semi_minor,
            } => write!(
                f,
                "semi_major ({semi_major}) must be >= semi_minor ({semi_minor})"
            ),
            KernelError::NonPositiveDimension { name, value } => {
                write!(
                    f,
                    "profile dimension {name} must be strictly positive, got {value}"
                )
            }
            KernelError::MajorDirNotInPlane { dot } => {
                write!(
                    f,
                    "major_dir must be perpendicular to normal (|dot| must be < angular tol), \
                     got |dot| = {dot}"
                )
            }
        }
    }
}

impl std::error::Error for KernelError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_zero_normal() {
        let msg = KernelError::ZeroNormal.to_string();
        assert!(msg.contains("normal"));
    }

    #[test]
    fn display_non_positive_radius() {
        let msg = KernelError::NonPositiveRadius { radius: -1.0_f64 }.to_string();
        assert!(msg.contains("-1"));
    }

    #[test]
    fn error_trait_is_implemented() {
        let e: &dyn std::error::Error = &KernelError::ZeroDirection;
        assert!(e.source().is_none());
    }
}
