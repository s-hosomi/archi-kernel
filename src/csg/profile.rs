//! 2-D cross-section profiles (minimal placeholder).
//!
//! Only the shape vocabulary and positive-dimension validation are defined
//! here; the full profile geometry (the actual boundary loops) is built in
//! Phase 2. Constructors are panic-free and validate that dimensions are
//! strictly positive.

use crate::error::KernelError;

/// A 2-D cross-section to be extruded.
///
/// Dimensions are half-extents / radii in metres. The variants cover the
/// building cases (`DESIGN.md` §6): rectangular and H-section members and round
/// columns. This is a minimal stand-in completed in Phase 2.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum Profile2d {
    /// A rectangle, given by its half-width and half-height.
    Rect {
        /// Half the rectangle's width (along the local x axis), in metres.
        half_w: f64,
        /// Half the rectangle's height (along the local y axis), in metres.
        half_h: f64,
    },
    /// An I/H section, given by overall half-width and half-height, web
    /// thickness and flange thickness.
    HSection {
        /// Half the overall section width, in metres.
        half_w: f64,
        /// Half the overall section height, in metres.
        half_h: f64,
        /// Web thickness, in metres.
        web: f64,
        /// Flange thickness, in metres.
        flange: f64,
    },
    /// A circle, given by its radius.
    Circle {
        /// Radius, in metres.
        radius: f64,
    },
}

impl Profile2d {
    /// Build a rectangle from its half-width and half-height.
    ///
    /// Returns [`KernelError::NonPositiveDimension`] if either is not strictly
    /// positive.
    pub fn rect(half_w: f64, half_h: f64) -> Result<Self, KernelError> {
        positive("half_w", half_w)?;
        positive("half_h", half_h)?;
        Ok(Profile2d::Rect { half_w, half_h })
    }

    /// Build an H-section.
    ///
    /// Returns [`KernelError::NonPositiveDimension`] if any dimension is not
    /// strictly positive.
    pub fn h_section(half_w: f64, half_h: f64, web: f64, flange: f64) -> Result<Self, KernelError> {
        positive("half_w", half_w)?;
        positive("half_h", half_h)?;
        positive("web", web)?;
        positive("flange", flange)?;
        Ok(Profile2d::HSection {
            half_w,
            half_h,
            web,
            flange,
        })
    }

    /// Build a circle from its radius.
    ///
    /// Returns [`KernelError::NonPositiveDimension`] if the radius is not
    /// strictly positive.
    pub fn circle(radius: f64) -> Result<Self, KernelError> {
        positive("radius", radius)?;
        Ok(Profile2d::Circle { radius })
    }
}

fn positive(name: &'static str, value: f64) -> Result<(), KernelError> {
    if value > 0.0 {
        Ok(())
    } else {
        Err(KernelError::NonPositiveDimension { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rect_rejects_non_positive() {
        assert!(matches!(
            Profile2d::rect(0.0_f64, 1.0_f64),
            Err(KernelError::NonPositiveDimension { .. })
        ));
    }

    #[test]
    fn circle_accepts_positive() {
        assert_eq!(
            Profile2d::circle(0.3_f64),
            Ok(Profile2d::Circle { radius: 0.3_f64 })
        );
    }
}
