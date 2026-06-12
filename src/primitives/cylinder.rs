use crate::error::KernelError;
use crate::math::Point3;
use crate::primitives::Line3;

/// Infinite circular cylinder surface defined by an axis and a radius.
///
/// Used for round columns, piles and circular voids. The surface is the set
/// of points at distance `radius` from `axis`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cylinder {
    axis: Line3,
    radius: f64,
}

impl Cylinder {
    /// Build a cylinder from an axis and a radius.
    ///
    /// Returns [`KernelError::NonPositiveRadius`] if `radius` is not strictly
    /// positive.
    pub fn new(axis: Line3, radius: f64) -> Result<Self, KernelError> {
        if radius <= 0.0 {
            return Err(KernelError::NonPositiveRadius { radius });
        }
        Ok(Self { axis, radius })
    }

    /// Build a cylinder without validating `radius`.
    ///
    /// Intended for internal use where the invariant is already established.
    #[allow(dead_code)]
    pub(crate) fn new_unchecked(axis: Line3, radius: f64) -> Self {
        Self { axis, radius }
    }

    /// The axis line of the cylinder.
    #[inline]
    pub fn axis(self) -> Line3 {
        self.axis
    }

    /// Radius of the cylinder in metres.
    #[inline]
    pub fn radius(self) -> f64 {
        self.radius
    }

    /// Signed distance from `p` to the surface (negative inside).
    pub fn signed_distance(self, p: Point3) -> f64 {
        self.axis.distance_to_point(p) - self.radius
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn signed_distance_negative_inside() {
        let cyl = Cylinder::new(Line3::new(Point3::origin(), Vec3::Z).unwrap(), 0.3_f64)
            .expect("valid cylinder");
        assert!(cyl.signed_distance(Point3::new(0.1_f64, 0.0_f64, 5.0_f64)) < 0.0_f64);
        assert!(cyl.signed_distance(Point3::new(1.0_f64, 0.0_f64, -2.0_f64)) > 0.0_f64);
        assert!(
            cyl.signed_distance(Point3::new(0.3_f64, 0.0_f64, 1.0_f64))
                .abs()
                < 1e-12_f64
        );
    }

    #[test]
    fn non_positive_radius_returns_error() {
        let axis = Line3::new(Point3::origin(), Vec3::Z).unwrap();
        assert!(matches!(
            Cylinder::new(axis, 0.0_f64),
            Err(KernelError::NonPositiveRadius { .. })
        ));
        assert!(matches!(
            Cylinder::new(axis, -1.0_f64),
            Err(KernelError::NonPositiveRadius { .. })
        ));
    }
}
