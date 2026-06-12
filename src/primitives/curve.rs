use crate::error::KernelError;
use crate::math::{Point3, Unit3, Vec3};

/// Circle embedded in 3-D space.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Circle3 {
    center: Point3,
    normal: Unit3,
    radius: f64,
}

impl Circle3 {
    /// Build a circle from a center, a unit normal, and a radius.
    ///
    /// Returns [`KernelError::ZeroNormal`] if `normal` has zero length, or
    /// [`KernelError::NonPositiveRadius`] if `radius` is not strictly positive.
    pub fn new(center: Point3, normal: Vec3, radius: f64) -> Result<Self, KernelError> {
        let unit = normal.try_unit().ok_or(KernelError::ZeroNormal)?;
        if radius <= 0.0 {
            return Err(KernelError::NonPositiveRadius { radius });
        }
        Ok(Self {
            center,
            normal: unit,
            radius,
        })
    }

    /// Build a `Circle3` without checking the invariants.
    ///
    /// Intended for internal use where the invariants are already established.
    pub(crate) fn new_unchecked(center: Point3, normal: Unit3, radius: f64) -> Self {
        Self {
            center,
            normal,
            radius,
        }
    }

    /// Center of the circle.
    #[inline]
    pub fn center(self) -> Point3 {
        self.center
    }

    /// Unit normal of the plane containing the circle.
    #[inline]
    pub fn normal(self) -> Unit3 {
        self.normal
    }

    /// Radius of the circle in metres.
    #[inline]
    pub fn radius(self) -> f64 {
        self.radius
    }
}

/// Ellipse embedded in 3-D space.
///
/// Produced by plane × cylinder intersections at oblique angles.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Ellipse3 {
    center: Point3,
    normal: Unit3,
    major_dir: Unit3,
    semi_major: f64,
    semi_minor: f64,
}

impl Ellipse3 {
    /// Build an ellipse.
    ///
    /// `semi_major` and `semi_minor` must both be strictly positive, and
    /// `semi_major >= semi_minor`. Returns an appropriate [`KernelError`]
    /// otherwise.
    pub fn new(
        center: Point3,
        normal: Vec3,
        major_dir: Vec3,
        semi_major: f64,
        semi_minor: f64,
    ) -> Result<Self, KernelError> {
        let unit_normal = normal.try_unit().ok_or(KernelError::ZeroNormal)?;
        let unit_major = major_dir.try_unit().ok_or(KernelError::ZeroDirection)?;
        if semi_minor <= 0.0 {
            return Err(KernelError::NonPositiveSemiAxis {
                semi_axis: semi_minor,
            });
        }
        if semi_major <= 0.0 {
            return Err(KernelError::NonPositiveSemiAxis {
                semi_axis: semi_major,
            });
        }
        if semi_major < semi_minor {
            return Err(KernelError::SemiMajorLessThanSemiMinor {
                semi_major,
                semi_minor,
            });
        }
        Ok(Self {
            center,
            normal: unit_normal,
            major_dir: unit_major,
            semi_major,
            semi_minor,
        })
    }

    /// Build an `Ellipse3` without checking the invariants.
    ///
    /// Intended for internal use where the invariants are already established.
    pub(crate) fn new_unchecked(
        center: Point3,
        normal: Unit3,
        major_dir: Unit3,
        semi_major: f64,
        semi_minor: f64,
    ) -> Self {
        Self {
            center,
            normal,
            major_dir,
            semi_major,
            semi_minor,
        }
    }

    /// Center of the ellipse.
    #[inline]
    pub fn center(self) -> Point3 {
        self.center
    }

    /// Unit normal of the plane containing the ellipse.
    #[inline]
    pub fn normal(self) -> Unit3 {
        self.normal
    }

    /// Unit direction of the major axis (lies in the ellipse plane).
    #[inline]
    pub fn major_dir(self) -> Unit3 {
        self.major_dir
    }

    /// Semi-major axis length in metres.
    #[inline]
    pub fn semi_major(self) -> f64 {
        self.semi_major
    }

    /// Semi-minor axis length in metres.
    #[inline]
    pub fn semi_minor(self) -> f64 {
        self.semi_minor
    }

    /// Point at parametric angle `t` (radians): `c + a·cos(t)·u + b·sin(t)·v`
    /// where `u` is the major direction and `v = normal × u`.
    pub fn point_at(self, t: f64) -> Point3 {
        let u = self.major_dir.as_vec();
        let v = self.normal.cross(self.major_dir.as_vec());
        self.center + u * (self.semi_major * t.cos()) + v * (self.semi_minor * t.sin())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circle3_zero_normal_error() {
        assert_eq!(
            Circle3::new(Point3::origin(), Vec3::ZERO, 1.0_f64),
            Err(KernelError::ZeroNormal)
        );
    }

    #[test]
    fn circle3_non_positive_radius_error() {
        assert!(matches!(
            Circle3::new(Point3::origin(), Vec3::Z, 0.0_f64),
            Err(KernelError::NonPositiveRadius { .. })
        ));
    }

    #[test]
    fn ellipse3_semi_major_less_than_semi_minor_error() {
        assert!(matches!(
            Ellipse3::new(Point3::origin(), Vec3::Z, Vec3::X, 0.5_f64, 1.0_f64),
            Err(KernelError::SemiMajorLessThanSemiMinor { .. })
        ));
    }

    #[test]
    fn ellipse3_point_at_lies_on_correct_axes() {
        let e = Ellipse3::new(Point3::origin(), Vec3::Z, Vec3::X, 2.0_f64, 1.0_f64)
            .expect("valid ellipse");
        // t = 0 → point along major axis at distance semi_major
        let p = e.point_at(0.0_f64);
        assert!((p.x - 2.0_f64).abs() < 1e-12_f64);
        assert!(p.y.abs() < 1e-12_f64);
        assert!(p.z.abs() < 1e-12_f64);
    }
}
