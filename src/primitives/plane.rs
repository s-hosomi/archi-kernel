use crate::error::KernelError;
use crate::math::{Point3, Unit3, Vec3};

/// Infinite plane defined by a point and a unit normal.
///
/// The side the normal points to is the "above" side for signed distances.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Plane {
    point: Point3,
    normal: Unit3,
}

impl Plane {
    /// Build a plane from a point and a (not necessarily unit) normal vector.
    ///
    /// Returns [`KernelError::ZeroNormal`] if `normal` has zero length.
    pub fn new(point: Point3, normal: Vec3) -> Result<Self, KernelError> {
        let unit = normal.try_unit().ok_or(KernelError::ZeroNormal)?;
        Ok(Self {
            point,
            normal: unit,
        })
    }

    /// Build a plane from a point and a pre-validated unit normal.
    ///
    /// Intended for internal use where the caller already guarantees `normal`
    /// is a unit vector (e.g. inside the intersect module).
    #[allow(dead_code)]
    pub(crate) fn new_unchecked(point: Point3, normal: Unit3) -> Self {
        Self { point, normal }
    }

    /// A point on the plane.
    #[inline]
    pub fn point(self) -> Point3 {
        self.point
    }

    /// Unit normal of the plane.
    #[inline]
    pub fn normal(self) -> Unit3 {
        self.normal
    }

    /// Signed distance from `p` to the plane (positive on the normal side).
    pub fn signed_distance(self, p: Point3) -> f64 {
        self.normal.dot(p - self.point)
    }

    /// Orthogonal projection of `p` onto the plane.
    pub fn project_point(self, p: Point3) -> Point3 {
        p - self.normal.as_vec() * self.signed_distance(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_distance_respects_normal_side() {
        let plane =
            Plane::new(Point3::new(0.0_f64, 0.0_f64, 2.0_f64), Vec3::Z).expect("valid plane");
        assert!(
            (plane.signed_distance(Point3::new(5.0_f64, 1.0_f64, 3.0_f64)) - 1.0_f64).abs()
                < 1e-12_f64
        );
        assert!(
            (plane.signed_distance(Point3::new(0.0_f64, 0.0_f64, 0.0_f64)) + 2.0_f64).abs()
                < 1e-12_f64
        );
    }

    #[test]
    fn project_point_lands_on_plane() {
        let plane = Plane::new(Point3::origin(), Vec3::new(1.0_f64, 1.0_f64, 1.0_f64))
            .expect("valid plane");
        let q = plane.project_point(Point3::new(1.0_f64, 2.0_f64, 3.0_f64));
        assert!(plane.signed_distance(q).abs() < 1e-12_f64);
    }

    #[test]
    fn zero_normal_returns_error() {
        assert_eq!(
            Plane::new(Point3::origin(), Vec3::ZERO),
            Err(KernelError::ZeroNormal)
        );
    }
}
