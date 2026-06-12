use crate::error::KernelError;
use crate::math::{Point3, Unit3, Vec3};

/// Infinite line defined by an origin point and a unit direction.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Line3 {
    origin: Point3,
    dir: Unit3,
}

impl Line3 {
    /// Build a line from a point and a (not necessarily unit) direction.
    ///
    /// Returns [`KernelError::ZeroDirection`] if `dir` has zero length.
    pub fn new(origin: Point3, dir: Vec3) -> Result<Self, KernelError> {
        let unit = dir.try_unit().ok_or(KernelError::ZeroDirection)?;
        Ok(Self { origin, dir: unit })
    }

    /// Build a line from a point and a pre-validated unit direction.
    ///
    /// Intended for internal use where the invariant is already established.
    pub(crate) fn new_unchecked(origin: Point3, dir: Unit3) -> Self {
        Self { origin, dir }
    }

    /// A point on the line (the parameterization origin).
    #[inline]
    pub fn origin(self) -> Point3 {
        self.origin
    }

    /// Unit direction of the line.
    #[inline]
    pub fn dir(self) -> Unit3 {
        self.dir
    }

    /// Point at parameter `t` (arc-length parameterization since `dir` is unit).
    pub fn point_at(self, t: f64) -> Point3 {
        self.origin + self.dir.as_vec() * t
    }

    /// Shortest distance from `p` to the line.
    pub fn distance_to_point(self, p: Point3) -> f64 {
        let v = p - self.origin;
        v.cross(self.dir.as_vec()).norm()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_to_point_is_perpendicular_distance() {
        let line = Line3::new(Point3::origin(), Vec3::Z).expect("valid line");
        let d = line.distance_to_point(Point3::new(3.0_f64, 4.0_f64, 10.0_f64));
        assert!((d - 5.0_f64).abs() < 1e-12_f64, "d = {d}");
    }

    #[test]
    fn zero_direction_returns_error() {
        assert_eq!(
            Line3::new(Point3::origin(), Vec3::ZERO),
            Err(KernelError::ZeroDirection)
        );
    }
}
