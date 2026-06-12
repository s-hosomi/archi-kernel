use nalgebra::{Point3, Unit, Vector3};

/// Infinite line defined by an origin point and a unit direction.
#[derive(Debug, Clone, PartialEq)]
pub struct Line3 {
    /// A point on the line.
    pub origin: Point3<f64>,
    /// Unit direction of the line.
    pub dir: Unit<Vector3<f64>>,
}

impl Line3 {
    /// Build a line from a point and a (not necessarily unit) direction.
    ///
    /// # Panics
    /// Panics if `dir` has zero length.
    pub fn new(origin: Point3<f64>, dir: Vector3<f64>) -> Self {
        assert!(dir.norm() > 0.0, "line direction must be non-zero");
        Self {
            origin,
            dir: Unit::new_normalize(dir),
        }
    }

    /// Point at parameter `t` (arc-length parameterization since `dir` is unit).
    pub fn point_at(&self, t: f64) -> Point3<f64> {
        self.origin + self.dir.into_inner() * t
    }

    /// Shortest distance from `p` to the line.
    pub fn distance_to_point(&self, p: &Point3<f64>) -> f64 {
        let v = p - self.origin;
        v.cross(&self.dir).norm()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_to_point_is_perpendicular_distance() {
        let line = Line3::new(Point3::origin(), Vector3::z());
        let d = line.distance_to_point(&Point3::new(3.0_f64, 4.0_f64, 10.0_f64));
        assert!((d - 5.0_f64).abs() < 1e-12_f64, "d = {d}");
    }
}
