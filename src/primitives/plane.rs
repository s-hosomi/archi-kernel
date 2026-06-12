use nalgebra::{Point3, Unit, Vector3};

/// Infinite plane defined by a point and a unit normal.
#[derive(Debug, Clone, PartialEq)]
pub struct Plane {
    /// A point on the plane.
    pub point: Point3<f64>,
    /// Unit normal. The side the normal points to is the "above" side for
    /// signed distances.
    pub normal: Unit<Vector3<f64>>,
}

impl Plane {
    /// Build a plane from a point and a (not necessarily unit) normal.
    ///
    /// # Panics
    /// Panics if `normal` has zero length.
    pub fn new(point: Point3<f64>, normal: Vector3<f64>) -> Self {
        assert!(normal.norm() > 0.0, "plane normal must be non-zero");
        Self {
            point,
            normal: Unit::new_normalize(normal),
        }
    }

    /// Signed distance from `p` to the plane (positive on the normal side).
    pub fn signed_distance(&self, p: &Point3<f64>) -> f64 {
        self.normal.dot(&(p - self.point))
    }

    /// Orthogonal projection of `p` onto the plane.
    pub fn project_point(&self, p: &Point3<f64>) -> Point3<f64> {
        p - self.normal.into_inner() * self.signed_distance(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_distance_respects_normal_side() {
        let plane = Plane::new(Point3::new(0.0_f64, 0.0_f64, 2.0_f64), Vector3::z());
        assert!(
            (plane.signed_distance(&Point3::new(5.0_f64, 1.0_f64, 3.0_f64)) - 1.0_f64).abs()
                < 1e-12_f64
        );
        assert!(
            (plane.signed_distance(&Point3::new(0.0_f64, 0.0_f64, 0.0_f64)) + 2.0_f64).abs()
                < 1e-12_f64
        );
    }

    #[test]
    fn project_point_lands_on_plane() {
        let plane = Plane::new(Point3::origin(), Vector3::new(1.0_f64, 1.0_f64, 1.0_f64));
        let q = plane.project_point(&Point3::new(1.0_f64, 2.0_f64, 3.0_f64));
        assert!(plane.signed_distance(&q).abs() < 1e-12_f64);
    }
}
