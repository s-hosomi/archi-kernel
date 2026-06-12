use nalgebra::Point3;

use super::Line3;

/// Infinite circular cylinder surface defined by an axis and a radius.
///
/// Used for round columns, piles and circular voids. The surface is the set
/// of points at distance `radius` from `axis`.
#[derive(Debug, Clone, PartialEq)]
pub struct Cylinder {
    /// Cylinder axis.
    pub axis: Line3,
    /// Radius in metres (must be positive).
    pub radius: f64,
}

impl Cylinder {
    /// Build a cylinder from an axis and a radius.
    ///
    /// # Panics
    /// Panics if `radius` is not strictly positive.
    pub fn new(axis: Line3, radius: f64) -> Self {
        assert!(radius > 0.0, "cylinder radius must be positive");
        Self { axis, radius }
    }

    /// Signed distance from `p` to the surface (negative inside).
    pub fn signed_distance(&self, p: &Point3<f64>) -> f64 {
        self.axis.distance_to_point(p) - self.radius
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn signed_distance_negative_inside() {
        let cyl = Cylinder::new(Line3::new(Point3::origin(), Vector3::z()), 0.3_f64);
        assert!(cyl.signed_distance(&Point3::new(0.1_f64, 0.0_f64, 5.0_f64)) < 0.0_f64);
        assert!(cyl.signed_distance(&Point3::new(1.0_f64, 0.0_f64, -2.0_f64)) > 0.0_f64);
        assert!(
            cyl.signed_distance(&Point3::new(0.3_f64, 0.0_f64, 1.0_f64))
                .abs()
                < 1e-12_f64
        );
    }
}
