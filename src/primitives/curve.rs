use nalgebra::{Point3, Unit, Vector3};

/// Circle embedded in 3D space.
#[derive(Debug, Clone, PartialEq)]
pub struct Circle3 {
    /// Center of the circle.
    pub center: Point3<f64>,
    /// Unit normal of the plane containing the circle.
    pub normal: Unit<Vector3<f64>>,
    /// Radius in metres.
    pub radius: f64,
}

/// Ellipse embedded in 3D space.
///
/// Produced by plane × cylinder intersections at oblique angles.
#[derive(Debug, Clone, PartialEq)]
pub struct Ellipse3 {
    /// Center of the ellipse.
    pub center: Point3<f64>,
    /// Unit normal of the plane containing the ellipse.
    pub normal: Unit<Vector3<f64>>,
    /// Unit direction of the major axis (lies in the ellipse plane).
    pub major_dir: Unit<Vector3<f64>>,
    /// Semi-major axis length in metres.
    pub semi_major: f64,
    /// Semi-minor axis length in metres.
    pub semi_minor: f64,
}

impl Ellipse3 {
    /// Point at parametric angle `t` (radians): `c + a·cos(t)·u + b·sin(t)·v`
    /// where `u` is the major direction and `v = normal × u`.
    pub fn point_at(&self, t: f64) -> Point3<f64> {
        let u = self.major_dir.into_inner();
        let v = self.normal.cross(&self.major_dir);
        self.center + u * (self.semi_major * t.cos()) + v * (self.semi_minor * t.sin())
    }
}
