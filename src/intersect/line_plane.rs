use crate::math::Point3;
use crate::primitives::{Line3, Plane};
use crate::tolerance::Tol;

/// Intersection point of a line and a plane.
///
/// Returns `None` if the line is parallel to the plane (within angular
/// tolerance), including the coincident case.
pub fn line_plane(line: &Line3, plane: &Plane, tol: &Tol) -> Option<Point3> {
    let denom = plane.normal().dot(line.dir().as_vec());
    if denom.abs() <= tol.angular {
        return None;
    }
    let t = plane.normal().dot(plane.point() - line.origin()) / denom;
    Some(line.point_at(t))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn line_hits_plane_at_expected_point() {
        let tol = Tol::default();
        let plane =
            Plane::new(Point3::new(0.0_f64, 0.0_f64, 3.0_f64), Vec3::Z).expect("valid plane");
        let line = Line3::new(Point3::new(1.0_f64, 2.0_f64, 0.0_f64), Vec3::Z).expect("valid line");
        let p = line_plane(&line, &plane, &tol).expect("must intersect");
        assert!(
            (p - Point3::new(1.0_f64, 2.0_f64, 3.0_f64)).norm() < 1e-12_f64,
            "p = {p:?}"
        );
    }

    #[test]
    fn parallel_line_returns_none() {
        let tol = Tol::default();
        let plane = Plane::new(Point3::origin(), Vec3::Z).expect("valid plane");
        let line = Line3::new(Point3::new(0.0_f64, 0.0_f64, 1.0_f64), Vec3::X).expect("valid line");
        assert!(line_plane(&line, &plane, &tol).is_none());
    }
}
