use nalgebra::Point3;

use crate::primitives::{Line3, Plane};
use crate::tolerance::Tol;

/// Result of intersecting two planes.
#[derive(Debug, Clone, PartialEq)]
pub enum PlanePlane {
    /// Planes coincide within tolerance.
    Coincident,
    /// Planes are parallel and distinct.
    Parallel,
    /// Planes intersect in a line.
    Line(Line3),
}

/// Intersect two planes in closed form.
///
/// The line direction is `n_a × n_b`; the line point is the standard
/// two-plane formula `((n_a·p_a)(n_b×d) + (n_b·p_b)(d×n_a)) / |d|²`.
pub fn plane_plane(a: &Plane, b: &Plane, tol: &Tol) -> PlanePlane {
    let dir = a.normal.cross(&b.normal);
    if dir.norm() <= tol.angular {
        return if a.signed_distance(&b.point).abs() <= tol.length {
            PlanePlane::Coincident
        } else {
            PlanePlane::Parallel
        };
    }
    let na = a.normal.into_inner();
    let nb = b.normal.into_inner();
    let da = na.dot(&a.point.coords);
    let db = nb.dot(&b.point.coords);
    let p = (nb.cross(&dir) * da + dir.cross(&na) * db) / dir.norm_squared();
    PlanePlane::Line(Line3::new(Point3::from(p), dir))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn offset_axis_planes_meet_on_expected_line() {
        let tol = Tol::default();
        // x = 1 and y = 2 → line {(1, 2, t)} along z
        let a = Plane::new(Point3::new(1.0_f64, 0.0_f64, 0.0_f64), Vector3::x());
        let b = Plane::new(Point3::new(0.0_f64, 2.0_f64, 0.0_f64), Vector3::y());
        match plane_plane(&a, &b, &tol) {
            PlanePlane::Line(line) => {
                assert!(line.dir.cross(&Vector3::z()).norm() < 1e-12_f64);
                assert!(a.signed_distance(&line.origin).abs() < 1e-12_f64);
                assert!(b.signed_distance(&line.origin).abs() < 1e-12_f64);
            }
            other => panic!("expected line, got {other:?}"),
        }
    }
}
