use crate::math::{Point3, Unit3};
use crate::primitives::{Line3, Plane};
use crate::tolerance::Tol;

/// Result of intersecting two planes.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
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
    let na = a.normal().as_vec();
    let nb = b.normal().as_vec();
    let dir = na.cross(nb);
    if dir.norm() <= tol.angular {
        return if a.signed_distance(b.point()).abs() <= tol.length {
            PlanePlane::Coincident
        } else {
            PlanePlane::Parallel
        };
    }
    let da = na.dot(a.point() - Point3::origin());
    let db = nb.dot(b.point() - Point3::origin());
    let p = (nb.cross(dir) * da + dir.cross(na) * db) * (1.0 / dir.norm_squared());
    let line = Line3::new_unchecked(
        Point3::new(p.x, p.y, p.z),
        Unit3::new_unchecked(dir * (1.0 / dir.norm())),
    );
    PlanePlane::Line(line)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn offset_axis_planes_meet_on_expected_line() {
        let tol = Tol::default();
        // x = 1 and y = 2 → line {(1, 2, t)} along z
        let a = Plane::new(Point3::new(1.0_f64, 0.0_f64, 0.0_f64), Vec3::X).expect("valid");
        let b = Plane::new(Point3::new(0.0_f64, 2.0_f64, 0.0_f64), Vec3::Y).expect("valid");
        match plane_plane(&a, &b, &tol) {
            PlanePlane::Line(line) => {
                assert!(line.dir().cross(Vec3::Z).norm() < 1e-12_f64);
                assert!(a.signed_distance(line.origin()).abs() < 1e-12_f64);
                assert!(b.signed_distance(line.origin()).abs() < 1e-12_f64);
            }
            other => panic!("expected line, got {other:?}"),
        }
    }
}
