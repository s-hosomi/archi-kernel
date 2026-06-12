//! Verification against closed-form analytic solutions.
//!
//! Every intersection result is checked against geometry that can be solved
//! by hand. Tolerances are explicit in each assertion.

use archi_kernel::intersect::{plane_cylinder, plane_plane, PlaneCylinder, PlanePlane};
use archi_kernel::{Cylinder, KernelError, Line3, Plane, Point3, Tol, Vec3};

const EPS: f64 = 1e-12;

#[test]
fn orthogonal_planes_meet_in_the_shared_axis() {
    let tol = Tol::default();
    // x = 0 and y = 0 → the z axis.
    let a = Plane::new(Point3::origin(), Vec3::X).expect("valid plane");
    let b = Plane::new(Point3::origin(), Vec3::Y).expect("valid plane");
    match plane_plane(&a, &b, &tol) {
        PlanePlane::Line(line) => {
            assert!(line.dir().cross(Vec3::Z).norm() < EPS);
            assert!(line.distance_to_point(Point3::origin()) < EPS);
        }
        other => panic!("expected line, got {other:?}"),
    }
}

#[test]
fn parallel_and_coincident_planes_are_distinguished() {
    let tol = Tol::default();
    let a = Plane::new(Point3::origin(), Vec3::Z).expect("valid plane");
    let b = Plane::new(Point3::new(0.0_f64, 0.0_f64, 0.5_f64), Vec3::Z).expect("valid plane");
    let c = Plane::new(Point3::new(7.0_f64, -3.0_f64, 0.0_f64), Vec3::Z * 2.0_f64)
        .expect("valid plane");
    assert_eq!(plane_plane(&a, &b, &tol), PlanePlane::Parallel);
    assert_eq!(plane_plane(&a, &c, &tol), PlanePlane::Coincident);
}

#[test]
fn perpendicular_section_of_column_is_a_circle() {
    let tol = Tol::default();
    // Round column, r = 0.3 m, axis = z. Slab plane at z = 2.8 m.
    let cyl = Cylinder::new(
        Line3::new(Point3::origin(), Vec3::Z).expect("valid line"),
        0.3_f64,
    )
    .expect("valid cylinder");
    let plane = Plane::new(Point3::new(0.0_f64, 0.0_f64, 2.8_f64), Vec3::Z).expect("valid plane");
    match plane_cylinder(&plane, &cyl, &tol) {
        PlaneCylinder::Circle(circle) => {
            assert!((circle.radius() - 0.3_f64).abs() < EPS);
            assert!(
                (circle.center() - Point3::new(0.0_f64, 0.0_f64, 2.8_f64)).norm() < EPS,
                "center = {:?}",
                circle.center()
            );
        }
        other => panic!("expected circle, got {other:?}"),
    }
}

#[test]
fn oblique_section_is_an_ellipse_with_known_axes() {
    let tol = Tol::default();
    // Plane normal at 45° to the axis: semi-major = r·√2, semi-minor = r.
    let r = 0.3_f64;
    let cyl = Cylinder::new(
        Line3::new(Point3::origin(), Vec3::Z).expect("valid line"),
        r,
    )
    .expect("valid cylinder");
    let plane =
        Plane::new(Point3::origin(), Vec3::new(0.0_f64, 1.0_f64, 1.0_f64)).expect("valid plane");
    match plane_cylinder(&plane, &cyl, &tol) {
        PlaneCylinder::Ellipse(e) => {
            assert!((e.semi_major() - r * 2.0_f64.sqrt()).abs() < EPS);
            assert!((e.semi_minor() - r).abs() < EPS);
            assert!((e.center() - Point3::origin()).norm() < EPS);
            // Every sampled point must lie on both surfaces.
            for i in 0..16 {
                let t = f64::from(i) * std::f64::consts::TAU / 16.0_f64;
                let p = e.point_at(t);
                assert!(plane.signed_distance(p).abs() < EPS, "off plane at t={t}");
                assert!(cyl.signed_distance(p).abs() < EPS, "off cylinder at t={t}");
            }
        }
        other => panic!("expected ellipse, got {other:?}"),
    }
}

#[test]
fn axis_parallel_plane_cuts_two_ruling_lines() {
    let tol = Tol::default();
    // Plane x = 0.1 through a column of r = 0.3 about the z axis:
    // lines at (0.1, ±√(0.09 − 0.01), t) = (0.1, ±√0.08, t).
    let r = 0.3_f64;
    let d = 0.1_f64;
    let cyl = Cylinder::new(
        Line3::new(Point3::origin(), Vec3::Z).expect("valid line"),
        r,
    )
    .expect("valid cylinder");
    let plane = Plane::new(Point3::new(d, 0.0_f64, 0.0_f64), Vec3::X).expect("valid plane");
    match plane_cylinder(&plane, &cyl, &tol) {
        PlaneCylinder::TwoLines(lines) => {
            let w = (r * r - d * d).sqrt();
            for line in &lines {
                assert!(line.dir().cross(Vec3::Z).norm() < EPS);
                assert!(plane.signed_distance(line.origin()).abs() < EPS);
                assert!(cyl.signed_distance(line.origin()).abs() < EPS);
            }
            let gap = (lines[0].origin() - lines[1].origin()).norm();
            assert!((gap - 2.0_f64 * w).abs() < EPS, "gap = {gap}");
        }
        other => panic!("expected two lines, got {other:?}"),
    }
}

#[test]
fn tangent_and_missing_planes_are_detected() {
    let tol = Tol::default();
    let r = 0.3_f64;
    let cyl = Cylinder::new(
        Line3::new(Point3::origin(), Vec3::Z).expect("valid line"),
        r,
    )
    .expect("valid cylinder");

    let tangent = Plane::new(Point3::new(r, 0.0_f64, 0.0_f64), Vec3::X).expect("valid plane");
    match plane_cylinder(&tangent, &cyl, &tol) {
        PlaneCylinder::TangentLine(line) => {
            assert!(cyl.signed_distance(line.origin()).abs() < EPS);
        }
        other => panic!("expected tangent line, got {other:?}"),
    }

    let missing =
        Plane::new(Point3::new(r + 0.5_f64, 0.0_f64, 0.0_f64), Vec3::X).expect("valid plane");
    assert_eq!(plane_cylinder(&missing, &cyl, &tol), PlaneCylinder::None);
}

// ── Constructor error tests ──────────────────────────────────────────────────

#[test]
fn plane_rejects_zero_normal() {
    assert_eq!(
        Plane::new(Point3::origin(), Vec3::ZERO),
        Err(KernelError::ZeroNormal)
    );
}

#[test]
fn line3_rejects_zero_direction() {
    assert_eq!(
        Line3::new(Point3::origin(), Vec3::ZERO),
        Err(KernelError::ZeroDirection)
    );
}

#[test]
fn cylinder_rejects_non_positive_radius() {
    let axis = Line3::new(Point3::origin(), Vec3::Z).expect("valid line");
    assert!(matches!(
        Cylinder::new(axis, 0.0_f64),
        Err(KernelError::NonPositiveRadius { .. })
    ));
    let axis2 = Line3::new(Point3::origin(), Vec3::Z).expect("valid line");
    assert!(matches!(
        Cylinder::new(axis2, -0.5_f64),
        Err(KernelError::NonPositiveRadius { .. })
    ));
}

// ── serde smoke tests ────────────────────────────────────────────────────────

#[cfg(feature = "serde")]
mod serde_tests {
    use archi_kernel::{Point3, Unit3, Vec3};

    #[test]
    fn point3_roundtrip() {
        let p = Point3::new(1.0_f64, 2.0_f64, 3.0_f64);
        let json = serde_json::to_string(&p).expect("serialize");
        let q: Point3 = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(p, q);
    }

    #[test]
    fn vec3_roundtrip() {
        let v = Vec3::new(4.0_f64, 5.0_f64, 6.0_f64);
        let json = serde_json::to_string(&v).expect("serialize");
        let w: Vec3 = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(v, w);
    }

    #[test]
    fn unit3_roundtrip() {
        let u = Unit3::Z;
        let json = serde_json::to_string(&u).expect("serialize");
        let v: Unit3 = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(u, v);
    }
}
