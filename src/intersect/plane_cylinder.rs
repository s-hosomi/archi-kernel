use crate::math::Unit3;

use super::line_plane;
use crate::primitives::{Circle3, Cylinder, Ellipse3, Line3, Plane};
use crate::tolerance::Tol;

/// Result of intersecting a plane with an infinite circular cylinder.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PlaneCylinder {
    /// Plane parallel to the axis, farther than the radius.
    None,
    /// Plane parallel to the axis, tangent to the surface.
    TangentLine(Line3),
    /// Plane parallel to the axis, cutting through: two ruling lines.
    TwoLines([Line3; 2]),
    /// Plane perpendicular to the axis: a circle of the cylinder radius.
    Circle(Circle3),
    /// Oblique plane: an ellipse with semi-minor = radius and
    /// semi-major = radius / |cos θ| where θ is the angle between the
    /// plane normal and the axis.
    Ellipse(Ellipse3),
}

/// Intersect a plane with a cylinder in closed form. Covers all five cases.
pub fn plane_cylinder(plane: &Plane, cyl: &Cylinder, tol: &Tol) -> PlaneCylinder {
    // Cosine of the angle between the plane normal and the cylinder axis.
    let cos = plane.normal().dot(cyl.axis().dir().as_vec());

    if cos.abs() <= tol.angular {
        // Plane is parallel to the axis.
        let sd = plane.signed_distance(cyl.axis().origin());
        let foot = plane.project_point(cyl.axis().origin());
        if tol.eq_length(sd.abs(), cyl.radius()) {
            return PlaneCylinder::TangentLine(Line3::new_unchecked(foot, cyl.axis().dir()));
        }
        if sd.abs() > cyl.radius() {
            return PlaneCylinder::None;
        }
        // Two ruling lines, offset from the axis projection by
        // w = sqrt(r² − d²) along u = axis × normal (unit, in the plane).
        let w = (cyl.radius() * cyl.radius() - sd * sd).sqrt();
        let u_raw = cyl.axis().dir().cross(plane.normal().as_vec());
        // u_raw is guaranteed non-zero because axis and plane normal are not
        // parallel in this branch.
        let u = Unit3::new_unchecked(u_raw * (1.0 / u_raw.norm()));
        return PlaneCylinder::TwoLines([
            Line3::new_unchecked(foot + u.as_vec() * w, cyl.axis().dir()),
            Line3::new_unchecked(foot - u.as_vec() * w, cyl.axis().dir()),
        ]);
    }

    // Axis crosses the plane: the section center is axis ∩ plane.
    let axis = cyl.axis();
    let center =
        line_plane(&axis, plane, tol).expect("axis is not parallel to the plane in this branch");

    if (cos.abs() - 1.0).abs() <= tol.angular {
        return PlaneCylinder::Circle(Circle3::new_unchecked(
            center,
            cyl.axis().dir(),
            cyl.radius(),
        ));
    }

    // Oblique section: ellipse. The major axis is the projection of the
    // cylinder axis direction onto the plane.
    let proj = cyl.axis().dir().as_vec() - plane.normal().as_vec() * cos;
    let major_dir = Unit3::new_unchecked(proj * (1.0 / proj.norm()));
    PlaneCylinder::Ellipse(Ellipse3::new_unchecked(
        center,
        plane.normal(),
        major_dir,
        cyl.radius() / cos.abs(),
        cyl.radius(),
    ))
}
