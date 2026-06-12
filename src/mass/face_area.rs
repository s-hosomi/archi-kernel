//! Per-face area in closed form (the shared core of surface area and formwork).
//!
//! A planar face's area is the net of its boundary loops oriented to the face's
//! outward normal: the outer loop contributes its positive planar area and each
//! hole loop its negative area, so the net is `outer − holes` (`DESIGN.md` §6-4;
//! openings show up as inner loops, so their area is removed automatically). A
//! cylinder patch's area is `r · Δφ · L` for a straight patch, or the oblique
//! arc-length integral for a cut patch.

use crate::brep::Brep;
use crate::geom::{CurveGeom, SurfaceGeom};
use crate::math::{Point3, Vec3};
use crate::primitives::{plane_basis, Plane};
use crate::topo::arena::Id;
use crate::topo::{Face, Loop, Sense};

/// A face whose area / orientation could not be computed in closed form.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum AreaError {
    /// A cylinder face whose boundary is not a recognised rim pair (mirrors the
    /// volume integral's limitation).
    UnsupportedCylinderFace,
    /// A surface kind the area computation does not handle.
    UnsupportedSurface,
}

/// The outward unit normal of a planar face (plane normal, flipped on
/// [`Sense::Reversed`]).
pub(crate) fn planar_outward_normal(plane: &Plane, sense: Sense) -> Vec3 {
    let n = plane.normal().as_vec();
    match sense {
        Sense::Same => n,
        Sense::Reversed => -n,
    }
}

/// The (positive) area of a face, in square metres, with holes removed.
///
/// # Errors
///
/// [`AreaError`] for an unsupported cylinder face or surface kind.
pub(crate) fn face_area(brep: &Brep, face: &Face) -> Result<f64, AreaError> {
    match brep.geom.surface(face.surface) {
        Some(SurfaceGeom::Plane(plane)) => {
            let n_out = planar_outward_normal(plane, face.sense);
            let mut area = planar_loop_area(brep, face.outer, n_out);
            for &inner in &face.inners {
                // Hole loops are wound opposite the outer loop, so their oriented
                // area is negative; adding removes them.
                area += planar_loop_area(brep, inner, n_out);
            }
            Ok(area.abs())
        }
        Some(SurfaceGeom::Cylinder(cyl)) => {
            cylinder_face_area(brep, cyl, face.outer).ok_or(AreaError::UnsupportedCylinderFace)
        }
        None => Err(AreaError::UnsupportedSurface),
    }
}

/// Signed planar area of a boundary loop oriented to `n_out` (the same arc-aware
/// area used by the volume integral, without the `d` plane-offset factor).
fn planar_loop_area(brep: &Brep, loop_id: Id<Loop>, n_out: Vec3) -> f64 {
    let Some(lp) = brep.topo.loops.get(loop_id) else {
        return 0.0;
    };
    // Polygon (chord) part via the shoelace sum projected onto the n_out frame.
    let verts: Vec<Point3> = lp
        .half_edges
        .iter()
        .filter_map(|&he_id| {
            let he = brep.topo.half_edges.get(he_id)?;
            let v = brep.topo.vertices.get(he.start)?;
            brep.geom.point(v.point).and_then(|g| g.as_point())
        })
        .collect();
    let mut area = 0.0_f64;
    if verts.len() >= 3 {
        // 2·Area·n̂ = Σ qᵢ × qᵢ₊₁ ; project onto n_out for the oriented area.
        let mut acc = Vec3::ZERO;
        let o = Point3::origin();
        for i in 0..verts.len() {
            let a = verts[i] - o;
            let b = verts[(i + 1) % verts.len()] - o;
            acc = acc + a.cross(b);
        }
        area = 0.5 * acc.dot(n_out);
    }
    // Arc corrections: circular and elliptical bulges off their chords.
    for &he_id in &lp.half_edges {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        match brep.geom.curve(he.curve) {
            Some(CurveGeom::Circle(c)) => {
                let r = c.radius();
                let dtheta = he.boundary[1] - he.boundary[0];
                let about = c.normal().as_vec().dot(n_out);
                let sign = if about >= 0.0 { 1.0 } else { -1.0 };
                area += sign * 0.5 * r * r * (dtheta - dtheta.sin());
            }
            Some(CurveGeom::Ellipse(e)) => {
                // An ellipse arc's segment area off its chord: in the ellipse's
                // own (a·cos t, b·sin t) frame the swept-minus-triangle area is
                // ½ab(Δt − sinΔt), the affine image of the circular segment.
                let a = e.semi_major();
                let b = e.semi_minor();
                let dt = he.boundary[1] - he.boundary[0];
                let about = e.normal().as_vec().dot(n_out);
                let sign = if about >= 0.0 { 1.0 } else { -1.0 };
                area += sign * 0.5 * a * b * (dt - dt.sin());
            }
            _ => {}
        }
    }
    area
}

/// Lateral area of a cylinder patch face: `r·Δφ·L` for a straight patch, the
/// oblique arc-length integral for a cut patch. Returns `None` for an
/// unrecognised rim pair.
fn cylinder_face_area(
    brep: &Brep,
    cyl: &crate::primitives::Cylinder,
    loop_id: Id<Loop>,
) -> Option<f64> {
    let lp = brep.topo.loops.get(loop_id)?;
    let mut circles = Vec::new();
    let mut ellipses = Vec::new();
    for &he_id in &lp.half_edges {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        match brep.geom.curve(he.curve) {
            Some(CurveGeom::Circle(c)) => circles.push((he.boundary[0], he.boundary[1], *c)),
            Some(CurveGeom::Ellipse(e)) => ellipses.push((he.boundary[0], he.boundary[1], *e)),
            _ => {}
        }
    }
    let axis = cyl.axis().dir();
    let axis_vec = axis.as_vec();
    let r = cyl.radius();
    match (circles.len(), ellipses.len()) {
        (2, 0) => {
            // Straight patch: constant height L, angular span Δφ ⇒ area r·Δφ·L.
            let (c0, c1) = (circles[0].2.center(), circles[1].2.center());
            let length = (c1 - c0).dot(axis_vec).abs();
            let dphi = (circles[0].1 - circles[0].0).abs();
            Some(r * dphi * length)
        }
        (1, 1) => {
            // Oblique patch: area = ∫ r·z₁(φ) dφ over the angular span, with
            // z₁(φ) the axial height of the cut plane (module docs of `volume`).
            let (phi0, phi1, circle) = circles[0];
            let bottom_centre = circle.center();
            let ell = ellipses[0].2;
            let cut_plane = Plane::new(ell.center(), ell.normal().as_vec()).ok()?;
            let (u, v) = plane_basis(axis);
            let n_p = cut_plane.normal().as_vec();
            let denom = n_p.dot(axis_vec);
            if denom.abs() <= f64::EPSILON {
                return None;
            }
            let k = n_p.dot(cut_plane.point() - bottom_centre) / denom;
            let p_u = n_p.dot(u) / denom;
            let p_v = n_p.dot(v) / denom;
            // ∫ r·(K − r(P_u cosφ + P_v sinφ)) dφ over [φ0, φ1].
            let integral = k * (phi1 - phi0)
                - r * p_u * (phi1.sin() - phi0.sin())
                - r * p_v * (phi0.cos() - phi1.cos());
            Some((r * integral).abs())
        }
        _ => None,
    }
}
