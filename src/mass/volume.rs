//! Signed volume by the divergence theorem.
//!
//! For a closed, outward-oriented surface `∂Ω`, the divergence theorem with the
//! vector field `F = x` (so `∇·F = 3`) gives
//!
//! ```text
//!   V(Ω) = (1/3) ∮∮_{∂Ω} x · n̂ dA.
//! ```
//!
//! Each face contributes its own surface integral, and the sign of the result
//! is the orientation test used in Phase 2: a solid built with outward-facing
//! normals has `V > 0`. The contribution is evaluated in closed form per
//! surface kind.
//!
//! # Planar face
//!
//! On a planar face `n̂` is constant and every boundary point satisfies
//! `x · n̂ = d` (the plane's signed origin distance), so
//! `∫ x · n̂ dA = d · Area`. We evaluate it without needing `d` or `n̂`
//! explicitly by fan-triangulating each boundary loop from its first vertex and
//! summing the signed tetra-to-origin volumes
//! `(1/6) Σ q0 · (qi × qi₊₁)`; this telescopes to the exact polygon integral for
//! any simple polygon (convex or not — the H-section is concave). Outer loops
//! add, inner (hole) loops subtract.
//!
//! # Cylindrical (half-)patch
//!
//! Parameterise a patch of a cylinder of radius `r`, bottom-centre `c` and axis
//! `ẑ` by angle `φ ∈ [φ₀, φ₁]` and height `z ∈ [0, L]`:
//!
//! ```text
//!   x(φ, z) = c + r (cosφ û + sinφ v̂) + z ẑ,     n̂(φ) = cosφ û + sinφ v̂,
//!   dA = r dφ dz.
//! ```
//!
//! Because `ẑ · n̂ = 0` and `n̂ · n̂ = 1`,
//!
//! ```text
//!   ∫∫ x · n̂ dA = r L ∫_{φ₀}^{φ₁} (c · n̂ + r) dφ
//!               = r L [ c · ( û (sinφ₁ − sinφ₀) + v̂ (cosφ₀ − cosφ₁) )
//!                       + r (φ₁ − φ₀) ].
//! ```
//!
//! The face contribution is `(1/3)` of that. Summing the two half-patches
//! (`φ₀..π` and `π..2π`) recovers the full-cylinder term `(2/3)πr²L`, which with
//! the top cap `(1/3)L·πr²` totals `πr²L`, the exact cylinder volume.

use crate::brep::Brep;
use crate::geom::{CurveGeom, SurfaceGeom};
use crate::math::{Point3, Vec3};
use crate::primitives::Plane;
use crate::topo::{Loop, Sense};

/// The outward unit normal of a planar face: the plane normal, flipped when the
/// face sense is [`Sense::Reversed`].
fn outward_normal(plane: &Plane, sense: Sense) -> Vec3 {
    let n = plane.normal().as_vec();
    match sense {
        Sense::Same => n,
        Sense::Reversed => -n,
    }
}

/// A face configuration the closed-form volume integral does not yet cover.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum VolumeError {
    /// A cylindrical face whose boundary is not the two circular rim arcs the
    /// closed-form patch integral expects (e.g. an oblique cut leaving an
    /// elliptical rim). The elliptical-rim integral is a Phase 5 item; until then
    /// such a face cannot be integrated and is reported rather than silently
    /// contributing zero.
    UnsupportedCylinderFace,
    /// A face with a surface kind the volume integral does not handle.
    UnsupportedSurface,
}

/// The signed volume of a [`Brep`], in cubic metres.
///
/// The surface is assumed closed and outward-oriented (as produced by the
/// extrusion builder); a positive result confirms that orientation. The
/// computation is exact closed-form for planar faces (polygonal or with circular
/// arcs / a chord — a circular segment) and for the half-cylinder patches the
/// extruder and the cut produce.
///
/// This is the lenient entry point used as an orientation check: a face the
/// integral cannot evaluate in closed form (an oblique elliptical cylinder rim,
/// a surface kind not yet supported) contributes **zero** rather than erroring.
/// Use [`signed_volume_checked`] when an unsupported face must be surfaced
/// instead of silently skipped.
pub fn signed_volume(brep: &Brep) -> f64 {
    let mut vol = 0.0_f64;
    for_each_face_contribution(brep, |c| {
        vol += c.unwrap_or(0.0);
    });
    vol / 3.0_f64
}

/// The signed volume of a [`Brep`], failing if any face cannot be integrated in
/// closed form.
///
/// Identical to [`signed_volume`] except that an unsupported face configuration
/// (an oblique elliptical cylinder rim, a surface kind not yet handled) yields
/// [`VolumeError`] instead of contributing zero.
///
/// # Errors
///
/// Returns [`VolumeError`] for the first unsupported face encountered.
pub fn signed_volume_checked(brep: &Brep) -> Result<f64, VolumeError> {
    let mut vol = 0.0_f64;
    let mut err: Option<VolumeError> = None;
    for_each_face_contribution(brep, |c| match c {
        Ok(v) => vol += v,
        Err(e) => {
            if err.is_none() {
                err = Some(e);
            }
        }
    });
    match err {
        Some(e) => Err(e),
        None => Ok(vol / 3.0_f64),
    }
}

/// Drive `f` with each face's `∫ x · n̂ dA` contribution (or the reason it could
/// not be computed), so the lenient and checked entry points share one walk.
fn for_each_face_contribution(brep: &Brep, mut f: impl FnMut(Result<f64, VolumeError>)) {
    for solid_id in &brep.solids {
        let Some(solid) = brep.topo.solids.get(*solid_id) else {
            continue;
        };
        for shell_id in &solid.shells {
            let Some(shell) = brep.topo.shells.get(*shell_id) else {
                continue;
            };
            for face_id in &shell.faces {
                let Some(face) = brep.topo.faces.get(*face_id) else {
                    continue;
                };
                match brep.geom.surface(face.surface) {
                    Some(SurfaceGeom::Plane(plane)) => {
                        let n_out = outward_normal(plane, face.sense);
                        let mut v = planar_loop_integral(brep, face.outer, n_out);
                        for inner in &face.inners {
                            // Inner (hole) loops are wound *opposite* to the outer
                            // boundary in a valid B-rep (outer CCW ⇒ hole CW for
                            // the same `n_out`), so their `n_out`-oriented planar
                            // integral is already negative and is *added* — summing
                            // every loop with its natural sign yields
                            // `outer − hole`. (Subtracting would instead require
                            // CCW holes, which break sibling pairing.)
                            v += planar_loop_integral(brep, *inner, n_out);
                        }
                        f(Ok(v));
                    }
                    Some(SurfaceGeom::Cylinder(_)) => {
                        f(cylinder_face_integral(brep, face.outer)
                            .ok_or(VolumeError::UnsupportedCylinderFace));
                    }
                    None => f(Err(VolumeError::UnsupportedSurface)),
                }
            }
        }
    }
}

/// `∫ x · n̂ dA` over a planar boundary loop (`= 3·V` contribution).
///
/// On a planar face `x · n̂ = d` is constant (the plane's signed origin distance),
/// so the integral is `d · Area`, with `Area` the loop's signed planar area
/// oriented to `n̂`. A loop bounded entirely by arcs of one circle (a round cap)
/// is a disk, `Area = πr²`. Otherwise we take the polygon area through the
/// boundary vertices (fan triangulation, exact for any simple polygon) and add a
/// **circular-segment correction** for every arc edge — the signed lens area
/// `½r²(Δθ − sinΔθ)` between the arc and its chord — so a planar face whose
/// boundary mixes straight edges and arcs (e.g. a chord-cut cylinder's
/// circular-segment cap) integrates exactly.
fn planar_loop_integral(brep: &Brep, loop_id: crate::topo::arena::Id<Loop>, n_out: Vec3) -> f64 {
    if let Some((centre, radius)) = disk_loop(brep, loop_id) {
        let d = (centre - Point3::origin()).dot(n_out);
        return d * std::f64::consts::PI * radius * radius;
    }
    let verts = loop_vertices(brep, loop_id);
    // Polygon part: Σ q0 · (qi × qi₊₁) / 2 = d · PolyArea (the planar integral
    // over the chord polygon). A loop with fewer than three vertices encloses no
    // polygon area, but its arc edges may still bound a circular segment (a
    // half-disk: one arc + one chord, only two vertices), so we still run the
    // arc-correction pass below rather than returning early.
    let mut integral = 0.0_f64;
    if verts.len() >= 3 {
        let q0 = verts[0];
        let q0v = q0 - Point3::origin();
        let mut acc = 0.0_f64;
        for i in 1..verts.len() - 1 {
            let a = verts[i] - q0;
            let b = verts[i + 1] - q0;
            acc += q0v.dot(a.cross(b));
        }
        integral = acc / 2.0_f64;
    }

    // Arc corrections: each arc edge bulges off its chord by a circular segment.
    // The segment's signed area in the n̂-oriented plane is ½r²(Δθ − sinΔθ) with
    // Δθ the arc's signed central angle measured about +n̂ (the arc's own
    // boundary params run about the circle normal). The planar integral gains
    // `d · segment_area`, with d the constant plane offset `x · n̂`.
    for &he_id in &lp_half_edges(brep, loop_id) {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        if let Some(CurveGeom::Circle(c)) = brep.geom.curve(he.curve) {
            let r = c.radius();
            // Signed central angle of the arc as traversed (boundary[0]→[1]),
            // about the circle's own normal.
            let dtheta = he.boundary[1] - he.boundary[0];
            // Orient about n̂: the circle normal may agree or oppose n_out.
            let about = c.normal().as_vec().dot(n_out);
            let sign = if about >= 0.0 { 1.0 } else { -1.0 };
            let seg = 0.5 * r * r * (dtheta - dtheta.sin());
            let d = (c.center() - Point3::origin()).dot(n_out);
            integral += d * sign * seg;
        }
    }
    integral
}

/// The half-edge ids of a loop (helper for the arc-correction pass).
fn lp_half_edges(
    brep: &Brep,
    loop_id: crate::topo::arena::Id<Loop>,
) -> Vec<crate::topo::arena::Id<crate::topo::HalfEdge>> {
    brep.topo
        .loops
        .get(loop_id)
        .map(|lp| lp.half_edges.clone())
        .unwrap_or_default()
}

/// If `loop_id` is bounded entirely by circular arcs of one shared circle (a
/// round cap), return its `(centre, radius)`; otherwise `None`.
fn disk_loop(brep: &Brep, loop_id: crate::topo::arena::Id<Loop>) -> Option<(Point3, f64)> {
    let lp = brep.topo.loops.get(loop_id)?;
    if lp.half_edges.is_empty() {
        return None;
    }
    let mut circle = None;
    for &he_id in &lp.half_edges {
        let he = brep.topo.half_edges.get(he_id)?;
        match brep.geom.curve(he.curve)? {
            CurveGeom::Circle(c) => {
                if let Some(prev) = circle {
                    if prev != *c {
                        return None;
                    }
                } else {
                    circle = Some(*c);
                }
            }
            _ => return None,
        }
    }
    circle.map(|c| (c.center(), c.radius()))
}

/// `∫ x · n̂ dA` over a cylinder (half-)patch face, in closed form.
///
/// Returns `None` for a cylinder face whose boundary is not the two circular rim
/// arcs this closed form expects (e.g. an oblique cut's elliptical rim), so the
/// caller can report it rather than silently treat it as zero.
fn cylinder_face_integral(brep: &Brep, loop_id: crate::topo::arena::Id<Loop>) -> Option<f64> {
    let lp = brep.topo.loops.get(loop_id)?;
    // Find the two circular-arc half-edges (the top and bottom rims) and the
    // height L (from the vertical extent between the two circle centres). Any
    // ellipse rim means an oblique cut this closed form cannot integrate.
    let mut arcs: Vec<(f64, f64, crate::primitives::Circle3)> = Vec::new();
    for &he_id in &lp.half_edges {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        match brep.geom.curve(he.curve) {
            Some(CurveGeom::Circle(circle)) => {
                arcs.push((he.boundary[0], he.boundary[1], *circle));
            }
            Some(CurveGeom::Ellipse(_)) => return None,
            _ => {}
        }
    }
    if arcs.len() != 2 {
        // Not a half-cylinder face as built by the extruder / cut: unsupported.
        return None;
    }
    // The two rims share the same radius and axis; their centres differ by L·ẑ.
    let (c0, c1) = (arcs[0].2.center(), arcs[1].2.center());
    let axis = arcs[0].2.normal().as_vec();
    let length = (c1 - c0).dot(axis).abs();
    let r = arcs[0].2.radius();

    // Use the lower rim (smaller projection onto the axis) as the bottom centre
    // and its angular interval, oriented as stored on the outward-facing face.
    let proj = |p: Point3| (p - Point3::origin()).dot(axis);
    let (bottom_centre, phi0, phi1) = if proj(c0) <= proj(c1) {
        (c0, arcs[0].0, arcs[0].1)
    } else {
        (c1, arcs[1].0, arcs[1].1)
    };

    // Reconstruct the same orthonormal basis the circle uses, so that û, v̂ here
    // match the angle parameter on the boundary.
    let (u, v) = crate::primitives::plane_basis(arcs[0].2.normal());
    let c = bottom_centre - Point3::origin();
    let term_c = c.dot(u * (phi1.sin() - phi0.sin()) + v * (phi0.cos() - phi1.cos()));
    let term_r = r * (phi1 - phi0);
    Some(r * length * (term_c + term_r))
}

/// Collect a loop's vertices (the start vertex of each half-edge) as points.
fn loop_vertices(brep: &Brep, loop_id: crate::topo::arena::Id<Loop>) -> Vec<Point3> {
    let Some(lp) = brep.topo.loops.get(loop_id) else {
        return Vec::new();
    };
    let mut pts = Vec::with_capacity(lp.half_edges.len());
    for &he_id in &lp.half_edges {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        let Some(vert) = brep.topo.vertices.get(he.start) else {
            continue;
        };
        if let Some(p) = brep.geom.point(vert.point).and_then(|g| g.as_point()) {
            pts.push(p);
        }
    }
    pts
}
