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
//!
//! # Cylindrical patch under an oblique rim
//!
//! When the patch is cut obliquely (one rim is a circle arc at the bottom, the
//! other an *ellipse* arc lying in the cut plane), the upper limit of the `z`
//! integral is no longer the constant `L` but a function `z₁(φ)` of the angle:
//! the axial coordinate at which the cut plane crosses the ruling line at angle
//! `φ`. Writing the cut plane as `n_p · (x − p₀) = 0` and substituting the
//! cylinder ruling `x(φ, z) = c + r(cosφ û + sinφ v̂) + z ẑ` gives
//!
//! ```text
//!   z₁(φ) = [ n_p · (p₀ − c) − r ( (n_p·û) cosφ + (n_p·v̂) sinφ ) ] / (n_p · ẑ).
//! ```
//!
//! With the lower rim at the cylinder bottom (`z = 0` in the `c`-anchored
//! frame), the patch integral becomes
//!
//! ```text
//!   ∫∫ x · n̂ dA = r ∫_{φ₀}^{φ₁} (c · n̂(φ) + r) · z₁(φ) dφ,
//! ```
//!
//! a product of two trigonometric polynomials of degree 1 in `(cosφ, sinφ)`,
//! integrated in closed form via the elementary primitives of `cosφ`, `sinφ`,
//! `cos²φ`, `sin²φ` and `cosφ sinφ` over `[φ₀, φ₁]` (see
//! [`oblique_patch_integral`]). The cut plane `(n_p, p₀)` is recovered from the
//! ellipse rim itself (its plane is the cut plane). When the rim is a circle
//! (`z₁(φ) ≡ L`) this reduces to the constant-height closed form above, so a
//! single routine serves both.

use crate::brep::Brep;
use crate::geom::{CurveGeom, SurfaceGeom};
use crate::math::{Point3, Vec3};
use crate::primitives::{Circle3, Ellipse3, Plane};
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
    /// A cylindrical face whose boundary the closed-form patch integral cannot
    /// interpret: it has neither the two circular rim arcs of a straight patch
    /// nor the circle-arc + ellipse-arc rims of an oblique cut (e.g. a cylinder
    /// face left in a state no builder produces). Reported rather than silently
    /// contributing zero (`DESIGN.md` §6-4: no silent zero).
    UnsupportedCylinderFace,
    /// A face with a surface kind the volume integral does not handle.
    UnsupportedSurface,
}

/// The signed volume of a [`Brep`], in cubic metres.
///
/// The surface is assumed closed and outward-oriented (as produced by the
/// extrusion builder); a positive result confirms that orientation. The
/// computation is exact closed-form for planar faces (polygonal or with circular
/// arcs / a chord — a circular segment) and for the cylinder patches the
/// extruder and the cut produce, including the *oblique* patch whose upper rim
/// is an ellipse (see the module docs).
///
/// This is the lenient entry point used as an orientation check: a face the
/// integral cannot evaluate in closed form (a cylinder face in a configuration
/// no builder produces, a surface kind not yet supported) contributes **zero**
/// rather than erroring. Use [`signed_volume_checked`] when an unsupported face
/// must be surfaced instead of silently skipped.
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
/// (a cylinder face in a state no builder produces, a surface kind not yet
/// handled) yields [`VolumeError`] instead of contributing zero.
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
                    Some(SurfaceGeom::Cylinder(cyl)) => {
                        f(cylinder_face_integral(brep, cyl, face.outer)
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
/// is a disk, `Area = πr²`; bounded by arcs of one ellipse (the oblique cut's
/// cap) it is an ellipse, `Area = πab`. Otherwise we take the polygon area
/// through the boundary vertices (fan triangulation, exact for any simple
/// polygon) and add a **circular-segment correction** for every arc edge — the
/// signed lens area `½r²(Δθ − sinΔθ)` between the arc and its chord — so a planar
/// face whose boundary mixes straight edges and arcs (e.g. a chord-cut
/// cylinder's circular-segment cap) integrates exactly.
fn planar_loop_integral(brep: &Brep, loop_id: crate::topo::arena::Id<Loop>, n_out: Vec3) -> f64 {
    if let Some((centre, radius)) = disk_loop(brep, loop_id) {
        let d = (centre - Point3::origin()).dot(n_out);
        // The disk area is **signed** by the loop's traversal about `n_out`: a CCW
        // round cap (outer boundary) encloses `+πr²`, a CW round hole encloses
        // `−πr²`. Read the sign from the total swept central angle of the arcs
        // (±2π) oriented to `n_out`, so a circular *hole* subtracts rather than
        // adds (the bug a fixed `+πr²` caused for a void cap).
        let total = disk_signed_sweep(brep, loop_id, n_out);
        let sign = if total >= 0.0 { 1.0 } else { -1.0 };
        return d * sign * std::f64::consts::PI * radius * radius;
    }
    if let Some((centre, a, b)) = ellipse_loop(brep, loop_id) {
        // Same as the disk case but with the ellipse area πab. The cap of an
        // oblique cylinder cut is a full ellipse (two half-ellipse arcs joined at
        // the seam); its planar integral is d·(±πab), signed by traversal about
        // `n_out` (read from the total swept ellipse-parameter, ±2π).
        let d = (centre - Point3::origin()).dot(n_out);
        let total = ellipse_signed_sweep(brep, loop_id, n_out);
        let sign = if total >= 0.0 { 1.0 } else { -1.0 };
        return d * sign * std::f64::consts::PI * a * b;
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

/// The total signed central angle swept by a disk loop's arcs, oriented to
/// `n_out` (positive ≈ +2π for a CCW round cap, negative ≈ −2π for a CW hole).
fn disk_signed_sweep(brep: &Brep, loop_id: crate::topo::arena::Id<Loop>, n_out: Vec3) -> f64 {
    let mut total = 0.0_f64;
    for &he_id in &lp_half_edges(brep, loop_id) {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        if let Some(CurveGeom::Circle(c)) = brep.geom.curve(he.curve) {
            let dtheta = he.boundary[1] - he.boundary[0];
            let about = c.normal().as_vec().dot(n_out);
            let sign = if about >= 0.0 { 1.0 } else { -1.0 };
            total += sign * dtheta;
        }
    }
    total
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

/// If `loop_id` is bounded entirely by arcs of one shared ellipse (an oblique
/// cut's cap), return its `(centre, semi_major, semi_minor)`; otherwise `None`.
fn ellipse_loop(brep: &Brep, loop_id: crate::topo::arena::Id<Loop>) -> Option<(Point3, f64, f64)> {
    let lp = brep.topo.loops.get(loop_id)?;
    if lp.half_edges.is_empty() {
        return None;
    }
    let mut ell = None;
    for &he_id in &lp.half_edges {
        let he = brep.topo.half_edges.get(he_id)?;
        match brep.geom.curve(he.curve)? {
            CurveGeom::Ellipse(e) => {
                if let Some(prev) = ell {
                    if prev != *e {
                        return None;
                    }
                } else {
                    ell = Some(*e);
                }
            }
            _ => return None,
        }
    }
    ell.map(|e| (e.center(), e.semi_major(), e.semi_minor()))
}

/// The total signed sweep of an ellipse loop's arcs, oriented to `n_out`
/// (positive ≈ +2π for a CCW cap, negative ≈ −2π for a hole). The ellipse angle
/// parameter `t` advances about the ellipse normal, so the sign of `n_out·normal`
/// orients it the same way [`disk_signed_sweep`] does for a circle.
fn ellipse_signed_sweep(brep: &Brep, loop_id: crate::topo::arena::Id<Loop>, n_out: Vec3) -> f64 {
    let mut total = 0.0_f64;
    for &he_id in &lp_half_edges(brep, loop_id) {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        if let Some(CurveGeom::Ellipse(e)) = brep.geom.curve(he.curve) {
            let dt = he.boundary[1] - he.boundary[0];
            let about = e.normal().as_vec().dot(n_out);
            let sign = if about >= 0.0 { 1.0 } else { -1.0 };
            total += sign * dt;
        }
    }
    total
}

/// `∫ x · n̂ dA` over a cylinder patch face, in closed form.
///
/// Handles both rim configurations the extruder and the cut produce:
///
/// * **Straight patch** — two circular rims (top and bottom), constant height
///   `L`. The closed form is `r·L·(c·n̂ integrated + r·Δφ)`.
/// * **Oblique patch** — one circular rim (the bottom, at the cylinder base) and
///   one *ellipse* rim (the cut). The ellipse's own plane is the cut plane, from
///   which the angle-dependent upper height `z₁(φ)` is recovered and the patch
///   integral [`oblique_patch_integral`] is evaluated.
///
/// Returns `None` for any other cylinder face (no recognisable rim pair), so the
/// caller can report it rather than silently treat it as zero.
fn cylinder_face_integral(
    brep: &Brep,
    cyl: &crate::primitives::Cylinder,
    loop_id: crate::topo::arena::Id<Loop>,
) -> Option<f64> {
    let lp = brep.topo.loops.get(loop_id)?;
    let mut circles: Vec<(f64, f64, Circle3)> = Vec::new();
    let mut ellipses: Vec<(f64, f64, Ellipse3)> = Vec::new();
    for &he_id in &lp.half_edges {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        match brep.geom.curve(he.curve) {
            Some(CurveGeom::Circle(circle)) => {
                circles.push((he.boundary[0], he.boundary[1], *circle));
            }
            Some(CurveGeom::Ellipse(ell)) => {
                ellipses.push((he.boundary[0], he.boundary[1], *ell));
            }
            _ => {}
        }
    }

    // The cylinder axis defines the patch's `ẑ` and the (û, v̂) the angle
    // parameter is written against (shared with `Circle3::point_at`).
    let axis = cyl.axis().dir();
    let axis_vec = axis.as_vec();
    let r = cyl.radius();
    let (u, v) = crate::primitives::plane_basis(axis);
    let proj = |p: Point3| (p - Point3::origin()).dot(axis_vec);

    match (circles.len(), ellipses.len()) {
        // ── Straight patch: two circular rims, constant height. ──────────────
        (2, 0) => {
            let (c0, c1) = (circles[0].2.center(), circles[1].2.center());
            let length = (c1 - c0).dot(axis_vec).abs();
            // Lower rim (smaller axial projection) is the bottom; use its angular
            // interval, oriented as stored on the outward-facing face.
            let (bottom_centre, phi0, phi1) = if proj(c0) <= proj(c1) {
                (c0, circles[0].0, circles[0].1)
            } else {
                (c1, circles[1].0, circles[1].1)
            };
            let c = bottom_centre - Point3::origin();
            let term_c = c.dot(u * (phi1.sin() - phi0.sin()) + v * (phi0.cos() - phi1.cos()));
            let term_r = r * (phi1 - phi0);
            Some(r * length * (term_c + term_r))
        }
        // ── Oblique patch: one circular bottom rim, one ellipse cut rim. ─────
        (1, 1) => {
            let (phi0, phi1, circle) = circles[0];
            let bottom_centre = circle.center();
            let ell = ellipses[0].2;
            // The ellipse rim lies in the cut plane; recover the plane from it.
            let cut_plane = Plane::new(ell.center(), ell.normal().as_vec()).ok()?;
            oblique_patch_integral(bottom_centre, u, v, axis_vec, r, phi0, phi1, &cut_plane)
        }
        _ => None,
    }
}

/// `∫ x · n̂ dA` over a cylinder patch with a circular bottom rim and an oblique
/// (cut-plane) top rim — the elliptical-rim case (module docs).
///
/// The patch is anchored at `bottom_centre` (the bottom circle's centre, taken
/// as `z = 0`), spans angles `[φ₀, φ₁]` about `(û, v̂)`, and rises to the cut
/// plane along the axis `ẑ`. With `c = bottom_centre − origin`, the integrand
/// `r (c·n̂(φ) + r) z₁(φ)` is a product of two degree-1 trig polynomials, where
/// `z₁(φ) = K − r(P_u cosφ + P_v sinφ)` is the axial height of the cut plane at
/// angle `φ`:
///
/// ```text
///   K   = n_p·(p₀ − bottom_centre) / (n_p·ẑ),
///   P_u = (n_p·û) / (n_p·ẑ),   P_v = (n_p·v̂) / (n_p·ẑ).
/// ```
///
/// Returns `None` if the cut plane is parallel to the axis (`n_p·ẑ ≈ 0`), which
/// is not an oblique cap (it would be an axis-parallel chord, handled by the
/// planar segment path instead).
#[allow(clippy::too_many_arguments)]
fn oblique_patch_integral(
    bottom_centre: Point3,
    u: Vec3,
    v: Vec3,
    axis: Vec3,
    r: f64,
    phi0: f64,
    phi1: f64,
    cut_plane: &Plane,
) -> Option<f64> {
    let n_p = cut_plane.normal().as_vec();
    let denom = n_p.dot(axis);
    if denom.abs() <= f64::EPSILON {
        return None;
    }
    let c = bottom_centre - Point3::origin();
    // g(φ) = c·n̂(φ) + r = A_g cosφ + B_g sinφ + r.
    let a_g = c.dot(u);
    let b_g = c.dot(v);
    // z₁(φ) = K − r(P_u cosφ + P_v sinφ).
    let k = n_p.dot(cut_plane.point() - bottom_centre) / denom;
    let p_u = n_p.dot(u) / denom;
    let p_v = n_p.dot(v) / denom;
    let a_z = -r * p_u;
    let b_z = -r * p_v;

    // Integrand / r = (A_g cosφ + B_g sinφ + r)(K + A_z cosφ + B_z sinφ).
    // Expand into the trig monomials and integrate each over [φ₀, φ₁].
    let const_term = r * k; // r·K
    let cos_term = a_g * k + r * a_z; // (A_g K + r A_z) cosφ
    let sin_term = b_g * k + r * b_z; // (B_g K + r B_z) sinφ
    let cos2_term = a_g * a_z; // A_g A_z cos²φ
    let sin2_term = b_g * b_z; // B_g B_z sin²φ
    let cossin_term = a_g * b_z + b_g * a_z; // (A_g B_z + B_g A_z) cosφ sinφ

    let integral = const_term * (phi1 - phi0)
        + cos_term * (phi1.sin() - phi0.sin())
        + sin_term * (phi0.cos() - phi1.cos())
        + cos2_term * int_cos2(phi0, phi1)
        + sin2_term * int_sin2(phi0, phi1)
        + cossin_term * int_cossin(phi0, phi1);

    Some(r * integral)
}

/// `∫_{a}^{b} cos²φ dφ = ½(φ + sinφ cosφ)|_a^b`.
fn int_cos2(a: f64, b: f64) -> f64 {
    0.5 * ((b - a) + (b.sin() * b.cos() - a.sin() * a.cos()))
}

/// `∫_{a}^{b} sin²φ dφ = ½(φ − sinφ cosφ)|_a^b`.
fn int_sin2(a: f64, b: f64) -> f64 {
    0.5 * ((b - a) - (b.sin() * b.cos() - a.sin() * a.cos()))
}

/// `∫_{a}^{b} cosφ sinφ dφ = ½ sin²φ|_a^b`.
fn int_cossin(a: f64, b: f64) -> f64 {
    0.5 * (b.sin() * b.sin() - a.sin() * a.sin())
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
