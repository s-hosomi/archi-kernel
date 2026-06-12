//! Higher mass properties built on the per-face area core: surface area,
//! centroid, and the formwork-area take-off (`DESIGN.md` §6-4).

use crate::brep::Brep;
use crate::geom::{CurveGeom, SurfaceGeom};
use crate::math::{Point3, Vec3};
use crate::primitives::{plane_basis, Plane};
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::{Face, Loop};

use super::face_area::{face_area, planar_outward_normal, AreaError};

/// Total surface area of a B-rep, in square metres.
///
/// Sums every face's net area (outer loop minus holes), exact in closed form for
/// planar (polygon + arc segments) and cylinder faces.
///
/// # Errors
///
/// [`AreaError`] for the first face whose area cannot be computed in closed form.
pub fn surface_area(brep: &Brep) -> Result<f64, AreaError> {
    let mut total = 0.0_f64;
    for_each_face(brep, |face| {
        total += face_area(brep, face)?;
        Ok(())
    })?;
    Ok(total)
}

/// A face configuration the centroid integral does not handle.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum CentroidError {
    /// The solid is empty (zero volume), so a centroid is undefined.
    EmptySolid,
    /// A non-planar face was present; the first-moment closed form here covers
    /// planar-faced solids exactly. (Cylinder first moments are a later item;
    /// reported rather than silently skipped, `DESIGN.md` §6-4.)
    UnsupportedCurvedFace,
}

/// Centroid of a B-rep, in metres (`DESIGN.md` §6-4).
///
/// Computed exactly for a planar-faced solid by signed-tetrahedron decomposition
/// from the origin: each planar face is fan-triangulated, every triangle forms a
/// tetrahedron with the origin whose signed volume and barycentre are summed,
/// and the centroid is the volume-weighted mean `Σ vᵢ cᵢ / Σ vᵢ`. This is the
/// divergence-theorem first moment specialised to planar faces, and is exact
/// regardless of face convexity (the H-section is concave).
///
/// # Errors
///
/// * [`CentroidError::UnsupportedCurvedFace`] if any face is non-planar.
/// * [`CentroidError::EmptySolid`] if the total volume is ~zero.
pub fn centroid(brep: &Brep, tol: &Tol) -> Result<Point3, CentroidError> {
    let mut vol = 0.0_f64;
    let mut moment = Vec3::ZERO;
    let mut err: Option<CentroidError> = None;
    for_each_face(brep, |face| {
        let Some(SurfaceGeom::Plane(_)) = brep.geom.surface(face.surface) else {
            err = Some(CentroidError::UnsupportedCurvedFace);
            return Ok(());
        };
        accumulate_planar_centroid(brep, face, &mut vol, &mut moment);
        Ok(())
    })
    .map_err(|_: AreaError| CentroidError::UnsupportedCurvedFace)?;
    if let Some(e) = err {
        return Err(e);
    }
    if vol.abs() <= tol.length * tol.length * tol.length {
        return Err(CentroidError::EmptySolid);
    }
    Ok(Point3::origin() + moment * (1.0 / vol))
}

/// Add a planar face's signed-tetrahedron volume and first moment to the running
/// totals. Outer and hole loops are both fanned; the B-rep stores each loop
/// already wound to its outward orientation (the same convention the volume
/// integral relies on), so no per-`sense` flip is applied — a hole loop's
/// reversed winding makes its tetra volumes negative and subtracts automatically.
fn accumulate_planar_centroid(brep: &Brep, face: &Face, vol: &mut f64, moment: &mut Vec3) {
    let mut loops = Vec::with_capacity(1 + face.inners.len());
    loops.push(face.outer);
    loops.extend(face.inners.iter().copied());
    for loop_id in loops {
        let verts = loop_points(brep, loop_id);
        if verts.len() < 3 {
            continue;
        }
        let o = Point3::origin();
        let p0 = verts[0];
        for i in 1..verts.len() - 1 {
            let a = verts[i];
            let b = verts[i + 1];
            // Signed volume of tetra (o, p0, a, b) = p0 · (a × b) / 6. The loop is
            // wound to the outward normal, so this signs correctly without a sense
            // flip; the barycentre is the tetra's vertex mean.
            let tv = (p0 - o).dot((a - o).cross(b - o)) / 6.0;
            let bary = o + ((p0 - o) + (a - o) + (b - o)) * 0.25;
            *vol += tv;
            *moment = *moment + (bary - o) * tv;
        }
    }
}

// ── formwork ─────────────────────────────────────────────────────────────────

/// The formwork (shuttering) area of a member, split by face orientation
/// (`DESIGN.md` §6-4, 公共建築数量積算基準).
///
/// Formwork is the area of concrete faces a form must be built against. The two
/// fields separate the two practical kinds because they are detailed and priced
/// differently:
///
/// * `side` — vertical faces (the member's sides), measured where the outward
///   normal is horizontal (within the angular tolerance of perpendicular to up).
/// * `bottom` — downward-facing faces (beam soffits, slab undersides), where the
///   outward normal points below horizontal.
///
/// Upward-facing (top) faces need no formwork and are excluded. Openings appear
/// as inner loops of the faces they pierce, so their area is already removed from
/// the net face area (the opening *reveals* — the form inside an opening — are a
/// separate, by-rule-often-omitted quantity and are not added here).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FormworkArea {
    /// Side (vertical-face) formwork area, in square metres.
    pub side: f64,
    /// Bottom (downward-face) formwork area, in square metres.
    pub bottom: f64,
}

/// Compute the [`FormworkArea`] of a B-rep (`DESIGN.md` §6-4).
///
/// Each face's outward normal is classified against world up (`+Z`):
///
/// * `|n·ẑ|` within `tol.angular` of zero ⇒ the face is vertical ⇒ `side`.
/// * `n·ẑ < −tol.angular` ⇒ the face faces down ⇒ `bottom`.
/// * `n·ẑ > +tol.angular` ⇒ the face faces up ⇒ no formwork.
///
/// A cylinder side face's outward normal is not constant, so it cannot be
/// classified by a single value: a vertical-axis column wall is all `side`, while
/// a horizontal-axis round beam's wall faces down on its lower half (a soffit ⇒
/// `bottom`) and up on its upper half (no formwork). Such a face is split by the
/// sign of the outward normal's `ẑ`-component, integrated over the patch (see
/// [`cylinder_formwork_split`]). The classification uses each face's *outward*
/// normal, so [`Sense`] is honoured.
///
/// # Errors
///
/// [`AreaError`] for the first face whose area / orientation cannot be computed.
pub fn formwork_area(brep: &Brep, tol: &Tol) -> Result<FormworkArea, AreaError> {
    let mut side = 0.0_f64;
    let mut bottom = 0.0_f64;
    for_each_face(brep, |face| {
        let area = face_area(brep, face)?;
        match face_up_component(brep, face) {
            Some(nz) => {
                if nz.abs() <= tol.angular {
                    side += area;
                } else if nz < 0.0 {
                    bottom += area;
                }
                // nz > 0: upward face, no formwork.
            }
            None => {
                // A cylinder side face. Its outward normal is not constant, so it
                // cannot be classified by a single `nz`: a vertical column wall is
                // all `side`, but a horizontal round beam's wall faces down on its
                // lower half (a soffit ⇒ `bottom`) and up on its upper half (no
                // formwork). Split the patch area by the sign of the outward
                // normal's `ẑ`-component, integrated over the face, rather than
                // silently dumping it all into `side`.
                let (s, b) = cylinder_formwork_split(brep, face, area, tol)
                    .ok_or(AreaError::UnsupportedCylinderFace)?;
                side += s;
                bottom += b;
            }
        }
        Ok(())
    })?;
    Ok(FormworkArea { side, bottom })
}

/// The `ẑ`-component of a face's outward normal, or `None` for a cylinder face
/// (whose normal is not constant — the caller handles it as a vertical wall when
/// the axis is vertical).
fn face_up_component(brep: &Brep, face: &Face) -> Option<f64> {
    match brep.geom.surface(face.surface)? {
        SurfaceGeom::Plane(plane) => Some(planar_outward_normal(plane, face.sense).z),
        SurfaceGeom::Cylinder(_) => None,
    }
}

/// Split a cylinder side face's area into `(side, bottom)` formwork by the sign of
/// its outward normal's `ẑ`-component, integrated over the patch.
///
/// The outward normal at angle `φ` about the cylinder's `(û, v̂)` basis is
/// `n̂(φ) = cosφ û + sinφ v̂` (always perpendicular to the axis), so its
/// `ẑ`-component is `nz(φ) = (û·ẑ) cosφ + (v̂·ẑ) sinφ`. The local axial height of
/// the patch at angle `φ` is `h(φ)`: a constant `L` for a straight patch, or the
/// cut-plane height `z₁(φ)` for an obliquely cut patch.
///
/// * **Vertical axis** — `û, v̂` are horizontal, so `nz ≡ 0`: the whole wall is
///   vertical and counts as `side` (a round column).
/// * **Non-vertical axis** — the down-facing part (`nz < 0`) of the wall is a
///   soffit and its area `r ∫_{nz<0} h(φ) dφ` is `bottom`; the up-facing part
///   needs no formwork; a tilted axis additionally contributes its horizontal
///   `(nz ≈ 0)` extremes, but those are a measure-zero seam and fold into the
///   dominant down/up split.
///
/// Returns `None` for a cylinder face whose rim pair is not one the closed form
/// recognises (mirroring the area/volume integrals), so the caller reports it
/// rather than silently misclassifying it.
fn cylinder_formwork_split(brep: &Brep, face: &Face, area: f64, tol: &Tol) -> Option<(f64, f64)> {
    let SurfaceGeom::Cylinder(cyl) = brep.geom.surface(face.surface)? else {
        return None;
    };
    let lp = brep.topo.loops.get(face.outer)?;
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
    let (u, v) = plane_basis(axis);
    // nz(φ) = a_n cosφ + b_n sinφ, the outward normal's ẑ-component.
    let a_n = u.z;
    let b_n = v.z;

    // Vertical axis: nz ≡ 0, the whole wall is vertical ⇒ all `side`. (A round
    // column.) `nz`'s amplitude is √(a_n²+b_n²) = |sin(angle between axis and ẑ)|,
    // so a near-zero amplitude within the angular tolerance means a vertical axis.
    let amp = (a_n * a_n + b_n * b_n).sqrt();
    if amp <= tol.angular {
        return Some((area, 0.0_f64));
    }

    // Height profile h(φ) of the patch at angle φ:
    //   straight (2 circular rims)     → constant L,
    //   oblique  (1 circle + 1 ellipse, possibly clipped into 2 arcs) → z₁(φ).
    enum Height {
        Constant(f64),
        /// z₁(φ) = k − r(p_u cosφ + p_v sinφ).
        Oblique {
            k: f64,
            p_u: f64,
            p_v: f64,
        },
    }
    let (height, arcs): (Height, Vec<(f64, f64)>) = match (circles.len(), ellipses.len()) {
        (2, 0) => {
            let (c0, c1) = (circles[0].2.center(), circles[1].2.center());
            let length = (c1 - c0).dot(axis_vec).abs();
            // A straight patch's two rims share the same angular span; use either.
            (Height::Constant(length), vec![(circles[0].0, circles[0].1)])
        }
        (1, 1) | (2, 1) => {
            let ell = ellipses[0].2;
            let cut_plane = Plane::new(ell.center(), ell.normal().as_vec()).ok()?;
            let n_p = cut_plane.normal().as_vec();
            let denom = n_p.dot(axis_vec);
            if denom.abs() <= f64::EPSILON {
                return None;
            }
            // All circle arcs here are arcs of the one bottom rim, so they share a
            // centre; z₁ is anchored there.
            let bottom_centre = circles[0].2.center();
            let k = n_p.dot(cut_plane.point() - bottom_centre) / denom;
            let p_u = n_p.dot(u) / denom;
            let p_v = n_p.dot(v) / denom;
            let spans = circles.iter().map(|&(p0, p1, _)| (p0, p1)).collect();
            (Height::Oblique { k, p_u, p_v }, spans)
        }
        _ => return None,
    };

    // ∫ r·h(φ) dφ over [a, b], the patch area between angles a and b.
    let height_integral = |a: f64, b: f64| -> f64 {
        match &height {
            Height::Constant(l) => r * l * (b - a),
            Height::Oblique { k, p_u, p_v } => {
                r * (k * (b - a) - r * (p_u * (b.sin() - a.sin()) + p_v * (a.cos() - b.cos())))
            }
        }
    };

    // nz(φ) = amp·cos(φ − ψ); it is < 0 exactly on (ψ + π/2, ψ + 3π/2). Integrate
    // the patch area over each arc clipped to that down-facing band to get the
    // soffit (`bottom`); the rest of the integrated area is up-facing (no
    // formwork). The seam where nz = 0 has measure zero, so it needs no `side`
    // bucket for a non-vertical axis.
    let psi = b_n.atan2(a_n);
    let down_lo = psi + std::f64::consts::FRAC_PI_2;
    let down_hi = psi + 3.0_f64 * std::f64::consts::FRAC_PI_2;

    let mut bottom = 0.0_f64;
    let mut swept = 0.0_f64;
    for &(phi0, phi1) in &arcs {
        bottom += integrate_over_down_band(phi0, phi1, down_lo, down_hi, &height_integral);
        swept += height_integral(phi0, phi1);
    }

    // For a non-vertical wall every patch point faces strictly up or down (the
    // horizontal seam is a measure-zero curve), so there is **no** `side`
    // contribution: the down-facing part is the soffit (`bottom`) and the
    // up-facing part needs no formwork. The integrated `swept` area equals the
    // exact `face_area` only up to closed-form rounding, so rescale `bottom` by
    // the exact `area / swept` to anchor it to the true patch area.
    if swept > 0.0_f64 {
        bottom *= area / swept;
    }
    Some((0.0_f64, bottom))
}

/// Integrate `f` (a function returning `∫ r·h dφ` over an interval) over the part
/// of `[phi0, phi1]` that lies in the down-facing angular band `[down_lo, down_hi]`
/// (an interval of length `π`), accounting for the `2π` periodicity of the band.
fn integrate_over_down_band(
    phi0: f64,
    phi1: f64,
    down_lo: f64,
    down_hi: f64,
    f: &impl Fn(f64, f64) -> f64,
) -> f64 {
    use std::f64::consts::PI;
    let two_pi = 2.0_f64 * PI;
    // Slide the arc start into [down_lo, down_lo + 2π) and the band to the same
    // origin, then clip [a, b] against the two band copies [down_lo, down_hi] and
    // [down_lo + 2π, down_hi + 2π] that can overlap a ≤ 2π-long arc.
    let shift = ((phi0 - down_lo) / two_pi).floor() * two_pi;
    let a = phi0 - shift;
    let b = phi1 - shift;
    let lo = down_lo - shift; // == down_lo's representative ≤ a
    let hi = down_hi - shift;
    let mut total = 0.0_f64;
    for band in [(lo, hi), (lo + two_pi, hi + two_pi)] {
        let s = a.max(band.0);
        let e = b.min(band.1);
        if e > s {
            // Map the clipped sub-interval back to the original angle frame for
            // the height integral (which depends on absolute φ).
            total += f(s + shift, e + shift);
        }
    }
    total
}

// ── shared face walk ─────────────────────────────────────────────────────────

/// Drive `f` over every face of every solid in the B-rep, short-circuiting on the
/// first error.
fn for_each_face<E>(brep: &Brep, mut f: impl FnMut(&Face) -> Result<(), E>) -> Result<(), E> {
    for &solid_id in &brep.solids {
        let Some(solid) = brep.topo.solids.get(solid_id) else {
            continue;
        };
        for &shell_id in &solid.shells {
            let Some(shell) = brep.topo.shells.get(shell_id) else {
                continue;
            };
            for &face_id in &shell.faces {
                if let Some(face) = brep.topo.faces.get(face_id) {
                    f(face)?;
                }
            }
        }
    }
    Ok(())
}

/// Collect a loop's start-vertex coordinates.
fn loop_points(brep: &Brep, loop_id: Id<Loop>) -> Vec<Point3> {
    let Some(lp) = brep.topo.loops.get(loop_id) else {
        return Vec::new();
    };
    lp.half_edges
        .iter()
        .filter_map(|&he_id| {
            let he = brep.topo.half_edges.get(he_id)?;
            let v = brep.topo.vertices.get(he.start)?;
            brep.geom.point(v.point).and_then(|g| g.as_point())
        })
        .collect()
}
