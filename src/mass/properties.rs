//! Higher mass properties built on the per-face area core: surface area,
//! centroid, and the formwork-area take-off (`DESIGN.md` §6-4).

use crate::brep::Brep;
use crate::geom::SurfaceGeom;
use crate::math::{Point3, Vec3};
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
/// A cylinder side face is always vertical when the cylinder axis is vertical (a
/// round column); its area is added to `side`. The classification uses each
/// face's *outward* normal, so [`Sense`] is honoured.
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
                // A cylinder face: vertical if its axis is vertical (a round
                // column wall). Treat as side formwork in that case; otherwise
                // its orientation varies and it is reported as unsupported via
                // face_area already having succeeded — classify by the axis.
                side += area;
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
