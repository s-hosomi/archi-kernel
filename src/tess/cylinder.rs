//! Cylinder-face tessellation by stitching its two rim polylines.
//!
//! A cylinder side face is a ruled patch bounded by two **rim** arcs (the bottom
//! and top edges) and two **seam** lines (the vertical rulings). Watertightness
//! requires the rim vertices to coincide with those of the adjacent cap faces,
//! so each rim is sampled **from its own boundary half-edge's shared curve** at
//! the same per-curve density the cap uses (`DESIGN.md` §6-5): the side patch's
//! bottom-rim half-edge is the sibling of the bottom cap's arc (same curve,
//! reversed boundary ⇒ identical sample points), and likewise the top rim with
//! the top cap — circular (straight patch) or elliptical (oblique cut).
//!
//! The two rim polylines run between the same two seam endpoints. They may have
//! different vertex counts (a circle rim and an ellipse rim discretise to
//! different numbers of segments), so they are stitched by a parameter-merge
//! triangle strip: advancing whichever rim has covered less of its arc keeps the
//! triangles well-shaped and uses every rim vertex, guaranteeing the rim
//! vertices match the caps exactly. The interior stitch diagonals are internal
//! to the face and do not affect cross-face watertightness.
//!
//! The triangle winding is chosen so the emitted normal points radially outward
//! (the cylinder's outward normal, folded through the face
//! [`Sense`](crate::topo::Sense)). The strip is built once with a single,
//! topologically consistent base winding; the patch is then oriented as a whole
//! by an **area-weighted** outward test (`Σ faceₙ · n_out`), so a sliver
//! triangle near the seam — whose own face normal is dominated by round-off —
//! cannot flip independently of its neighbours and break watertightness.

use crate::brep::Brep;
use crate::geom::CurveGeom;
use crate::math::{Point3, Vec3};
use crate::primitives::Cylinder;
use crate::tolerance::Tol;
use crate::topo::{Face, Sense};

use super::intern::MeshBuilder;
use super::{arc_segment_count, TessError, TessOptions};

/// A rim of the patch: a polyline of 3-D points plus the parameter fraction
/// `0..=1` at each point (used to merge two rims of unequal length).
struct Rim {
    pts: Vec<Point3>,
    frac: Vec<f64>,
}

/// Tessellate one cylinder side face into the builder.
pub(crate) fn tessellate_cylinder_face(
    brep: &Brep,
    face: &Face,
    cyl: &Cylinder,
    builder: &mut MeshBuilder,
    face_tag: u32,
    opts: &TessOptions,
    _tol: &Tol,
) -> Result<(), TessError> {
    let axis = cyl.axis().dir().as_vec();
    let lp = brep
        .topo
        .loops
        .get(face.outer)
        .ok_or(TessError::DanglingReference)?;

    // Collect the rim arcs (circle / ellipse half-edges). A straight patch has
    // two circular rims; an oblique patch one circle + one ellipse. The seam
    // lines (straight half-edges) are implied by the rim endpoints and need no
    // explicit handling — their endpoints are shared rim/cap vertices.
    let mut rims: Vec<(f64, Rim)> = Vec::new(); // (axial projection of a rim point, rim)
    for &he_id in &lp.half_edges {
        let he = brep
            .topo
            .half_edges
            .get(he_id)
            .ok_or(TessError::DanglingReference)?;
        let curve = brep
            .geom
            .curve(he.curve)
            .ok_or(TessError::DanglingReference)?;
        let [a, b] = he.boundary;
        match curve {
            CurveGeom::Circle(_) | CurveGeom::Ellipse(_) => {
                let radius = match curve {
                    CurveGeom::Circle(c) => c.radius(),
                    CurveGeom::Ellipse(e) => e.semi_major(),
                    CurveGeom::Line(_) => unreachable!(),
                };
                let segs = arc_segment_count(radius, b - a, opts.chord_tolerance);
                // Sample the rim curve inclusive of both endpoints (segs + 1
                // points), so the seam endpoints are present and shared.
                let mut pts = Vec::with_capacity(segs + 1);
                let mut frac = Vec::with_capacity(segs + 1);
                for s in 0..=segs {
                    let f = (s as f64) / (segs as f64);
                    let t = a + (b - a) * f;
                    pts.push(curve.point_at(t));
                    frac.push(f);
                }
                let axial = (pts[0] - Point3::origin()).dot(axis);
                rims.push((axial, Rim { pts, frac }));
            }
            CurveGeom::Line(_) => { /* seam ruling */ }
        }
    }

    if rims.len() != 2 {
        return Err(TessError::UnsupportedSurface);
    }

    // The lower rim (smaller axial projection of its start) is the "bottom".
    // Orient both rims to run between the same two seam endpoints: the top rim
    // is traversed opposite the bottom (the loop walks bottom one way, top the
    // other), so reverse the top so both start at the same seam point.
    let (mut a, mut b) = (rims.pop().unwrap(), rims.pop().unwrap());
    let (bottom, top) = if a.0 <= b.0 {
        (&mut a.1, &mut b.1)
    } else {
        (&mut b.1, &mut a.1)
    };
    // Reverse the top rim so both polylines start at the same seam point.
    top.pts.reverse();
    for f in &mut top.frac {
        *f = 1.0 - *f;
    }
    top.frac.reverse();

    // Stitch the two polylines into a triangle strip. Both run seamA → seamB and
    // share endpoints; advance whichever has covered less of its arc. Every
    // triangle is emitted with the strip's single consistent base winding, then
    // the whole patch is oriented once to the cylinder's outward normal (the
    // radial direction folded through the face sense), so winding is correct
    // regardless of how the rim half-edges were stored — and no individual
    // triangle can flip against its neighbours.
    let reversed = matches!(face.sense, Sense::Reversed);
    let axis_origin = cyl.axis().origin();
    let orient = Orienter {
        axis_origin,
        axis_dir: axis,
        reversed,
    };
    stitch(bottom, top, builder, face_tag, &orient);

    Ok(())
}

/// Decides triangle winding by the cylinder's outward normal.
struct Orienter {
    axis_origin: Point3,
    axis_dir: Vec3,
    reversed: bool,
}

impl Orienter {
    /// The outward normal at point `p`: the radial direction from the axis to
    /// `p`, negated when the face sense reverses the surface normal.
    fn outward_at(&self, p: Point3) -> Vec3 {
        let rel = p - self.axis_origin;
        let along = self.axis_dir * rel.dot(self.axis_dir);
        let radial = rel - along;
        let n = radial
            .try_unit()
            .map(|u| u.as_vec())
            .unwrap_or(self.axis_dir);
        if self.reversed {
            -n
        } else {
            n
        }
    }
}

/// Stitch two rim polylines (same endpoints, possibly unequal lengths) into a
/// triangle strip.
///
/// `bottom.pts[0] == top.pts[0]` and `bottom.pts[last] == top.pts[last]` (the
/// two seam points). We walk both from index 0, each step forming one triangle
/// with the current bottom point `bi`, current top point `ti`, and the next
/// point on whichever rim has the smaller parameter fraction — a parameter
/// merge that consumes every vertex of both rims.
///
/// All triangles are collected with one consistent base winding (`(b0, b1, t)`
/// for a bottom advance, `(b, t1, t0)` for a top advance — both wind the same
/// way around the strip). The whole patch is then oriented to the outward
/// normal by a single **area-weighted** test: `Σ faceₙ · n_out` summed over all
/// triangles, where `faceₙ = (p1−p0)×(p2−p0)` has magnitude proportional to
/// twice the triangle area. A sliver near the seam contributes a vanishing
/// term, so the sound, area-bearing triangles set the sign; if it is negative
/// every triangle is reversed together. This keeps the strip internally
/// consistent (no neighbour can disagree) while still folding through the face
/// sense correctly.
fn stitch(bottom: &Rim, top: &Rim, builder: &mut MeshBuilder, face_tag: u32, orient: &Orienter) {
    let mut i = 0usize; // bottom index
    let mut j = 0usize; // top index
    let nb = bottom.pts.len();
    let nt = top.pts.len();

    // Collect the strip in its base winding, accumulating the area-weighted
    // agreement with the outward normal as we go.
    let mut tris: Vec<(Point3, Point3, Point3)> = Vec::with_capacity(nb + nt);
    let mut signed = 0.0_f64;

    while i + 1 < nb || j + 1 < nt {
        // Decide which rim to advance: the one whose next point has the smaller
        // parameter fraction (so the strip stays monotone). If one rim is
        // exhausted, advance the other.
        let advance_bottom = if i + 1 >= nb {
            false
        } else if j + 1 >= nt {
            true
        } else {
            bottom.frac[i + 1] <= top.frac[j + 1]
        };

        let (p0, p1, p2) = if advance_bottom {
            let tri = (bottom.pts[i], bottom.pts[i + 1], top.pts[j]);
            i += 1;
            tri
        } else {
            let tri = (bottom.pts[i], top.pts[j + 1], top.pts[j]);
            j += 1;
            tri
        };

        let centroid = Point3::new(
            (p0.x + p1.x + p2.x) / 3.0,
            (p0.y + p1.y + p2.y) / 3.0,
            (p0.z + p1.z + p2.z) / 3.0,
        );
        let n_out = orient.outward_at(centroid);
        let face_n = (p1 - p0).cross(p2 - p0);
        // Area-weighted because `face_n` magnitude is 2·area: a sliver's tiny,
        // round-off-dominated contribution cannot swing the patch's sign.
        signed += face_n.dot(n_out);
        tris.push((p0, p1, p2));
    }

    // Orient the whole patch once. `signed == 0` only for a degenerate (zero
    // area) patch, where either winding is equivalent — keep the base winding.
    let flip = signed < 0.0;
    for (p0, p1, p2) in tris {
        let a = builder.vertex(p0);
        let b = builder.vertex(p1);
        let c = builder.vertex(p2);
        if flip {
            builder.triangle(a, c, b, face_tag);
        } else {
            builder.triangle(a, b, c, face_tag);
        }
    }
}
