//! Quantity take-off: concrete volume and formwork area from an evaluated
//! member (`DESIGN.md` §6-2, §6-4, §7).
//!
//! The take-off runs on the member's *clipped* B-rep — the geometry after the
//! priority deductions have been applied (a girder already trimmed to its column
//! inner-face length). Concrete volume is then the plain signed volume.
//!
//! Formwork is the side / bottom face-area split, but with one domain rule: the
//! girder's end faces where a column trimmed it are *contact* faces (the girder
//! is embedded in the column there) and carry no form (公共建築数量積算基準).
//! They are excluded by the geometric test "stepping a hair outward from the
//! face lands inside a clipper prism", so what remains is exactly the exposed
//! web sides and soffit.

use crate::boolean::prismatic::ExtrudeLeaf;
use crate::csg::{CsgNode, EvalError, Profile2d, StableId};
use crate::geom::SurfaceGeom;
use crate::mass::signed_volume_checked;
use crate::math::{Point3, Vec3};
use crate::primitives::plane_basis;
use crate::tolerance::Tol;
use crate::topo::{Face, Sense};

use super::{occupancy_leaf, Model};

/// The formwork area of a take-off, split by face orientation.
///
/// Mirrors [`crate::mass::FormworkArea`] but is re-exposed through the take-off
/// so callers import one quantity-take-off vocabulary.
pub use crate::mass::FormworkArea;

/// The concrete quantity take-off of one member (`DESIGN.md` §7).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct QuantityTakeoff {
    /// Concrete volume after priority deductions, in cubic metres.
    pub concrete_volume: f64,
    /// Side (vertical-face) formwork area, in square metres.
    pub formwork_side: f64,
    /// Bottom (downward-face) formwork area, in square metres.
    pub formwork_bottom: f64,
}

impl QuantityTakeoff {
    /// Total formwork area (side + bottom), in square metres.
    pub fn formwork_total(&self) -> f64 {
        self.formwork_side + self.formwork_bottom
    }
}

/// Compute the [`QuantityTakeoff`] of the member `id` in `model`.
///
/// The member is evaluated to its clipped B-rep (resolving its priority
/// deductions against the other members), then the concrete volume and formwork
/// split are read off it. The B-rep is built with the checked volume integral,
/// so an unsupported face surfaces as an error rather than a wrong quantity.
///
/// # Formwork and column contact
///
/// Formwork is only counted on *exposed* concrete faces. A clipped girder's end
/// faces — created where a column trimmed the girder — abut the column and carry
/// no form; the take-off excludes them by the geometric test "stepping a hair
/// outward from the face lands inside a clipper" (公共建築数量積算基準: 柱との
/// 接触面は型枠から控除). Side formwork is then the two vertical web faces over
/// the inner-clear length, and bottom formwork the soffit, exactly the hand
/// calculation.
///
/// The signature takes `&mut Model` to match the lazy-evaluation contract (a
/// take-off may trigger a re-evaluation), though the current implementation does
/// not mutate the model; the `&mut` keeps the API stable if per-member caching
/// is later threaded through the model.
///
/// # Errors
///
/// * Any [`EvalError`] the member's evaluation raises (unknown clipper, cyclic
///   dependency, unsupported boolean, …).
/// * A volume/area failure is carried as an [`EvalError::Construction`] string
///   when an evaluated face cannot be integrated.
pub fn takeoff(model: &mut Model, id: StableId, tol: &Tol) -> Result<QuantityTakeoff, EvalError> {
    let brep = model.evaluate(id, tol)?;
    let concrete_volume = signed_volume_checked(&brep)
        .map_err(|e| EvalError::Construction(format!("volume integral failed: {e:?}")))?;

    // The clippers this member deducts against; their prisms identify contact
    // faces to exclude from formwork.
    let clippers = model.clipper_leaves_of(id);

    let mut side = 0.0_f64;
    let mut bottom = 0.0_f64;
    let mut err: Option<String> = None;
    for_each_face(&brep, |face| {
        // A contact face abuts a clipper: stepping outward lands inside one.
        if let Some((nz, centroid)) = face_outward(&brep, face) {
            if is_contact_face(centroid, nz_dir(&brep, face), &clippers, tol) {
                return;
            }
            match crate::mass::face_area_of(&brep, face) {
                Ok(area) => {
                    if nz.abs() <= tol.angular {
                        side += area;
                    } else if nz < 0.0 {
                        bottom += area;
                    }
                }
                Err(e) => {
                    if err.is_none() {
                        err = Some(format!("formwork area failed: {e:?}"));
                    }
                }
            }
        } else {
            // Cylinder (round column) wall: vertical, all side formwork, unless
            // its whole patch abuts a clipper (rare for a round column).
            match crate::mass::face_area_of(&brep, face) {
                Ok(area) => side += area,
                Err(e) => {
                    if err.is_none() {
                        err = Some(format!("formwork area failed: {e:?}"));
                    }
                }
            }
        }
    });
    if let Some(e) = err {
        return Err(EvalError::Construction(e));
    }

    Ok(QuantityTakeoff {
        concrete_volume,
        formwork_side: side,
        formwork_bottom: bottom,
    })
}

/// Walk every face of a B-rep.
fn for_each_face(brep: &crate::brep::Brep, mut f: impl FnMut(&Face)) {
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
                    f(face);
                }
            }
        }
    }
}

/// The `(ẑ-component of outward normal, face centroid)` for a planar face, or
/// `None` for a cylinder face.
fn face_outward(brep: &crate::brep::Brep, face: &Face) -> Option<(f64, Point3)> {
    let SurfaceGeom::Plane(plane) = brep.geom.surface(face.surface)? else {
        return None;
    };
    let n = match face.sense {
        Sense::Same => plane.normal().as_vec(),
        Sense::Reversed => -plane.normal().as_vec(),
    };
    Some((n.z, face_centroid(brep, face)))
}

/// The outward normal direction of a planar face (unit), or zero for cylinders.
fn nz_dir(brep: &crate::brep::Brep, face: &Face) -> Vec3 {
    match brep.geom.surface(face.surface) {
        Some(SurfaceGeom::Plane(plane)) => match face.sense {
            Sense::Same => plane.normal().as_vec(),
            Sense::Reversed => -plane.normal().as_vec(),
        },
        _ => Vec3::ZERO,
    }
}

/// The centroid of a face's outer loop (vertex average — adequate for the
/// step-outward contact probe).
fn face_centroid(brep: &crate::brep::Brep, face: &Face) -> Point3 {
    let Some(lp) = brep.topo.loops.get(face.outer) else {
        return Point3::origin();
    };
    let mut acc = Vec3::ZERO;
    let mut n = 0.0_f64;
    for &he_id in &lp.half_edges {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        let Some(v) = brep.topo.vertices.get(he.start) else {
            continue;
        };
        if let Some(p) = brep.geom.point(v.point).and_then(|g| g.as_point()) {
            acc = acc + (p - Point3::origin());
            n += 1.0;
        }
    }
    if n == 0.0 {
        Point3::origin()
    } else {
        Point3::origin() + acc * (1.0 / n)
    }
}

/// `true` if stepping a hair (`100·tol.length`) outward from `centroid` along
/// `outward` lands inside any clipper prism — the signature of a contact face.
fn is_contact_face(centroid: Point3, outward: Vec3, clippers: &[ExtrudeLeaf], tol: &Tol) -> bool {
    let step = 100.0 * tol.length;
    let probe = centroid + outward * step;
    clippers.iter().any(|c| point_in_leaf(probe, c, tol))
}

/// `true` if `p` lies inside the (rectangular) extrusion `leaf`, within `tol`.
///
/// Round (circular) clippers use the radial test. The point is projected onto
/// the extrusion axis (must be within `[−tol, length+tol]`) and onto the profile
/// frame (must be within the half-extents).
fn point_in_leaf(p: Point3, leaf: &ExtrudeLeaf, tol: &Tol) -> bool {
    let Some(axis) = leaf.axis.try_unit() else {
        return false;
    };
    let av = axis.as_vec();
    let rel = p - leaf.origin;
    let t = rel.dot(av);
    if t < -tol.length || t > leaf.length + tol.length {
        return false;
    }
    let (u, v) = plane_basis(axis);
    let pu = rel.dot(u);
    let pv = rel.dot(v);
    match leaf.profile {
        Profile2d::Rect { half_w, half_h } => {
            // Profile frame: u ↔ half_w, v ↔ half_h (matches the extruder).
            pu.abs() <= half_w + tol.length && pv.abs() <= half_h + tol.length
        }
        Profile2d::Circle { radius } => (pu * pu + pv * pv).sqrt() <= radius + tol.length,
        // An H-section clipper for contact exclusion is out of scope here; treat
        // its bounding box conservatively (no false contact beyond the box).
        Profile2d::HSection { half_w, half_h, .. } => {
            pu.abs() <= half_w + tol.length && pv.abs() <= half_h + tol.length
        }
    }
}

impl Model {
    /// The gross extrusion leaves of the clippers that member `id` deducts.
    fn clipper_leaves_of(&self, id: StableId) -> Vec<ExtrudeLeaf> {
        let mut out = Vec::new();
        if let Some(m) = self.get(id) {
            collect_clipper_leaves(self, m.csg(), &mut out);
        }
        out
    }
}

/// Collect, from a member's CSG tree, the occupancy leaves of every clipper it
/// references.
fn collect_clipper_leaves(model: &Model, node: &CsgNode, out: &mut Vec<ExtrudeLeaf>) {
    if let CsgNode::Clip { base, clippers, .. } = node {
        for &cid in clippers {
            if let Some(cm) = model.get(cid) {
                if let Some(leaf) = occupancy_leaf(cm.csg()) {
                    out.push(leaf);
                }
            }
        }
        collect_clipper_leaves(model, base, out);
    }
}
