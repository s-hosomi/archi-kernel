//! Axis-aligned bounding boxes for the coarse clash phase.
//!
//! The coarse phase rejects member pairs whose bounding boxes do not overlap
//! before any (expensive) exact intersection is attempted (`DESIGN.md` §6-6,
//! `docs/research/06-domain.md` §6.1: hard-clash detection is an AABB pre-filter
//! followed by a volume test). The box is built from a member's *evaluated*
//! B-rep so it reflects the real placed geometry, openings and clips included.
//!
//! Two subtleties make a naive "min/max of the vertices" wrong for this kernel:
//!
//! * **Cylinder walls carry their extent in the surface, not the vertices.** A
//!   round column's B-rep has explicit vertices only at the cap seams; the radial
//!   bulge of the wall is implicit in the [`Cylinder`](crate::primitives::Cylinder)
//!   surface. Taking only the vertices would under-bound the column by almost its
//!   whole radius. Each cylinder face therefore contributes the exact AABB of its
//!   finite axial segment expanded by the radius (`DESIGN.md` §6-6: "円筒は
//!   中心±半径の拡張").
//! * **An empty B-rep has no box.** A member that evaluated to nothing yields
//!   `None` rather than a degenerate point box, so the caller skips it instead of
//!   reporting a spurious clash at the origin.

use crate::brep::Brep;
use crate::geom::{CurveGeom, SurfaceGeom};
use crate::math::Point3;
use crate::primitives::Cylinder;
use crate::tolerance::Tol;

/// An axis-aligned bounding box in world coordinates (metres).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Aabb {
    /// The minimum corner `(x, y, z)`.
    pub min: [f64; 3],
    /// The maximum corner `(x, y, z)`.
    pub max: [f64; 3],
}

impl Aabb {
    /// Start an empty (inverted) box that any point expands.
    fn empty() -> Self {
        Self {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }

    /// `true` if no point has been added (the box is still inverted).
    fn is_unset(&self) -> bool {
        self.min[0] > self.max[0]
    }

    /// Grow the box to include the world point `p`.
    fn add_point(&mut self, p: Point3) {
        self.extend([p.x, p.y, p.z], [p.x, p.y, p.z]);
    }

    /// Grow the box to include the axis-aligned interval `[lo, hi]` per axis.
    fn extend(&mut self, lo: [f64; 3], hi: [f64; 3]) {
        for ((mn, mx), (l, h)) in self
            .min
            .iter_mut()
            .zip(self.max.iter_mut())
            .zip(lo.into_iter().zip(hi))
        {
            if l < *mn {
                *mn = l;
            }
            if h > *mx {
                *mx = h;
            }
        }
    }

    /// `true` if this box and `other` overlap, treating a separation of up to
    /// `tol.length` as still overlapping.
    ///
    /// The tolerant comparison means a face-to-face touch (zero true overlap)
    /// counts as an overlap in the coarse phase, so the fine phase still gets the
    /// chance to classify it as [`Touching`](super::ClashKind::Touching) rather
    /// than the coarse phase discarding it (`DESIGN.md` §6-6: no silent "no
    /// clash").
    pub(crate) fn overlaps(&self, other: &Aabb, tol: &Tol) -> bool {
        let sep = |a_max: f64, b_min: f64| a_max + tol.length < b_min;
        !self
            .max
            .iter()
            .zip(other.min.iter())
            .zip(other.max.iter().zip(self.min.iter()))
            .any(|((&smax, &omin), (&omax, &smin))| sep(smax, omin) || sep(omax, smin))
    }
}

/// The world-space AABB of an evaluated B-rep, or `None` if it is empty.
///
/// Every explicit vertex is included, and every cylinder face additionally
/// contributes the box of its finite axial segment expanded by the cylinder
/// radius (see the module docs). Returns `None` when the B-rep contributes no
/// geometry at all, so the caller skips empty members.
pub(crate) fn aabb_of(brep: &Brep) -> Option<Aabb> {
    let mut bb = Aabb::empty();

    // Every explicit vertex.
    for (_, vert) in brep.topo.vertices.iter() {
        if let Some(p) = brep.geom.point(vert.point).and_then(|g| g.as_point()) {
            bb.add_point(p);
        }
    }

    // Every cylinder face: expand by the radius over its axial segment.
    for (_, face) in brep.topo.faces.iter() {
        if let Some(SurfaceGeom::Cylinder(cyl)) = brep.geom.surface(face.surface) {
            if let Some((c0, c1)) = cylinder_face_segment(brep, face) {
                add_cylinder_segment(&mut bb, cyl, c0, c1);
            }
        }
    }

    if bb.is_unset() {
        None
    } else {
        Some(bb)
    }
}

/// The two cap-circle centres bounding a cylinder face's axial segment, taken
/// from the circular rim arcs of its outer loop. `None` if the face has no
/// recognisable circular rims.
fn cylinder_face_segment(brep: &Brep, face: &crate::topo::Face) -> Option<(Point3, Point3)> {
    let lp = brep.topo.loops.get(face.outer)?;
    let mut centres: Vec<Point3> = Vec::new();
    for &he_id in &lp.half_edges {
        let Some(he) = brep.topo.half_edges.get(he_id) else {
            continue;
        };
        if let Some(CurveGeom::Circle(circle)) = brep.geom.curve(he.curve) {
            let c = circle.center();
            if !centres.iter().any(|p| (*p - c).norm() <= f64::EPSILON) {
                centres.push(c);
            }
        }
    }
    match centres.as_slice() {
        [a, b] => Some((*a, *b)),
        [a] => Some((*a, *a)),
        _ => None,
    }
}

/// Expand `bb` to bound a finite cylinder segment from `c0` to `c1` (cap
/// centres) of the given radius.
///
/// For a cylinder of unit axis `a` and radius `r`, the half-extent of an
/// infinite cylinder along world axis `eᵢ` is `r·√(1 − (a·eᵢ)²)` — the radius of
/// the silhouette circle projected onto that axis. Adding this around both cap
/// centres bounds the finite segment exactly (caps included).
fn add_cylinder_segment(bb: &mut Aabb, cyl: &Cylinder, c0: Point3, c1: Point3) {
    let r = cyl.radius();
    let a = cyl.axis().dir().as_vec();
    // Per-axis silhouette half-extent: r·√(1 − aᵢ²); clamp tiny negative rounding.
    let ext = [
        r * (1.0 - a.x * a.x).max(0.0).sqrt(),
        r * (1.0 - a.y * a.y).max(0.0).sqrt(),
        r * (1.0 - a.z * a.z).max(0.0).sqrt(),
    ];
    for centre in [c0, c1] {
        let cc = [centre.x, centre.y, centre.z];
        let lo = [cc[0] - ext[0], cc[1] - ext[1], cc[2] - ext[2]];
        let hi = [cc[0] + ext[0], cc[1] + ext[1], cc[2] + ext[2]];
        bb.extend(lo, hi);
    }
}
