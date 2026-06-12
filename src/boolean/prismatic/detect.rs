//! Common prismatic direction detection and reprofiling (`DESIGN.md` §4.2).
//!
//! The 2.5-D reduction does **not** ask "are the two extrusion axes parallel?".
//! It asks "does a *common prismatic direction* exist?" A solid is prismatic
//! along a direction `d` when its cross-section perpendicular to `d` is constant
//! along `d`. For an extruded leaf this gives the candidate set:
//!
//! * **any profile**: the extrusion axis (the section is constant along it);
//! * **a rectangular profile (a box)**: *all three* of its axes — the two
//!   in-plane profile axes are prismatic too, because a box is a prism along
//!   each of its three edge directions. This is exactly what lets a rectangular
//!   column crossing a rectangular beam at a right angle fall to 2.5-D: take the
//!   beam axis as `d` and view the column as a prism along the beam axis
//!   (`DESIGN.md` §4.2, `synthesis.md` §2-6).
//!
//! A circular profile is prismatic only along its own axis; viewed along any
//! other direction its section is not constant. When a common direction is found
//! and a circular section is involved, the 2-D side would carry an arc, which is
//! Phase 3c — we report it as such so the caller maps it to
//! [`EvalError::Unsupported3dBoolean`](crate::csg::EvalError).
//!
//! Direction agreement is in **length units, not angle** (`DESIGN.md` §2.2): two
//! unit directions `a`, `b` agree when `sinθ · L_max ≤ Tol::length`, i.e.
//! `|a × b| · L_max ≤ tol.length`, where `L_max` is the largest dimension in
//! play. The cosine / `Tol::angular` shortcut is avoided because it is quadratic
//! near ±1 and lets millimetre projection errors through at building spans.

use crate::boolean::poly2d::{Contour, Region};
use crate::csg::Profile2d;
use crate::math::{Point3, Vec3};
use crate::primitives::plane_basis;
use crate::profile::ProfileGeom;
use crate::tolerance::Tol;

use super::error::{Operand, PrismError};

/// A deterministic right-handed 2-D frame for a prismatic direction `d`.
///
/// `d` is the prism direction (unit); `e1`, `e2` span the perpendicular plane
/// with `e1 × e2 = d`, taken from [`plane_basis`] so the frame is the same one
/// the extruder and cut use. A 3-D point `p` reprofiles to
/// `(p·e1, p·e2)` in 2-D with axial coordinate `p·d`; the inverse lift is
/// `origin + e1·x + e2·y + d·t` with `origin` the world origin, so **both
/// operands share one frame** (`DESIGN.md` §4.2 "同一の 2D フレーム").
#[derive(Debug, Clone, Copy)]
pub(crate) struct Frame {
    /// The prism direction (unit).
    pub d: Vec3,
    /// First in-plane basis vector (unit).
    pub e1: Vec3,
    /// Second in-plane basis vector (unit), with `e1 × e2 = d`.
    pub e2: Vec3,
}

impl Frame {
    /// Build the frame for a unit prism direction.
    fn new(d: Vec3) -> Self {
        let unit = d.try_unit().expect("prism direction is non-zero");
        let (e1, e2) = plane_basis(unit);
        Self {
            d: unit.as_vec(),
            e1,
            e2,
        }
    }

    /// Reprofile a 3-D point: `(x, y) = (p·e1, p·e2)`.
    #[inline]
    pub(crate) fn project(&self, p: Point3) -> [f64; 2] {
        let v = p - Point3::origin();
        [v.dot(self.e1), v.dot(self.e2)]
    }

    /// Axial coordinate of a 3-D point along `d`.
    #[inline]
    pub(crate) fn axial(&self, p: Point3) -> f64 {
        (p - Point3::origin()).dot(self.d)
    }

    /// Lift a 2-D frame point at axial height `t` back to 3-D.
    #[inline]
    pub(crate) fn lift(&self, xy: [f64; 2], t: f64) -> Point3 {
        Point3::origin() + self.e1 * xy[0] + self.e2 * xy[1] + self.d * t
    }
}

/// One operand reduced to a 2-D region plus an axial interval along the shared
/// prism direction.
#[derive(Debug, Clone)]
pub(crate) struct PrismOperand {
    /// The cross-section perpendicular to `d`, in the shared frame's 2-D coords.
    pub region: Region,
    /// Lower axial bound along `d`.
    pub t0: f64,
    /// Upper axial bound along `d`.
    pub t1: f64,
}

/// An extruded leaf described as an oriented box or a general prism, plus the
/// span dimensions used for length-based direction agreement.
#[derive(Debug, Clone)]
pub(crate) struct Leaf {
    /// World origin of the extrusion (bottom-cap centre).
    origin: Point3,
    /// Extrusion direction (unit).
    axis: Vec3,
    /// Extrusion length.
    length: f64,
    /// The profile.
    profile: Profile2d,
    /// In-plane profile basis `(u, v)` with `u × v = axis`.
    u: Vec3,
    v: Vec3,
}

impl Leaf {
    /// Describe an extruded leaf from its CSG fields.
    pub(crate) fn new(profile: Profile2d, origin: Point3, axis: Vec3, length: f64) -> Self {
        let unit = axis.try_unit().map(|u| u.as_vec()).unwrap_or(Vec3::Z);
        let (u, v) = plane_basis(unit.try_unit().unwrap_or(crate::math::Unit3::Z));
        Self {
            origin,
            axis: unit,
            length,
            profile,
            u,
            v,
        }
    }

    /// `true` if the profile is circular (round column / void).
    fn is_circular(&self) -> bool {
        matches!(self.profile, Profile2d::Circle { .. })
    }

    /// The candidate prismatic directions of this leaf (unit vectors).
    ///
    /// Always the extrusion axis; additionally the two in-plane profile axes for
    /// a rectangular profile (a box is a prism along all three of its axes).
    fn candidate_directions(&self) -> Vec<Vec3> {
        let mut dirs = vec![self.axis];
        if matches!(self.profile, Profile2d::Rect { .. }) {
            dirs.push(self.u);
            dirs.push(self.v);
        }
        dirs
    }

    /// The largest dimension of this leaf, for length-based direction agreement.
    fn max_dimension(&self) -> f64 {
        let prof = match self.profile {
            Profile2d::Rect { half_w, half_h } => (2.0 * half_w).max(2.0 * half_h),
            Profile2d::HSection { half_w, half_h, .. } => (2.0 * half_w).max(2.0 * half_h),
            Profile2d::Circle { radius } => 2.0 * radius,
        };
        prof.max(self.length)
    }

    /// The eight (box) or general corner points of this solid in 3-D, used to
    /// bound the axial interval along an arbitrary prism direction.
    fn corner_points(&self) -> Vec<Point3> {
        let outline = match self.profile.outline() {
            Ok(o) => o,
            Err(_) => return Vec::new(),
        };
        let ring: Vec<[f64; 2]> = match outline {
            ProfileGeom::Polygon(r) => r,
            // For a circle the bounding interval along `d == axis` only needs the
            // two cap centres ± radius along d, but circular operands never reach
            // the interval computation along a non-axis direction (rejected
            // earlier). Use the axis endpoints; the radius extent is handled by
            // the axis-aligned interval directly.
            ProfileGeom::Circle { .. } => Vec::new(),
        };
        let mut pts = Vec::with_capacity(ring.len() * 2);
        for p in &ring {
            let base = self.origin + self.u * p[0] + self.v * p[1];
            pts.push(base);
            pts.push(base + self.axis * self.length);
        }
        pts
    }

    /// Reprofile this leaf onto `frame`: produce its 2-D region and axial
    /// interval along the frame direction `d`.
    ///
    /// Requires `d` to be one of this leaf's prismatic directions (the caller
    /// guarantees it via [`common_direction`]).
    fn reprofile(&self, frame: &Frame, tol: &Tol) -> Result<PrismOperand, PrismError> {
        // Axial interval: project every corner onto d and take min/max. For a
        // circular section viewed along its own axis the corner list is empty, so
        // fall back to the two cap centres.
        let corners = self.corner_points();
        let (t0, t1) = if corners.is_empty() {
            let a = frame.axial(self.origin);
            let b = frame.axial(self.origin + self.axis * self.length);
            (a.min(b), a.max(b))
        } else {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for c in &corners {
                let t = frame.axial(*c);
                lo = lo.min(t);
                hi = hi.max(t);
            }
            (lo, hi)
        };

        // Cross-section region perpendicular to d. Take the face of the solid at
        // the t0 end and project it into the frame: the prism cross-section is
        // invariant along d, so the t0 face is the section. We build it as the
        // convex hull-free explicit ring of the solid face whose outward normal is
        // most anti-aligned with d.
        let region = self.section_region(frame, tol)?;
        Ok(PrismOperand { region, t0, t1 })
    }

    /// The cross-section region of this solid perpendicular to `frame.d`,
    /// expressed in the frame's 2-D coordinates.
    fn section_region(&self, frame: &Frame, _tol: &Tol) -> Result<Region, PrismError> {
        match self.profile.outline() {
            Ok(ProfileGeom::Polygon(ring)) => {
                // The solid's eight (box) / 2N (prism) corners; project them and
                // take the section as the polygon of the face perpendicular to d.
                // When d == axis, the section is the profile ring lifted by the
                // profile basis. When d is a box axis, the section is one of the
                // box's side faces. In both cases the perpendicular face is the
                // rectangle / polygon spanned by the two non-d axes; we obtain it
                // by projecting the solid's corner set and tracing the extreme
                // rectangle. For the supported operands (boxes and axis-aligned
                // prisms) the section is exactly the set of projected corners with
                // duplicate (collapsed-along-d) points merged, traced as a convex
                // ring. We compute it as the 2-D outline of the projected corners.
                let pts3 = self.section_face_points(frame, &ring);
                let ring2: Vec<[f64; 2]> = pts3.iter().map(|&p| frame.project(p)).collect();
                let contour = Contour::from_points(
                    &ring2
                        .iter()
                        .map(|p| crate::boolean::poly2d::Point2::new(p[0], p[1]))
                        .collect::<Vec<_>>(),
                );
                Ok(Region::new(vec![contour]))
            }
            Ok(ProfileGeom::Circle { .. }) => Err(PrismError::CircularInvolved {
                operand: Operand::A, // overwritten by the caller with the real side
            }),
            Err(e) => Err(PrismError::Poly2(
                crate::boolean::poly2d::Poly2Error::Internal {
                    what: leak_profile_error(e),
                },
            )),
        }
    }

    /// The 3-D points of the section face perpendicular to `frame.d`.
    ///
    /// * If `d` is (anti)parallel to the extrusion axis, the section is the
    ///   profile ring lifted at the `t0` cap — i.e. the profile polygon itself,
    ///   which preserves concavity (the H-section is concave).
    /// * Otherwise `d` is a box axis (only rectangles reach here with a non-axis
    ///   `d`); the section is the box's side face perpendicular to `d`, whose four
    ///   corners are two profile-ring corners at both extrusion ends.
    fn section_face_points(&self, frame: &Frame, ring: &[[f64; 2]]) -> Vec<Point3> {
        let cross = self.axis.cross(frame.d).norm();
        let along_axis = cross <= 1e-9_f64 * self.max_dimension().max(1.0);
        if along_axis {
            // Section == profile ring lifted at the bottom cap.
            ring.iter()
                .map(|&p| self.origin + self.u * p[0] + self.v * p[1])
                .collect()
        } else {
            // d is a profile axis (rectangle only). The section is the box face
            // perpendicular to d: pick the profile edge whose direction is
            // perpendicular to d, and sweep it along the extrusion. Concretely,
            // the rectangle's four section corners are: take the two ring corners
            // that are extreme along the projection that is *not* d, at both
            // extrusion ends. Simpler and robust: the section is the rectangle
            // spanned by the extrusion direction and the in-plane profile axis
            // perpendicular to d. We construct it from the box's eight corners by
            // projecting to the plane perpendicular to d and taking the
            // axis-aligned (in u',v') bounding rectangle, which for an oriented box
            // perpendicular cut is exact.
            let mut corners = Vec::with_capacity(ring.len() * 2);
            for &p in ring {
                let base = self.origin + self.u * p[0] + self.v * p[1];
                corners.push(base);
                corners.push(base + self.axis * self.length);
            }
            // Project to (e1, e2), then trace the convex outline (a rectangle).
            let proj: Vec<[f64; 2]> = corners.iter().map(|&c| frame.project(c)).collect();
            let hull = convex_hull(&proj);
            hull.iter().map(|&xy| frame.lift(xy, 0.0)).collect()
        }
    }
}

/// Convert a profile-construction error into a leaked `&'static str` reason. The
/// strings are a closed set, so leaking is bounded and keeps `Poly2Error`'s
/// `&'static str` contract without storing the full `KernelError`.
fn leak_profile_error(_e: crate::error::KernelError) -> &'static str {
    "profile outline could not be built"
}

/// 2-D convex hull (monotone chain) of a small point set. Used to recover the
/// rectangular section of a box cut perpendicular to a box axis; the input is a
/// handful of points so the `O(n log n)` cost is irrelevant.
fn convex_hull(pts: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut p: Vec<[f64; 2]> = pts.to_vec();
    p.sort_by(|a, b| {
        a[0].partial_cmp(&b[0])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a[1].partial_cmp(&b[1]).unwrap_or(std::cmp::Ordering::Equal))
    });
    p.dedup_by(|a, b| (a[0] - b[0]).abs() < 1e-12 && (a[1] - b[1]).abs() < 1e-12);
    let n = p.len();
    if n < 3 {
        return p;
    }
    let cross = |o: [f64; 2], a: [f64; 2], b: [f64; 2]| {
        (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    };
    let mut lower: Vec<[f64; 2]> = Vec::new();
    for &pt in &p {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], pt) <= 0.0 {
            lower.pop();
        }
        lower.push(pt);
    }
    let mut upper: Vec<[f64; 2]> = Vec::new();
    for &pt in p.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], pt) <= 0.0 {
            upper.pop();
        }
        upper.push(pt);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// Find a common prismatic direction of two leaves, build the shared frame, and
/// reprofile both operands onto it.
///
/// Returns the frame and the two [`PrismOperand`]s, or a [`PrismError`] when no
/// common direction exists or a circular section would land on the 2-D side.
pub(crate) fn detect(
    a: &Leaf,
    b: &Leaf,
    tol: &Tol,
) -> Result<(Frame, PrismOperand, PrismOperand), PrismError> {
    let l_max = a.max_dimension().max(b.max_dimension()).max(1.0);
    let cands_a = a.candidate_directions();
    let cands_b = b.candidate_directions();

    // Search for the first agreeing direction pair (length-based agreement).
    let mut chosen: Option<Vec3> = None;
    'outer: for &da in &cands_a {
        for &db in &cands_b {
            // sinθ · L_max ≤ tol.length, with sinθ = |da × db| for unit vectors.
            if da.cross(db).norm() * l_max <= tol.length {
                chosen = Some(da);
                break 'outer;
            }
        }
    }
    let d = chosen.ok_or(PrismError::NoCommonDirection)?;

    // A circular section is only prismatic along its own axis; if a circle is
    // involved at all, the 2-D side carries an arc → Phase 3c.
    if a.is_circular() {
        return Err(PrismError::CircularInvolved {
            operand: Operand::A,
        });
    }
    if b.is_circular() {
        return Err(PrismError::CircularInvolved {
            operand: Operand::B,
        });
    }

    let frame = Frame::new(d);
    let pa = a
        .reprofile(&frame, tol)
        .map_err(|e| relabel(e, Operand::A))?;
    let pb = b
        .reprofile(&frame, tol)
        .map_err(|e| relabel(e, Operand::B))?;
    Ok((frame, pa, pb))
}

/// Re-tag a [`PrismError::CircularInvolved`] with the correct operand side
/// (`section_region` cannot know which side it is building).
fn relabel(e: PrismError, operand: Operand) -> PrismError {
    match e {
        PrismError::CircularInvolved { .. } => PrismError::CircularInvolved { operand },
        other => other,
    }
}
