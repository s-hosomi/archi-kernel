//! Solid × half-space cut.
//!
//! [`cut`] slices a solid by a cutting plane, keeps one side, and seals the
//! opening with cap faces on the cut plane. The output is a fresh watertight
//! [`Brep`] that passes [`ValidateLevel::Full`](crate::topo::ValidateLevel).
//!
//! # Pipeline (`docs/research/05-boolean.md`, the cut is Imprint → Classify → Cap)
//!
//! Faces are processed independently. All shared geometry (vertices, the new
//! section curves on the cut plane) is interned in coordinate-keyed caches, so
//! when two adjacent faces split the same edge they reach the *same*
//! [`Id<Vertex>`](crate::topo::arena::Id) and the *same*
//! [`CurveId`](crate::geom::CurveId) — which is exactly what makes the result
//! watertight (sibling half-edges pair by shared curve + reversed boundary).
//!
//! ## Vertex / edge classification and splitting
//!
//! Each vertex is classified `Below` / `On` / `Above` once
//! ([`predicate::side_of_plane`](crate::predicate::side_of_plane)). A half-edge
//! whose endpoints straddle the plane is split at its closed-form intersection
//! (line: [`line_plane`]; circle / ellipse: the analytic solution of
//! `A cos t + B sin t = C`, the conic parameterisation substituted into the
//! plane equation — see [`split_conic_param`]). The split point is `On`.
//!
//! ## Face splitting and `On` degeneracy
//!
//! For a straddling face the kept-side boundary is walked from the per-loop
//! point sequence: boundary fragments on the kept side are kept verbatim, and
//! consecutive `On` vertices that bound the *interior* of the kept region are
//! joined by a new **section edge** along the surface ∩ plane curve
//! ([`plane_plane`] / [`plane_cylinder`]). When the cut lands exactly on an
//! existing edge (`On`–`On`, the building common case) that edge is reused as
//! the section edge — no new vertex, no split, no T-junction.
//!
//! ## Coplanar (`On`) faces
//!
//! A face lying wholly in the cut plane is governed by the rule in
//! `DESIGN.md` §4.3 as specialised to a cut: for `KeepSide::Below`, an `On` face
//! whose outward normal agrees with the plane normal (`dot > 0`) is the cap of
//! the kept material and is kept; the opposite-facing one is dropped. `Above` is
//! symmetric. Such a kept `On` face's loops feed the cap pool directly.
//!
//! ## Cap generation
//!
//! The `On` segments form open chains that are stitched end-to-end (by shared
//! `On` vertex) into closed loops. In the cut plane's 2-D frame the loops are
//! nested by exact signed area / containment ([`crate::boolean::support`]),
//! giving cap faces with an outer loop and hole loops — a through-hole solid cut
//! across the hole therefore yields an annulus cap. The cap surface is the cut
//! plane inserted canonically; its outward normal points away from the kept
//! material.

use std::collections::HashMap;

use crate::brep::Brep;
use crate::csg::EvalError;
use crate::geom::{CurveGeom, CurveId, SurfaceGeom, VertexGeom};
use crate::intersect::{line_plane, plane_cylinder, plane_plane, PlaneCylinder, PlanePlane};
use crate::math::{Point3, Vec3};
use crate::predicate::side_of_plane;
use crate::primitives::{plane_basis, Circle3, Cylinder, Ellipse3, Line3, Plane};
use crate::tolerance::{Sign3, Tol};
use crate::topo::arena::Id;
use crate::topo::validate::ValidateLevel;
use crate::topo::{Face, HalfEdge, Loop, Sense, Shell, Solid, Vertex};

use super::support::{key, point_in_polygon, signed_area_2d, CoordKey, PlaneFrame};

/// A section connector leaving a vertex: `(far vertex, far point, near point)`.
type SegLink = (Id<Vertex>, Point3, Point3);
/// Adjacency from a vertex to the section connectors leaving it.
type SegFrom = HashMap<Id<Vertex>, Vec<SegLink>>;
/// A conic's centre and its two in-plane semi-axis vectors.
type ConicAxes = (Point3, Vec3, Vec3);

/// Which side of the cutting plane to keep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum KeepSide {
    /// Keep the side the plane normal points *away* from (signed distance ≤ 0).
    Below,
    /// Keep the side the plane normal points *toward* (signed distance ≥ 0).
    Above,
}

/// The outcome of a cut.
///
/// The cut is total: it always produces a [`Brep`] (possibly empty), together
/// with the ids of the cap faces created on the cut plane. The cap-face list is
/// what the Phase 4 section drawing consumes directly.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum CutResult {
    /// The plane separates the solid; `brep` is the kept side and `caps` are the
    /// new faces sealing the opening.
    Cut {
        /// The kept solid.
        brep: Brep,
        /// The cap faces created on the cut plane.
        caps: Vec<Id<Face>>,
    },
    /// The whole solid lies on the kept side; the input is reproduced unchanged
    /// (no caps).
    AllKept {
        /// A copy of the kept solid.
        brep: Brep,
    },
    /// Nothing of the solid lies on the kept side; the result is empty.
    Empty,
}

impl CutResult {
    /// The resulting B-rep, or an empty B-rep when the cut removed everything.
    pub fn brep(&self) -> Brep {
        match self {
            CutResult::Cut { brep, .. } | CutResult::AllKept { brep } => brep.clone(),
            CutResult::Empty => Brep::new(),
        }
    }

    /// The cap faces created on the cut plane (empty unless the plane cut the
    /// solid).
    pub fn caps(&self) -> &[Id<Face>] {
        match self {
            CutResult::Cut { caps, .. } => caps,
            CutResult::AllKept { .. } | CutResult::Empty => &[],
        }
    }
}

/// Cut `solid` of `brep` by `plane`, keeping the side selected by `keep`.
///
/// The input `brep` is read-only; the result is a brand-new [`Brep`] with its
/// own geometry store (`DESIGN.md` §5.2). The opening on the cut plane is sealed
/// with cap faces, so the result is always watertight, and it is validated at
/// [`ValidateLevel::Full`].
///
/// # Errors
///
/// Returns [`EvalError::InvalidResult`] if the cut produced a structurally
/// invalid B-rep (a self-check that should never fire for valid input — it is
/// the cut's fail-safe, `DESIGN.md` §7).
pub fn cut(
    brep: &Brep,
    solid: Id<Solid>,
    plane: &Plane,
    keep: KeepSide,
    tol: &Tol,
) -> Result<CutResult, EvalError> {
    let mut cutter = Cutter::new(brep, *plane, keep, *tol);
    cutter.run(solid)
}

/// The whole cut state: the output brep plus the de-dup caches.
struct Cutter<'a> {
    input: &'a Brep,
    plane: Plane,
    keep: KeepSide,
    tol: Tol,
    out: Brep,
    /// New / surviving vertices, keyed on quantised coordinate.
    vert_by_key: HashMap<CoordKey, Id<Vertex>>,
    /// Resolved coordinates for each vertex key.
    coord_by_key: HashMap<CoordKey, Point3>,
    /// Surviving original line curves, keyed on quantised coordinate so the two
    /// faces sharing an edge resolve to one [`CurveId`].
    line_by_key: HashMap<(CoordKey, CoordKey), CurveId>,
    /// Surviving original conic curves, keyed on the input [`CurveId`] (their
    /// geometry is unchanged by trimming).
    conic_by_src: HashMap<CurveId, CurveId>,
    /// New section line curves on the cut plane, keyed on endpoint coordinates.
    section_line_by_key: HashMap<(CoordKey, CoordKey), CurveId>,
    /// Collected section edges (on the cut plane) for cap assembly: endpoints +
    /// curve + the boundary as seen from the *kept material* side (so the cap's
    /// half-edge is the reversed sibling).
    section_edges: Vec<SectionEdge>,
    /// Collected section arcs (on the cut plane) for a disk / annular cap when
    /// the section curve is a circle or ellipse (cylinder cuts).
    section_arcs: Vec<SectionArc>,
    /// The shared output curve for the section conic (circle / ellipse), so the
    /// half-cylinder faces and the cap reference one curve.
    section_conic_curve: Option<CurveId>,
    /// Faces accumulated into the output shell.
    faces: Vec<Id<Face>>,
}

/// The conic the cut plane makes with a cylinder (the cap-edge geometry).
#[derive(Debug, Clone, Copy)]
enum SectionConic {
    Circle(Circle3),
    Ellipse(Ellipse3),
}

impl SectionConic {
    /// The angle parameter of the conic at a 3-D point on it.
    fn param_at(&self, p: Point3) -> f64 {
        match self {
            SectionConic::Circle(c) => {
                let (u, v) = plane_basis(c.normal());
                let d = p - c.center();
                d.dot(v).atan2(d.dot(u))
            }
            SectionConic::Ellipse(e) => {
                // X(t) = c + a cos t · û + b sin t · v̂ ⇒ recover t from the
                // projections onto û and v̂ scaled by the semi-axes.
                let u = e.major_dir().as_vec();
                let v = e.normal().cross(e.major_dir().as_vec());
                let d = p - e.center();
                let cx = d.dot(u) / e.semi_major();
                let cy = d.dot(v) / e.semi_minor();
                cy.atan2(cx)
            }
        }
    }
}

/// A section arc on the cut plane, recorded for a disk / annular cap.
#[derive(Debug, Clone, Copy)]
struct SectionArc {
    /// Wall-side start vertex. Retained for the directed-edge model.
    #[allow(dead_code)]
    a: Id<Vertex>,
    /// Wall-side end vertex. The cap arc starts here and reverses the boundary.
    b: Id<Vertex>,
    /// The shared section conic curve.
    curve: CurveId,
    /// Angular boundary `[ta, tb]` on the wall side; the cap arc uses `[tb, ta]`.
    ta: f64,
    tb: f64,
}

/// A section edge lying in the cut plane, recorded for cap stitching.
///
/// The wall-side and cap-side half-edges resolve to the same section
/// [`CurveId`](crate::geom::CurveId) through [`Cutter::section_line`] (keyed on
/// the endpoint coordinates), so only the endpoints are stored here.
#[derive(Debug, Clone, Copy)]
struct SectionEdge {
    /// Wall-side start vertex (at `pa`). Retained for the directed-edge model
    /// even though cap stitching keys on coordinates.
    #[allow(dead_code)]
    a: Id<Vertex>,
    /// Wall-side end vertex (at `pb`). The cap half-edge starts here and runs
    /// back to `pa`, the reversed sibling of the wall edge.
    b: Id<Vertex>,
    pa: Point3,
    pb: Point3,
}

impl<'a> Cutter<'a> {
    fn new(input: &'a Brep, plane: Plane, keep: KeepSide, tol: Tol) -> Self {
        Self {
            input,
            plane,
            keep,
            tol,
            out: Brep::new(),
            vert_by_key: HashMap::new(),
            coord_by_key: HashMap::new(),
            line_by_key: HashMap::new(),
            conic_by_src: HashMap::new(),
            section_line_by_key: HashMap::new(),
            section_edges: Vec::new(),
            section_arcs: Vec::new(),
            section_conic_curve: None,
            faces: Vec::new(),
        }
    }

    fn run(&mut self, solid_id: Id<Solid>) -> Result<CutResult, EvalError> {
        let Some(solid) = self.input.topo.solids.get(solid_id) else {
            return Ok(CutResult::Empty);
        };

        // Classify all reachable vertices once.
        let class = self.classify_vertices(solid);

        // Fast paths: entirely on one side.
        let (any_kept, any_dropped, any_on) = side_summary(&class, self.keep);
        if !any_kept && !any_on {
            return Ok(CutResult::Empty);
        }
        if !any_dropped && !any_on {
            // Whole solid kept by vertex classes. A cylindrical face can still be
            // crossed between its vertices (a chord / oblique cut), so only take
            // the AllKept fast path when no cylinder face is actually crossed.
            if !self.any_cylinder_face_crossed(solid, &class) {
                return Ok(CutResult::AllKept {
                    brep: self.clone_input_solid(solid_id),
                });
            }
        }

        // Process each face of each shell.
        let shells: Vec<Id<Shell>> = solid.shells.clone();
        for shell_id in shells {
            let Some(shell) = self.input.topo.shells.get(shell_id) else {
                continue;
            };
            let face_ids: Vec<Id<Face>> = shell.faces.clone();
            for face_id in face_ids {
                self.process_face(face_id, &class);
            }
        }

        // Build caps from the collected section edges.
        let caps = self.build_caps();

        // If after processing there are no faces, nothing survived.
        if self.faces.is_empty() {
            return Ok(CutResult::Empty);
        }

        // Assemble the shell + solid.
        let shell = self.out.topo.add_shell(Shell {
            faces: self.faces.clone(),
        });
        let solid = self.out.topo.add_solid(Solid {
            shells: vec![shell],
        });
        self.out.solids = vec![solid];

        let brep = std::mem::take(&mut self.out);
        brep.validate(&self.tol, ValidateLevel::Full)
            .map_err(EvalError::InvalidResult)?;

        Ok(CutResult::Cut { brep, caps })
    }

    // ── classification ──────────────────────────────────────────────────────

    /// Classify every vertex reachable from `solid` against the cutting plane.
    fn classify_vertices(&self, solid: &Solid) -> HashMap<Id<Vertex>, Sign3> {
        let mut class = HashMap::new();
        for &shell_id in &solid.shells {
            let Some(shell) = self.input.topo.shells.get(shell_id) else {
                continue;
            };
            for &face_id in &shell.faces {
                let Some(face) = self.input.topo.faces.get(face_id) else {
                    continue;
                };
                for loop_id in self.face_loops(face) {
                    let Some(lp) = self.input.topo.loops.get(loop_id) else {
                        continue;
                    };
                    for &he_id in &lp.half_edges {
                        let Some(he) = self.input.topo.half_edges.get(he_id) else {
                            continue;
                        };
                        if class.contains_key(&he.start) {
                            continue;
                        }
                        let Some(vert) = self.input.topo.vertices.get(he.start) else {
                            continue;
                        };
                        let Some(g) = self.input.geom.point(vert.point) else {
                            continue;
                        };
                        let sign = side_of_plane(&self.plane, g, &self.input.geom, &self.tol);
                        class.insert(he.start, sign);
                    }
                }
            }
        }
        class
    }

    fn face_loops(&self, face: &Face) -> Vec<Id<Loop>> {
        let mut v = Vec::with_capacity(1 + face.inners.len());
        v.push(face.outer);
        v.extend(face.inners.iter().copied());
        v
    }

    // ── output interning ──────────────────────────────────────────────────

    fn vertex_at(&mut self, p: Point3) -> Id<Vertex> {
        let k = key(p);
        if let Some(&v) = self.vert_by_key.get(&k) {
            return v;
        }
        let pid = self.out.geom.insert_point(VertexGeom::Explicit(p));
        let v = self.out.topo.add_vertex(Vertex { point: pid });
        self.vert_by_key.insert(k, v);
        self.coord_by_key.insert(k, p);
        v
    }

    /// Resolve the coordinate of an *input* vertex.
    fn input_point(&self, v: Id<Vertex>) -> Option<Point3> {
        let vert = self.input.topo.vertices.get(v)?;
        self.input.geom.point(vert.point)?.as_point()
    }

    /// Get-or-create a shared output line curve for a surviving straight edge.
    fn out_line(&mut self, line: Line3, a: Point3, b: Point3) -> CurveId {
        let (ka, kb) = (key(a), key(b));
        let unordered = if ka <= kb { (ka, kb) } else { (kb, ka) };
        if let Some(&c) = self.line_by_key.get(&unordered) {
            return c;
        }
        let cid = self.out.geom.insert_curve(CurveGeom::Line(line));
        self.line_by_key.insert(unordered, cid);
        cid
    }

    /// Get-or-create a shared output curve for a surviving conic edge (circle /
    /// ellipse), keyed on the source curve so both faces share it.
    fn out_conic(&mut self, src: CurveId, geom: CurveGeom) -> CurveId {
        if let Some(&c) = self.conic_by_src.get(&src) {
            return c;
        }
        let cid = self.out.geom.insert_curve(geom);
        self.conic_by_src.insert(src, cid);
        cid
    }

    /// Get-or-create a shared output section line on the cut plane between two
    /// points, returning the shared curve id **and the line as actually stored**
    /// so every caller parameterises boundaries against the identical origin /
    /// direction. This is what makes the wall-side and cap-side half-edges along
    /// the same section edge exact reversed siblings.
    fn section_line(&mut self, a: Point3, b: Point3) -> Option<(CurveId, Line3)> {
        let (ka, kb) = (key(a), key(b));
        let unordered = if ka <= kb { (ka, kb) } else { (kb, ka) };
        if let Some(&cid) = self.section_line_by_key.get(&unordered) {
            if let Some(CurveGeom::Line(l)) = self.out.geom.curve(cid) {
                return Some((cid, *l));
            }
        }
        // Create the line from the canonical (lexicographically-smaller) endpoint
        // so the parameterisation is independent of call order.
        let (o, other) = if ka <= kb { (a, b) } else { (b, a) };
        let line = Line3::new(o, other - o).ok()?;
        let cid = self.out.geom.insert_curve(CurveGeom::Line(line));
        self.section_line_by_key.insert(unordered, cid);
        Some((cid, line))
    }

    /// Parameter of a point on a stored line.
    #[inline]
    fn line_param(line: &Line3, p: Point3) -> f64 {
        (p - line.origin()).dot(line.dir().as_vec())
    }

    fn clone_input_solid(&self, solid_id: Id<Solid>) -> Brep {
        // The whole solid is kept verbatim: copy the input but keep only this
        // solid as the top-level. The geometry store is shared structurally
        // (cloned), satisfying "output is a fresh store".
        let mut b = self.input.clone();
        b.solids = vec![solid_id];
        b
    }

    // ── face processing ───────────────────────────────────────────────────

    fn process_face(&mut self, face_id: Id<Face>, class: &HashMap<Id<Vertex>, Sign3>) {
        let Some(face) = self.input.topo.faces.get(face_id).cloned() else {
            return;
        };
        let surface = self.input.geom.surface(face.surface).copied();
        let Some(surface) = surface else {
            return;
        };

        // Determine the face's status from its vertex classes.
        let loops = self.face_loops(&face);
        let (kept, dropped, on) = self.face_side_summary(&loops, class);

        // A cylindrical face can be crossed by the plane *between* its vertices
        // (the surface bulges past them), so vertex classes alone can miss the
        // crossing. Detect a geometric crossing and force the split path.
        if let SurfaceGeom::Cylinder(cyl) = &surface {
            if self.cylinder_face_crossed(&face, *cyl, &loops, class) {
                self.split_cylinder_face(&face, *cyl, &loops, class);
                return;
            }
        }

        // Coplanar (On) face: governed by the coplanar rule.
        if !kept && !dropped && on {
            self.process_coplanar_face(&face, &surface, &loops);
            return;
        }

        // A planar face with curved (arc) edges can also be crossed between its
        // vertices — e.g. a disk cap whose circular rim bulges past the chord of
        // an axis-parallel cut. Detect an interior arc crossing and split.
        let arc_crossed = self.face_has_arc_crossing(&loops);

        if dropped && !kept && !arc_crossed {
            // Whole face on the dropped side: drop it. Its On boundary edges (if
            // any) are produced by the neighbouring kept face.
            return;
        }

        if kept && !dropped && !arc_crossed {
            // Whole face kept verbatim.
            self.copy_face_verbatim(&face, &surface, &loops);
            return;
        }

        // Straddling (or arc-crossed) face: split it.
        self.split_face(&face, &surface, &loops, class);
    }

    /// `true` if any arc edge of the face is crossed by the cut plane in its
    /// interior (so the face straddles even when its vertices do not).
    fn face_has_arc_crossing(&self, loops: &[Id<Loop>]) -> bool {
        for &loop_id in loops {
            let Some(lp) = self.input.topo.loops.get(loop_id) else {
                continue;
            };
            for &he_id in &lp.half_edges {
                let Some(he) = self.input.topo.half_edges.get(he_id) else {
                    continue;
                };
                if let Some(curve @ (CurveGeom::Circle(_) | CurveGeom::Ellipse(_))) =
                    self.input.geom.curve(he.curve)
                {
                    if !self.conic_interior_crossings(curve, he.boundary).is_empty() {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// `true` if any cylindrical face of the solid is crossed by the cut plane.
    fn any_cylinder_face_crossed(&self, solid: &Solid, class: &HashMap<Id<Vertex>, Sign3>) -> bool {
        for &shell_id in &solid.shells {
            let Some(shell) = self.input.topo.shells.get(shell_id) else {
                continue;
            };
            for &face_id in &shell.faces {
                let Some(face) = self.input.topo.faces.get(face_id) else {
                    continue;
                };
                if let Some(SurfaceGeom::Cylinder(cyl)) = self.input.geom.surface(face.surface) {
                    let loops = self.face_loops(face);
                    if self.cylinder_face_crossed(face, *cyl, &loops, class) {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// `true` if the cut plane crosses a cylindrical face — either its vertices
    /// straddle, or the plane meets the cylinder surface within the face's
    /// angular span even though all vertices are on one side (a chord / oblique
    /// cut whose intersection lies between the seam edges).
    fn cylinder_face_crossed(
        &self,
        _face: &Face,
        cyl: Cylinder,
        loops: &[Id<Loop>],
        class: &HashMap<Id<Vertex>, Sign3>,
    ) -> bool {
        let (kept, dropped, _on) = self.face_side_summary(loops, class);
        if kept && dropped {
            return true; // vertices already straddle
        }
        // Sample the surface across the face's angular span: if the signed
        // distance to the cut plane changes sign, the plane crosses the face.
        let axis = cyl.axis();
        let (u, v) = plane_basis(axis.dir());
        let axis_dir = axis.dir().as_vec();
        let r = cyl.radius();
        // Collect the face's vertex angles and axial heights to bound the span.
        let mut angles: Vec<f64> = Vec::new();
        let mut zs: Vec<f64> = Vec::new();
        for &loop_id in loops {
            let Some(lp) = self.input.topo.loops.get(loop_id) else {
                continue;
            };
            for &he_id in &lp.half_edges {
                let Some(he) = self.input.topo.half_edges.get(he_id) else {
                    continue;
                };
                if let Some(p) = self.input_point(he.start) {
                    let d = p - axis.origin();
                    angles.push(d.dot(v).atan2(d.dot(u)));
                    zs.push(d.dot(axis_dir));
                }
            }
        }
        if angles.is_empty() {
            return false;
        }
        let z_lo = zs.iter().cloned().fold(f64::INFINITY, f64::min);
        let z_hi = zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        // The true angular span comes from the face's arc edges (their boundary
        // is the φ interval); diametrically-opposite vertices alone are
        // ambiguous. Take the union of arc intervals.
        let (mut span_lo, mut span_hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &loop_id in loops {
            let Some(lp) = self.input.topo.loops.get(loop_id) else {
                continue;
            };
            for &he_id in &lp.half_edges {
                let Some(he) = self.input.topo.half_edges.get(he_id) else {
                    continue;
                };
                if let Some(CurveGeom::Circle(_)) = self.input.geom.curve(he.curve) {
                    span_lo = span_lo.min(he.boundary[0]).min(he.boundary[1]);
                    span_hi = span_hi.max(he.boundary[0]).max(he.boundary[1]);
                }
            }
        }
        if !span_lo.is_finite() || !span_hi.is_finite() {
            // No arc edges: fall back to the vertex-angle span.
            let mut a: Vec<f64> = angles.clone();
            a.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            span_lo = a[0];
            span_hi = a[a.len() - 1];
        }
        let mut seen_pos = false;
        let mut seen_neg = false;
        let steps = 24;
        for i in 0..=steps {
            let phi = span_lo + (span_hi - span_lo) * (i as f64) / (steps as f64);
            for &z in &[z_lo, 0.5 * (z_lo + z_hi), z_hi] {
                let p = axis.origin() + u * (r * phi.cos()) + v * (r * phi.sin()) + axis_dir * z;
                let sd = self.plane.signed_distance(p);
                if sd > self.tol.length {
                    seen_pos = true;
                } else if sd < -self.tol.length {
                    seen_neg = true;
                }
            }
        }
        seen_pos && seen_neg
    }

    /// Summarise whether a face's vertices include kept / dropped / on classes.
    fn face_side_summary(
        &self,
        loops: &[Id<Loop>],
        class: &HashMap<Id<Vertex>, Sign3>,
    ) -> (bool, bool, bool) {
        let (mut kept, mut dropped, mut on) = (false, false, false);
        for &loop_id in loops {
            let Some(lp) = self.input.topo.loops.get(loop_id) else {
                continue;
            };
            for &he_id in &lp.half_edges {
                let Some(he) = self.input.topo.half_edges.get(he_id) else {
                    continue;
                };
                match class.get(&he.start).copied().unwrap_or(Sign3::On) {
                    Sign3::On => on = true,
                    s if is_strictly_kept(s, self.keep) => kept = true,
                    _ => dropped = true,
                }
            }
        }
        (kept, dropped, on)
    }

    /// Copy a face (all loops) verbatim into the output, sharing vertices/curves.
    fn copy_face_verbatim(&mut self, face: &Face, surface: &SurfaceGeom, loops: &[Id<Loop>]) {
        let (surf_id, sense) = self.intern_surface(surface, face.sense);
        let mut out_loops: Vec<Id<Loop>> = Vec::with_capacity(loops.len());
        for &loop_id in loops {
            if let Some(out_loop) = self.copy_loop_verbatim(loop_id) {
                out_loops.push(out_loop);
            }
        }
        if out_loops.is_empty() {
            return;
        }
        let outer = out_loops[0];
        let inners = out_loops[1..].to_vec();
        let f = self.out.topo.add_face(Face {
            surface: surf_id,
            sense,
            outer,
            inners,
        });
        self.faces.push(f);
    }

    /// Copy a loop verbatim, interning vertices and curves into the output.
    fn copy_loop_verbatim(&mut self, loop_id: Id<Loop>) -> Option<Id<Loop>> {
        let lp = self.input.topo.loops.get(loop_id)?.clone();
        let mut hes = Vec::with_capacity(lp.half_edges.len());
        for &he_id in &lp.half_edges {
            let he = *self.input.topo.half_edges.get(he_id)?;
            let start_pt = self.input_point(he.start)?;
            let curve = *self.input.geom.curve(he.curve)?;
            // The end vertex is the start of the next half-edge; evaluate it via
            // the curve at boundary[1] so we can intern both endpoints.
            let end_pt = curve.point_at(he.boundary[1]);
            let _start_v = self.vertex_at(start_pt);
            let _end_v = self.vertex_at(end_pt);
            let out_curve = match curve {
                CurveGeom::Line(l) => self.out_line(l, start_pt, end_pt),
                conic => self.out_conic(he.curve, conic),
            };
            let start = self.vertex_at(start_pt);
            hes.push(self.out.topo.add_half_edge(HalfEdge {
                start,
                curve: out_curve,
                boundary: he.boundary,
            }));
        }
        Some(self.out.topo.add_loop(Loop { half_edges: hes }))
    }

    /// Intern a surface into the output store, returning its id and the adjusted
    /// face sense (folding any canonicalisation flip).
    fn intern_surface(&mut self, surface: &SurfaceGeom, sense: Sense) -> (SurfaceId, Sense) {
        match surface {
            SurfaceGeom::Plane(plane) => {
                let (id, flipped) = self.out.geom.insert_plane(*plane, &self.tol);
                let sense = fold_flip(sense, flipped);
                (id, sense)
            }
            SurfaceGeom::Cylinder(cyl) => {
                let id = self.out.geom.insert_surface(SurfaceGeom::Cylinder(*cyl));
                (id, sense)
            }
        }
    }

    // ── coplanar face rule ────────────────────────────────────────────────

    /// Handle a face lying wholly in the cut plane (`DESIGN.md` §4.3 for a cut).
    ///
    /// For `KeepSide::Below` the kept material is on the `−normal` side, so the
    /// cap is the face whose outward normal points along `+normal`: such a face
    /// is the lid of the kept material and is kept (and feeds the cap). The
    /// opposite-facing coplanar face is dropped. `Above` is symmetric.
    fn process_coplanar_face(&mut self, face: &Face, surface: &SurfaceGeom, loops: &[Id<Loop>]) {
        let SurfaceGeom::Plane(face_plane) = surface else {
            // Only planar faces can be coplanar with the cut plane.
            return;
        };
        let face_normal = match face.sense {
            Sense::Same => face_plane.normal().as_vec(),
            Sense::Reversed => -face_plane.normal().as_vec(),
        };
        let plane_n = self.plane.normal().as_vec();
        let dot = face_normal.dot(plane_n);
        // The cap lid faces away from kept material:
        //   Below keeps −n side ⇒ lid normal = +n ⇒ dot > 0 kept.
        //   Above keeps +n side ⇒ lid normal = −n ⇒ dot < 0 kept.
        let keep_this = match self.keep {
            KeepSide::Below => dot > self.tol.angular,
            KeepSide::Above => dot < -self.tol.angular,
        };
        if !keep_this {
            return;
        }
        // This coplanar face is the lid of the kept material and is already part
        // of the watertight boundary — keep it verbatim. There is no separate
        // opening to cap (the lid *is* the cap), so it is not routed to the cap
        // pool; doing so would duplicate its boundary edges.
        self.copy_face_verbatim(face, surface, loops);
    }

    // ── straddling face split (planar) ────────────────────────────────────

    fn split_face(
        &mut self,
        face: &Face,
        surface: &SurfaceGeom,
        loops: &[Id<Loop>],
        class: &HashMap<Id<Vertex>, Sign3>,
    ) {
        match surface {
            SurfaceGeom::Plane(face_plane) => {
                self.split_planar_face(face, *face_plane, loops, class)
            }
            SurfaceGeom::Cylinder(cyl) => self.split_cylinder_face(face, *cyl, loops, class),
        }
    }

    /// Split a straddling planar face by the cut plane.
    fn split_planar_face(
        &mut self,
        face: &Face,
        face_plane: Plane,
        loops: &[Id<Loop>],
        class: &HashMap<Id<Vertex>, Sign3>,
    ) {
        // The section curve on this face is face_plane ∩ cut_plane (a line).
        let section_line = match plane_plane(&face_plane, &self.plane, &self.tol) {
            PlanePlane::Line(l) => l,
            // Parallel/Coincident shouldn't reach here for a straddling face.
            _ => return,
        };

        // 2-D frame on the *face* plane for the polygon walk / containment.
        let frame = PlaneFrame::new(&face_plane);

        // Build augmented boundary fragments per loop (with split points), and
        // collect the full 2-D outline (kept+dropped) for containment tests.
        let mut kept_fragments: Vec<Vec<KeptNode>> = Vec::new();
        let mut outline_2d: Vec<Vec<[f64; 2]>> = Vec::new();
        let mut portal_vertices: Vec<(Id<Vertex>, Point3)> = Vec::new();

        for &loop_id in loops {
            let aug = self.augment_loop(loop_id, class);
            if aug.is_empty() {
                continue;
            }
            // Record the full 2-D outline of this loop for point-in-face tests.
            outline_2d.push(aug.iter().map(|n| frame.project(n.point)).collect());

            // Extract kept boundary fragments (maximal runs of kept-or-on edges
            // whose midpoint is kept), and the On portal vertices.
            let frags = extract_kept_fragments(&aug, self.keep);
            for f in frags {
                if let (Some(first), Some(last)) = (f.first(), f.last()) {
                    if first.sign == Sign3::On {
                        portal_vertices.push((first.vertex_out, first.point));
                    }
                    if last.sign == Sign3::On {
                        portal_vertices.push((last.vertex_out, last.point));
                    }
                }
                kept_fragments.push(f);
            }
        }

        if kept_fragments.is_empty() {
            return;
        }

        // Pair portals along the section line and decide which segments lie
        // inside the face material (even-odd along the line).
        let section_edges =
            self.section_segments_for_face(&section_line, &portal_vertices, &frame, &outline_2d);

        // Build the kept loops by walking fragments + section edges.
        let (surf_id, sense) = self.intern_surface_plane(face_plane, face.sense);
        let out_loops = self.assemble_kept_loops(&kept_fragments, &section_edges, &frame);
        if out_loops.is_empty() {
            return;
        }
        // Nest the loops (outer vs holes) by 2-D signed area / containment.
        let nested = nest_loops(&out_loops);
        for group in nested {
            let outer = group.outer.loop_id;
            let inners: Vec<Id<Loop>> = group.inners.iter().map(|l| l.loop_id).collect();
            let f = self.out.topo.add_face(Face {
                surface: surf_id,
                sense,
                outer,
                inners,
            });
            self.faces.push(f);
        }
    }

    fn intern_surface_plane(&mut self, plane: Plane, sense: Sense) -> (SurfaceId, Sense) {
        let (id, flipped) = self.out.geom.insert_plane(plane, &self.tol);
        (id, fold_flip(sense, flipped))
    }

    /// Augment one loop into a closed sequence of boundary nodes, inserting a
    /// split node wherever an edge straddles the cut plane. Each node carries an
    /// *output* vertex id and the original curve so kept fragments can rebuild
    /// half-edges. Returns the nodes in loop order.
    fn augment_loop(
        &mut self,
        loop_id: Id<Loop>,
        class: &HashMap<Id<Vertex>, Sign3>,
    ) -> Vec<BoundaryNode> {
        let Some(lp) = self.input.topo.loops.get(loop_id).cloned() else {
            return Vec::new();
        };
        let n = lp.half_edges.len();
        let mut nodes: Vec<BoundaryNode> = Vec::new();
        for i in 0..n {
            let he_id = lp.half_edges[i];
            let next_id = lp.half_edges[(i + 1) % n];
            let (Some(he), Some(next)) = (
                self.input.topo.half_edges.get(he_id).copied(),
                self.input.topo.half_edges.get(next_id).copied(),
            ) else {
                continue;
            };
            let Some(curve) = self.input.geom.curve(he.curve).copied() else {
                continue;
            };
            let start_pt = self
                .input_point(he.start)
                .unwrap_or_else(|| curve.point_at(he.boundary[0]));
            let end_pt = self
                .input_point(next.start)
                .unwrap_or_else(|| curve.point_at(he.boundary[1]));
            let s_start = class.get(&he.start).copied().unwrap_or(Sign3::On);
            let s_end = class.get(&next.start).copied().unwrap_or(Sign3::On);

            let out_start = self.vertex_at(start_pt);
            // Append the start node.
            nodes.push(BoundaryNode {
                vertex_out: out_start,
                point: start_pt,
                sign: s_start,
                // Edge geometry leaving this node:
                edge_curve: he.curve,
                edge_geom: curve,
                edge_b0: he.boundary[0],
                edge_b1: he.boundary[1],
                edge_kept: None,
                split: false,
            });

            // Collect interior crossing parameters on this edge. For a straight
            // edge that straddles there is one; for a conic the plane may cut it
            // 0, 1 or 2 times — including when both endpoints are on one side but
            // the arc bulges across (a chord cut of a cylinder rim).
            let crossings: Vec<(Point3, f64)> = match curve {
                CurveGeom::Line(_) => {
                    if straddles(s_start, s_end) {
                        self.split_point_on_edge(&curve, he.boundary, start_pt, end_pt)
                            .into_iter()
                            .collect()
                    } else {
                        Vec::new()
                    }
                }
                _ => self.conic_interior_crossings(&curve, he.boundary),
            };

            // Insert a split node at each crossing, trimming the running edge.
            let mut prev_b0 = he.boundary[0];
            let crossing_count = crossings.len();
            for (mid_pt, t_mid) in crossings {
                let out_mid = self.vertex_at(mid_pt);
                // Trim the previous node's edge to end at this split.
                if let Some(last) = nodes.last_mut() {
                    last.edge_b1 = t_mid;
                }
                nodes.push(BoundaryNode {
                    vertex_out: out_mid,
                    point: mid_pt,
                    sign: Sign3::On,
                    edge_curve: he.curve,
                    edge_geom: curve,
                    edge_b0: t_mid,
                    edge_b1: he.boundary[1],
                    edge_kept: None,
                    split: true,
                });
                prev_b0 = t_mid;
            }
            let _ = (prev_b0, crossing_count);
        }
        // Classify each node's leaving edge by its midpoint side.
        for node in nodes.iter_mut() {
            let mid_param = 0.5 * (node.edge_b0 + node.edge_b1);
            let mid_pt = node.edge_geom.point_at(mid_param);
            let sd = self.plane.signed_distance(mid_pt);
            let sign = self.tol.classify_length(sd);
            node.edge_kept = Some(match sign {
                Sign3::On => true, // an On edge bounds the cap; treat as kept.
                s => !is_dropped(s, self.keep),
            });
        }
        nodes
    }

    /// All interior crossings of a conic edge with the cut plane, as
    /// `(point, param)` sorted along the edge's directed boundary interval.
    fn conic_interior_crossings(
        &self,
        curve: &CurveGeom,
        boundary: [f64; 2],
    ) -> Vec<(Point3, f64)> {
        let (axes, point_at): (ConicAxes, Box<dyn Fn(f64) -> Point3 + '_>) = match curve {
            CurveGeom::Circle(c) => {
                let (u, v) = plane_basis(c.normal());
                (
                    (c.center(), u * c.radius(), v * c.radius()),
                    Box::new(move |t| c.point_at(t)),
                )
            }
            CurveGeom::Ellipse(e) => {
                let u = e.major_dir().as_vec();
                let v = e.normal().cross(e.major_dir().as_vec());
                (
                    (e.center(), u * e.semi_major(), v * e.semi_minor()),
                    Box::new(move |t| e.point_at(t)),
                )
            }
            CurveGeom::Line(_) => return Vec::new(),
        };
        let (c, p_vec, q_vec) = axes;
        let roots = split_conic_param(c, p_vec, q_vec, &self.plane);
        let (lo, hi) = (boundary[0], boundary[1]);
        let two_pi = std::f64::consts::TAU;
        let mut hits: Vec<(Point3, f64)> = Vec::new();
        for root in roots {
            let t = align_into_interval(root, lo, hi, two_pi);
            let inside = if lo <= hi {
                t > lo + 1e-9_f64 && t < hi - 1e-9_f64
            } else {
                t < lo - 1e-9_f64 && t > hi + 1e-9_f64
            };
            if inside {
                let pt = point_at(t);
                if self.plane.signed_distance(pt).abs() <= self.tol.length * 10.0 {
                    hits.push((pt, t));
                }
            }
        }
        // Sort along the directed interval.
        if lo <= hi {
            hits.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        } else {
            hits.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        }
        hits
    }

    /// Compute the intersection of an edge curve with the cut plane within the
    /// half-edge's boundary interval. Returns `(point, param)`.
    fn split_point_on_edge(
        &self,
        curve: &CurveGeom,
        boundary: [f64; 2],
        start_pt: Point3,
        end_pt: Point3,
    ) -> Option<(Point3, f64)> {
        match curve {
            CurveGeom::Line(line) => {
                let p = line_plane(line, &self.plane, &self.tol).or_else(|| {
                    // Defensive fallback for a straddling segment.
                    let n = self.plane.normal().as_vec();
                    let denom = n.dot(line.dir().as_vec());
                    if denom.abs() < 1e-300_f64 {
                        return None;
                    }
                    let t = n.dot(self.plane.point() - line.origin()) / denom;
                    Some(line.point_at(t))
                })?;
                let t = (p - line.origin()).dot(line.dir().as_vec());
                Some((p, t))
            }
            CurveGeom::Circle(circle) => self.split_conic_within(
                boundary,
                start_pt,
                end_pt,
                |t| circle.point_at(t),
                || {
                    let (u, v) = plane_basis(circle.normal());
                    (circle.center(), u * circle.radius(), v * circle.radius())
                },
            ),
            CurveGeom::Ellipse(ellipse) => self.split_conic_within(
                boundary,
                start_pt,
                end_pt,
                |t| ellipse.point_at(t),
                || {
                    let u = ellipse.major_dir().as_vec();
                    let v = ellipse.normal().cross(ellipse.major_dir().as_vec());
                    (
                        ellipse.center(),
                        u * ellipse.semi_major(),
                        v * ellipse.semi_minor(),
                    )
                },
            ),
        }
    }

    /// Shared conic-split helper: solve the plane equation for the angle, pick
    /// the root lying strictly inside `[boundary]`, and return its point/param.
    fn split_conic_within(
        &self,
        boundary: [f64; 2],
        _start_pt: Point3,
        _end_pt: Point3,
        point_at: impl Fn(f64) -> Point3,
        axes: impl Fn() -> (Point3, Vec3, Vec3),
    ) -> Option<(Point3, f64)> {
        let (c, p_vec, q_vec) = axes();
        let roots = split_conic_param(c, p_vec, q_vec, &self.plane);
        // Map each root into the boundary branch and keep interior ones.
        let (lo, hi) = (boundary[0], boundary[1]);
        let two_pi = std::f64::consts::TAU;
        let mut best: Option<(Point3, f64)> = None;
        for root in roots {
            // Bring `root` into the directed interval [lo, hi].
            let t = align_into_interval(root, lo, hi, two_pi);
            let inside = if lo <= hi {
                t > lo + 1e-9_f64 && t < hi - 1e-9_f64
            } else {
                t < lo - 1e-9_f64 && t > hi + 1e-9_f64
            };
            if inside {
                let pt = point_at(t);
                // Confirm it lies on the plane.
                if self.plane.signed_distance(pt).abs() <= self.tol.length * 10.0 {
                    best = Some((pt, t));
                    break;
                }
            }
        }
        best
    }

    /// Add a section edge (on the cut plane) between two points, interning the
    /// vertices and a shared section line.
    fn add_section_edge(&mut self, pa: Point3, pb: Point3) {
        if (pb - pa).norm() <= self.tol.length {
            return;
        }
        let a = self.vertex_at(pa);
        let b = self.vertex_at(pb);
        // Intern the shared section line so wall- and cap-side half-edges pair.
        let _ = self.section_line(pa, pb);
        self.section_edges.push(SectionEdge { a, b, pa, pb });
    }

    /// Decide the section segments for a face: pair the On portals along the
    /// section line and keep those whose midpoint lies inside the face material.
    fn section_segments_for_face(
        &mut self,
        section_line: &Line3,
        portals: &[(Id<Vertex>, Point3)],
        frame: &PlaneFrame,
        outline_2d: &[Vec<[f64; 2]>],
    ) -> Vec<SectionSeg> {
        // Deduplicate portals by coordinate and sort along the line parameter.
        let mut uniq: HashMap<CoordKey, (Id<Vertex>, Point3, f64)> = HashMap::new();
        for &(v, p) in portals {
            let t = (p - section_line.origin()).dot(section_line.dir().as_vec());
            uniq.entry(key(p)).or_insert((v, p, t));
        }
        let mut sorted: Vec<(Id<Vertex>, Point3, f64)> = uniq.into_values().collect();
        sorted.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        let mut segs = Vec::new();
        // Connect consecutive portal pairs whose midpoint is inside the face.
        for w in sorted.windows(2) {
            let (va, pa, _) = w[0];
            let (vb, pb, _) = w[1];
            let mid = midpoint(pa, pb);
            if self.point_inside_outline(frame.project(mid), outline_2d) {
                segs.push(SectionSeg {
                    a: va,
                    b: vb,
                    pa,
                    pb,
                });
            }
        }
        segs
    }

    /// Point-in-face test against the projected outline loops (outer adds, holes
    /// subtract via the even-odd rule).
    fn point_inside_outline(&self, p: [f64; 2], outline_2d: &[Vec<[f64; 2]>]) -> bool {
        let mut count = 0;
        for poly in outline_2d {
            if point_in_polygon(p, poly) {
                count += 1;
            }
        }
        count % 2 == 1
    }

    /// Walk the kept fragments and section edges into closed output loops.
    fn assemble_kept_loops(
        &mut self,
        fragments: &[Vec<KeptNode>],
        section_segs: &[SectionSeg],
        _frame: &PlaneFrame,
    ) -> Vec<BuiltLoop> {
        // Adjacency by output vertex id: from each vertex, the half-edge chain
        // segments and section segments leaving it.
        // We treat each kept fragment as a directed polyline a→…→b, and each
        // section seg as an undirected connector that we orient as needed.
        // Build a map: start vertex → fragment index, and section connectors.

        // Materialise fragment half-edges (boundary geometry) first.
        struct Frag {
            start: Id<Vertex>,
            end: Id<Vertex>,
            // (start_vertex, curve, boundary) for each half-edge in order.
            hes: Vec<(Id<Vertex>, CurveId, [f64; 2])>,
            points: Vec<Point3>, // for area/winding (start of each segment + final end)
        }

        let mut frags: Vec<Frag> = Vec::new();
        for nodes in fragments {
            if nodes.len() < 2 {
                continue;
            }
            let mut hes = Vec::new();
            let mut points = Vec::new();
            for win in nodes.windows(2) {
                let a = &win[0];
                let out_curve = match a.edge_geom {
                    CurveGeom::Line(l) => self.out_line(l, a.point, win[1].point),
                    conic => self.out_conic(a.edge_curve, conic),
                };
                hes.push((a.vertex_out, out_curve, [a.edge_b0, a.edge_b1]));
                points.push(a.point);
            }
            points.push(nodes.last().unwrap().point);
            frags.push(Frag {
                start: nodes.first().unwrap().vertex_out,
                end: nodes.last().unwrap().vertex_out,
                hes,
                points,
            });
        }

        // Build connectivity: at a vertex, a fragment ends and either another
        // fragment or a section seg continues. We record section connectors as
        // (vertex → (other vertex, pa, pb)).
        let mut seg_from: SegFrom = HashMap::new();
        for s in section_segs {
            seg_from.entry(s.a).or_default().push((s.b, s.pa, s.pb));
            seg_from.entry(s.b).or_default().push((s.a, s.pb, s.pa));
        }
        let mut frag_from: HashMap<Id<Vertex>, Vec<usize>> = HashMap::new();
        for (i, f) in frags.iter().enumerate() {
            frag_from.entry(f.start).or_default().push(i);
        }

        let mut used_frag = vec![false; frags.len()];
        let mut used_seg: std::collections::HashSet<(CoordKey, CoordKey)> =
            std::collections::HashSet::new();
        let mut built = Vec::new();

        for start_idx in 0..frags.len() {
            if used_frag[start_idx] {
                continue;
            }
            // Walk a loop.
            let mut he_ids: Vec<Id<HalfEdge>> = Vec::new();
            let mut poly: Vec<Point3> = Vec::new();
            let loop_start = frags[start_idx].start;
            let mut cur = start_idx;
            let mut guard = 0usize;
            let mut closed = false;
            loop {
                guard += 1;
                if guard > frags.len() * 4 + 8 {
                    break;
                }
                if used_frag[cur] {
                    // Already consumed: loop closed if we are back to start.
                    break;
                }
                used_frag[cur] = true;
                // Emit this fragment's half-edges.
                for &(sv, curve, boundary) in &frags[cur].hes {
                    he_ids.push(self.out.topo.add_half_edge(HalfEdge {
                        start: sv,
                        curve,
                        boundary,
                    }));
                }
                poly.extend(frags[cur].points.iter().copied());
                let frag_end = frags[cur].end;
                if frag_end == loop_start && !he_ids.is_empty() {
                    closed = true;
                    break;
                }
                // Continue: try a section seg from frag_end, then a fragment.
                if let Some(next) = self.continue_via_section(
                    frag_end,
                    &seg_from,
                    &mut used_seg,
                    &mut he_ids,
                    &mut poly,
                ) {
                    if next == loop_start {
                        closed = true;
                        break;
                    }
                    // Find a fragment starting at `next`.
                    if let Some(idxs) = frag_from.get(&next) {
                        if let Some(&ni) = idxs.iter().find(|&&i| !used_frag[i]) {
                            cur = ni;
                            continue;
                        }
                    }
                    // No fragment continues; if next == loop_start we'd be closed.
                    break;
                } else if let Some(idxs) = frag_from.get(&frag_end) {
                    if let Some(&ni) = idxs.iter().find(|&&i| !used_frag[i]) {
                        cur = ni;
                        continue;
                    }
                    break;
                } else {
                    break;
                }
            }

            if closed && he_ids.len() >= 3 {
                let loop_id = self.out.topo.add_loop(Loop { half_edges: he_ids });
                built.push(BuiltLoop {
                    loop_id,
                    points: poly,
                });
            }
        }
        built
    }

    /// Try to continue the walk through a section segment from `from`, emitting
    /// the connector half-edge. Returns the far vertex if one was taken.
    fn continue_via_section(
        &mut self,
        from: Id<Vertex>,
        seg_from: &SegFrom,
        used_seg: &mut std::collections::HashSet<(CoordKey, CoordKey)>,
        he_ids: &mut Vec<Id<HalfEdge>>,
        poly: &mut Vec<Point3>,
    ) -> Option<Id<Vertex>> {
        let candidates = seg_from.get(&from)?;
        for &(other, pa, pb) in candidates {
            let k = (key(pa), key(pb));
            if used_seg.contains(&k) {
                continue;
            }
            used_seg.insert(k);
            // Record the section edge for cap building (kept-material side
            // boundary is from→other, the cap is the reversed sibling).
            self.add_section_edge(pa, pb);
            // Emit the wall-side connector half-edge from `from` to `other`,
            // parameterised against the shared stored section line.
            let (curve, line) = self.section_line(pa, pb)?;
            let ta = Self::line_param(&line, pa);
            let tb = Self::line_param(&line, pb);
            he_ids.push(self.out.topo.add_half_edge(HalfEdge {
                start: from,
                curve,
                boundary: [ta, tb],
            }));
            poly.push(pb);
            return Some(other);
        }
        None
    }

    // ── cylinder face split ───────────────────────────────────────────────

    /// Split a straddling cylindrical (half-cylinder) face by the cut plane.
    ///
    /// The extruder builds a cylinder side as two half-cylinder faces, each a
    /// loop of `bottom arc → seam up → top arc → seam down`. A cut that crosses
    /// the face splits the two seam (vertical) edges, keeps the wanted-side arc,
    /// and joins the two split points with the **section arc** — an arc of the
    /// `plane ∩ cylinder` curve (a circle for a perpendicular cut, an ellipse for
    /// an oblique one, from [`plane_cylinder`]). That section arc is recorded so
    /// the disk cap is rebuilt from the arcs exactly as the extruder builds its
    /// circular caps.
    fn split_cylinder_face(
        &mut self,
        face: &Face,
        cyl: Cylinder,
        loops: &[Id<Loop>],
        class: &HashMap<Id<Vertex>, Sign3>,
    ) {
        // The section curve on the cut plane (the cap-edge geometry).
        let section = match plane_cylinder(&self.plane, &cyl, &self.tol) {
            PlaneCylinder::Circle(c) => SectionConic::Circle(c),
            PlaneCylinder::Ellipse(e) => SectionConic::Ellipse(e),
            // Axis-parallel cuts (ruling lines) and non-crossing cases are
            // handled by the generic straddle walk with straight section edges.
            _ => {
                self.split_cylinder_face_ruling(face, cyl, loops, class);
                return;
            }
        };

        let surf_id = self.out.geom.insert_surface(SurfaceGeom::Cylinder(cyl));
        let sense = face.sense;

        for &loop_id in loops {
            let Some(lp) = self.input.topo.loops.get(loop_id).cloned() else {
                continue;
            };
            // Augmented nodes with split points on straddling seam edges.
            let aug = self.augment_loop(loop_id, class);
            if aug.is_empty() {
                continue;
            }
            // Collect the kept fragments and the two On split points (portals).
            let frags = extract_kept_fragments(&aug, self.keep);
            let mut portals: Vec<Point3> = Vec::new();
            for f in &frags {
                if let (Some(first), Some(last)) = (f.first(), f.last()) {
                    if first.sign == Sign3::On {
                        portals.push(first.point);
                    }
                    if last.sign == Sign3::On {
                        portals.push(last.point);
                    }
                }
            }
            if frags.is_empty() {
                continue;
            }
            // Dedup portals by coordinate.
            let mut seen = std::collections::HashSet::new();
            portals.retain(|&p| seen.insert(key(p)));
            // A half-cylinder face split by a single plane has exactly two
            // portals (the two seam split points). If not, fall back.
            if portals.len() != 2 {
                let _ = lp;
                continue;
            }

            // Build the kept sub-face. The section arc must follow the same
            // angular half this face occupies. Disambiguate the arc interval by a
            // representative interior angle: a non-portal fragment point's angle
            // on the section conic (which shares the cylinder axis, so the angle
            // matches). We pick the [ta, tb] interval (±2π) that contains it.
            let (pa, pb) = (portals[0], portals[1]);
            let rep_angle = frags
                .iter()
                .flat_map(|f| f.iter())
                .find(|n| n.sign != Sign3::On)
                .map(|n| section.param_at(self.plane.project_point(n.point)))
                .unwrap_or_else(|| 0.5 * (section.param_at(pa) + section.param_at(pb)));
            let (ta, tb) =
                arc_interval_containing(section.param_at(pa), section.param_at(pb), rep_angle);

            if let Some(out_loop) = self.assemble_cyl_kept_loop(&frags, pa, pb, &section, ta, tb) {
                let f = self.out.topo.add_face(Face {
                    surface: surf_id,
                    sense,
                    outer: out_loop,
                    inners: Vec::new(),
                });
                self.faces.push(f);
            }
        }
    }

    /// Split a half-cylinder face by an **axis-parallel** cut (a chord cut).
    ///
    /// The cut crosses the face along straight **ruling lines** (vertical, from
    /// [`plane_cylinder`]'s `TwoLines`/`TangentLine`). The face's arc edges are
    /// split where the cut plane meets them; the kept boundary fragments are then
    /// joined by straight ruling-line section edges through the face interior, so
    /// the rectangular bow-cap is built by the same straight-edge cap machinery
    /// as the planar faces (the disk caps split likewise, planar). The face is
    /// kept whole / dropped when it does not straddle.
    fn split_cylinder_face_ruling(
        &mut self,
        face: &Face,
        cyl: Cylinder,
        loops: &[Id<Loop>],
        class: &HashMap<Id<Vertex>, Sign3>,
    ) {
        // Always run the split walk: an axis-parallel cut can cross a half-
        // cylinder face along its bulge even when no vertex straddles, so the
        // vertex summary alone is not enough.
        let surf_id = self.out.geom.insert_surface(SurfaceGeom::Cylinder(cyl));
        let sense = face.sense;
        let mut kept_fragments: Vec<Vec<KeptNode>> = Vec::new();
        let mut portals: Vec<(Id<Vertex>, Point3)> = Vec::new();
        let mut any_dropped = false;

        for &loop_id in loops {
            let aug = self.augment_loop(loop_id, class);
            if aug.is_empty() {
                continue;
            }
            if aug.iter().any(|n| n.edge_kept == Some(false)) {
                any_dropped = true;
            }
            let frags = extract_kept_fragments(&aug, self.keep);
            for f in frags {
                if let (Some(first), Some(last)) = (f.first(), f.last()) {
                    if first.sign == Sign3::On {
                        portals.push((first.vertex_out, first.point));
                    }
                    if last.sign == Sign3::On {
                        portals.push((last.vertex_out, last.point));
                    }
                }
                kept_fragments.push(f);
            }
        }
        if kept_fragments.is_empty() {
            return;
        }

        // Nothing dropped ⇒ the whole face is kept; copy verbatim (the closed
        // no-portal case).
        if !any_dropped {
            let surf = SurfaceGeom::Cylinder(cyl);
            self.copy_face_verbatim(face, &surf, loops);
            return;
        }

        // Connect portals by straight ruling segments. For a single half-cylinder
        // chord cut there are exactly two portals; pair them directly.
        let mut uniq: HashMap<CoordKey, (Id<Vertex>, Point3)> = HashMap::new();
        for &(v, p) in &portals {
            uniq.entry(key(p)).or_insert((v, p));
        }
        let pts: Vec<(Id<Vertex>, Point3)> = uniq.into_values().collect();
        let mut segs: Vec<SectionSeg> = Vec::new();
        if pts.len() == 2 {
            segs.push(SectionSeg {
                a: pts[0].0,
                b: pts[1].0,
                pa: pts[0].1,
                pb: pts[1].1,
            });
        }

        // Reuse the planar loop assembler (it walks fragments + straight section
        // segments and records cap edges).
        let frame = PlaneFrame::new(&self.plane);
        let built = self.assemble_kept_loops(&kept_fragments, &segs, &frame);
        for bl in built {
            let f = self.out.topo.add_face(Face {
                surface: surf_id,
                sense,
                outer: bl.loop_id,
                inners: Vec::new(),
            });
            self.faces.push(f);
        }
    }

    /// Assemble a half-cylinder kept loop from its boundary fragments plus the
    /// section arc joining the two portals.
    fn assemble_cyl_kept_loop(
        &mut self,
        frags: &[Vec<KeptNode>],
        pa: Point3,
        pb: Point3,
        section: &SectionConic,
        ta: f64,
        tb: f64,
    ) -> Option<Id<Loop>> {
        // Emit the boundary fragment half-edges in order, then close with the
        // section arc. The fragments together with the arc form one closed loop.
        let mut hes: Vec<Id<HalfEdge>> = Vec::new();
        // Chain fragments by shared end/start vertex.
        let mut remaining: Vec<&Vec<KeptNode>> = frags.iter().collect();
        // Start with any fragment; greedily append.
        let mut ordered: Vec<&Vec<KeptNode>> = Vec::new();
        if let Some(first) = remaining.pop() {
            ordered.push(first);
        }
        while !remaining.is_empty() {
            let tail = ordered.last().unwrap().last().unwrap().vertex_out;
            if let Some(pos) = remaining
                .iter()
                .position(|f| f.first().map(|n| n.vertex_out) == Some(tail))
            {
                ordered.push(remaining.remove(pos));
            } else {
                // No continuation: append the rest in arbitrary order.
                ordered.push(remaining.remove(0));
            }
        }
        for nodes in &ordered {
            for win in nodes.windows(2) {
                let a = &win[0];
                let out_curve = match a.edge_geom {
                    CurveGeom::Line(l) => self.out_line(l, a.point, win[1].point),
                    conic => self.out_conic(a.edge_curve, conic),
                };
                hes.push(self.out.topo.add_half_edge(HalfEdge {
                    start: a.vertex_out,
                    curve: out_curve,
                    boundary: [a.edge_b0, a.edge_b1],
                }));
            }
        }
        // The fragment chain runs from `frag_start` to `frag_end`. The section
        // arc must connect `frag_end` back to `frag_start`. Determine which
        // portal is which by matching the chain endpoints.
        let chain_start = ordered.first()?.first()?;
        let chain_end = ordered.last()?.last()?;
        let start_pt = chain_start.point;
        let end_pt = chain_end.point;
        // The arc goes from chain_end's point to chain_start's point.
        let (from_pt, to_pt, t_from, t_to) = if key(end_pt) == key(pa) {
            (pa, pb, ta, tb)
        } else if key(end_pt) == key(pb) {
            (pb, pa, tb, ta)
        } else {
            // Fallback orientation.
            (
                end_pt,
                start_pt,
                section.param_at(end_pt),
                section.param_at(start_pt),
            )
        };
        let _ = to_pt;
        // Record the wall-side arc (from_pt → to_pt) for the cap (reversed).
        self.record_section_arc(section, from_pt, to_pt, t_from, t_to);
        // Emit the wall-side arc half-edge.
        let arc_curve = self.section_arc_curve(section);
        let from_v = self.vertex_at(from_pt);
        hes.push(self.out.topo.add_half_edge(HalfEdge {
            start: from_v,
            curve: arc_curve,
            boundary: [t_from, t_to],
        }));
        if hes.len() < 2 {
            return None;
        }
        Some(self.out.topo.add_loop(Loop { half_edges: hes }))
    }

    /// Get-or-create the output curve for the section conic (shared across the
    /// half-cylinder faces and the cap).
    fn section_arc_curve(&mut self, section: &SectionConic) -> CurveId {
        if let Some(c) = self.section_conic_curve {
            return c;
        }
        let geom = match section {
            SectionConic::Circle(c) => CurveGeom::Circle(*c),
            SectionConic::Ellipse(e) => CurveGeom::Ellipse(*e),
        };
        let cid = self.out.geom.insert_curve(geom);
        self.section_conic_curve = Some(cid);
        cid
    }

    /// Record a section arc for the disk cap (wall side `pa → pb`).
    fn record_section_arc(
        &mut self,
        section: &SectionConic,
        pa: Point3,
        pb: Point3,
        ta: f64,
        tb: f64,
    ) {
        let curve = self.section_arc_curve(section);
        let a = self.vertex_at(pa);
        let b = self.vertex_at(pb);
        self.section_arcs.push(SectionArc {
            a,
            b,
            curve,
            ta,
            tb,
        });
    }

    // ── cap building ──────────────────────────────────────────────────────

    /// Stitch the collected section edges into closed loops and emit cap faces.
    ///
    /// Each section edge was recorded as a **directed** wall-side edge `a → b`
    /// (the material boundary). The cap half-edge along it is its reversed
    /// sibling `b → a`. Those reversed edges form one outgoing edge per vertex,
    /// so they chain into directed cycles uniquely — no winding guesswork, and
    /// every cap half-edge is the exact reversed sibling of a wall half-edge.
    fn build_caps(&mut self) -> Vec<Id<Face>> {
        let mut cap_faces = self.build_arc_caps();
        if self.section_edges.is_empty() {
            return cap_faces;
        }
        cap_faces.extend(self.build_straight_caps());
        cap_faces
    }

    /// Build a disk / annular cap from the recorded section arcs (cylinder cut).
    ///
    /// The arcs form a closed conic (one circle / ellipse); the cap is a single
    /// disk face whose boundary is the section arcs, each the reversed sibling of
    /// a half-cylinder wall arc — mirroring the extruder's two-arc disk cap.
    fn build_arc_caps(&mut self) -> Vec<Id<Face>> {
        if self.section_arcs.is_empty() {
            return Vec::new();
        }
        let arcs = self.section_arcs.clone();
        // Cap surface: cut plane inserted canonically.
        let (surf_id, _flip) = self.out.geom.insert_plane(self.plane, &self.tol);
        let desired = self.cap_outward_normal();
        let canon = match self.out.geom.surface(surf_id) {
            Some(SurfaceGeom::Plane(p)) => *p,
            _ => self.plane,
        };
        let canon_normal = canon.normal().as_vec();

        // Emit the cap arcs as reversed siblings of the wall arcs: wall arc was
        // pa→pb on `curve` with boundary [ta, tb]; the cap arc is pb→pa with
        // boundary [tb, ta] on the same curve.
        let mut hes = Vec::with_capacity(arcs.len());
        for arc in &arcs {
            hes.push(self.out.topo.add_half_edge(HalfEdge {
                start: arc.b,
                curve: arc.curve,
                boundary: [arc.tb, arc.ta],
            }));
        }
        let loop_id = self.out.topo.add_loop(Loop { half_edges: hes });

        // Sense: choose so the cap's outward normal equals `desired`. A disk cap
        // bounded by arcs traversed pb→pa winds about the section curve's normal;
        // compare that normal against `desired`.
        let arc_normal = match self.section_arcs.first() {
            Some(a) => match self.out.geom.curve(a.curve) {
                Some(CurveGeom::Circle(c)) => c.normal().as_vec(),
                Some(CurveGeom::Ellipse(e)) => e.normal().as_vec(),
                _ => canon_normal,
            },
            None => canon_normal,
        };
        let _ = arc_normal;
        // The cap outward must equal `desired`; the canonical surface normal is
        // `canon_normal`. Pick sense to match.
        let sense = if canon_normal.dot(desired) > 0.0 {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let f = self.out.topo.add_face(Face {
            surface: surf_id,
            sense,
            outer: loop_id,
            inners: Vec::new(),
        });
        self.faces.push(f);
        vec![f]
    }

    /// Stitch the straight section edges into closed loops and emit cap faces.
    fn build_straight_caps(&mut self) -> Vec<Id<Face>> {
        if self.section_edges.is_empty() {
            return Vec::new();
        }
        // Resolve directed section edges. The cut boundary on the plane is a set
        // of oriented loops: a real cap edge is recorded in exactly one
        // direction (material on one side, cap on the other). If an edge was
        // recorded in *both* directions it is an internal edge shared by two
        // kept sub-faces — it cancels and is not a cap boundary. We tally signed
        // multiplicity per undirected edge and keep only the net direction.
        // Resolve the directed cap boundary. Each recorded section edge `pa→pb`
        // mirrors a wall half-edge that was already emitted with that exact
        // direction, so the cap half-edge along it must run `pb→pa` (the reversed
        // sibling). We must therefore preserve each survivor's **original**
        // direction — never re-canonicalise it. An edge recorded in *both*
        // directions is internal (shared by two kept sub-faces) and cancels.
        let mut signed: HashMap<(CoordKey, CoordKey), i32> = HashMap::new();
        for &e in &self.section_edges {
            let (ka, kb) = (key(e.pa), key(e.pb));
            let (uk, dir) = if ka <= kb {
                ((ka, kb), 1)
            } else {
                ((kb, ka), -1)
            };
            *signed.entry(uk).or_insert(0) += dir;
        }
        let mut emitted: std::collections::HashSet<(CoordKey, CoordKey)> =
            std::collections::HashSet::new();
        let mut edges: Vec<SectionEdge> = Vec::new();
        for &e in &self.section_edges {
            let (ka, kb) = (key(e.pa), key(e.pb));
            let uk = if ka <= kb { (ka, kb) } else { (kb, ka) };
            if signed.get(&uk).copied().unwrap_or(0) == 0 {
                continue; // internal edge: cancels, not a cap boundary.
            }
            // Keep one survivor per directed edge, preserving its direction.
            if emitted.insert((ka, kb)) {
                edges.push(e);
            }
        }

        // Directed cap edges b → a: index by their start vertex key (pb).
        // The cap half-edge starts at vertex `b` (point pb) and ends at `a`
        // (point pa) — the reversed wall edge leaving pb.
        let mut next_cap: HashMap<CoordKey, (Id<Vertex>, Point3, Point3)> = HashMap::new();
        for e in &edges {
            next_cap.insert(key(e.pb), (e.b, e.pb, e.pa));
        }

        // Walk directed cycles.
        let mut visited: std::collections::HashSet<CoordKey> = std::collections::HashSet::new();
        let mut cycles: Vec<Vec<(Id<Vertex>, Point3)>> = Vec::new();
        for e in &edges {
            let start_key = key(e.pb);
            if visited.contains(&start_key) {
                continue;
            }
            let mut ring: Vec<(Id<Vertex>, Point3)> = Vec::new();
            let mut cur = start_key;
            let mut guard = 0usize;
            loop {
                guard += 1;
                if guard > edges.len() * 2 + 8 {
                    break;
                }
                if !visited.insert(cur) {
                    break;
                }
                let Some(&(v, p_from, _p_to)) = next_cap.get(&cur) else {
                    break;
                };
                ring.push((v, p_from));
                cur = key(_p_to);
                if cur == start_key {
                    break;
                }
            }
            if ring.len() >= 3 {
                cycles.push(ring);
            }
        }
        if cycles.is_empty() {
            return Vec::new();
        }

        // Nest the cycles by exact 2-D containment.
        let frame = PlaneFrame::new(&self.plane);
        let cap_cycles: Vec<CapCycle> = cycles
            .into_iter()
            .map(|ring| {
                let proj = ring.iter().map(|&(_, p)| frame.project(p)).collect();
                CapCycle { ring, proj }
            })
            .collect();
        let nested = nest_cap_cycles(&cap_cycles);

        // Cap surface: cut plane inserted canonically; sense so its outward
        // normal faces away from the kept material.
        let (surf_id, _flipped) = self.out.geom.insert_plane(self.plane, &self.tol);
        let desired = self.cap_outward_normal();
        let canon = match self.out.geom.surface(surf_id) {
            Some(SurfaceGeom::Plane(p)) => *p,
            _ => self.plane,
        };
        let canon_normal = canon.normal().as_vec();
        // The cap's outer ring is built from the reversed wall edges; its winding
        // in the canonical (u × v = canon_normal) frame fixes the face's natural
        // outward normal. Choose `sense` so the actual outward equals `desired`.
        let canon_frame = PlaneFrame::new(&canon);

        let mut cap_faces = Vec::new();
        for group in nested {
            // Natural outward of the built outer loop (right-hand rule about its
            // winding in the canonical frame).
            let outer_area = signed_area_2d(
                &group
                    .outer
                    .ring
                    .iter()
                    .map(|&(_, p)| canon_frame.project(p))
                    .collect::<Vec<_>>(),
            );
            let natural_outward = if outer_area >= 0.0 {
                canon_normal
            } else {
                -canon_normal
            };
            let sense = if natural_outward.dot(desired) > 0.0 {
                Sense::Same
            } else {
                Sense::Reversed
            };

            let outer = self.build_cap_loop(&group.outer);
            let inners: Vec<Id<Loop>> = group
                .inners
                .iter()
                .map(|c| self.build_cap_loop(c))
                .collect();
            let f = self.out.topo.add_face(Face {
                surface: surf_id,
                sense,
                outer,
                inners,
            });
            self.faces.push(f);
            cap_faces.push(f);
        }
        cap_faces
    }

    /// The cap's outward normal (away from kept material).
    fn cap_outward_normal(&self) -> Vec3 {
        let n = self.plane.normal().as_vec();
        match self.keep {
            // Kept material on −n ⇒ cap faces +n.
            KeepSide::Below => n,
            // Kept material on +n ⇒ cap faces −n.
            KeepSide::Above => -n,
        }
    }

    /// Build a cap loop's half-edges directly from its directed ring (each edge
    /// is the reversed sibling of a wall half-edge, sharing the section curve).
    fn build_cap_loop(&mut self, cycle: &CapCycle) -> Id<Loop> {
        let ring = &cycle.ring;
        let n = ring.len();
        let mut hes = Vec::with_capacity(n);
        for i in 0..n {
            let (va, pa) = ring[i];
            let (_vb, pb) = ring[(i + 1) % n];
            let Some((curve, line)) = self.section_line(pa, pb) else {
                continue;
            };
            let ta = Self::line_param(&line, pa);
            let tb = Self::line_param(&line, pb);
            hes.push(self.out.topo.add_half_edge(HalfEdge {
                start: va,
                curve,
                boundary: [ta, tb],
            }));
        }
        self.out.topo.add_loop(Loop { half_edges: hes })
    }
}

use crate::geom::SurfaceId;

/// `true` for a sign strictly on the kept side (On excluded).
fn is_strictly_kept(sign: Sign3, keep: KeepSide) -> bool {
    matches!(
        (keep, sign),
        (KeepSide::Below, Sign3::Below) | (KeepSide::Above, Sign3::Above)
    )
}

/// Fold a canonicalisation flip into a face sense.
fn fold_flip(sense: Sense, flipped: bool) -> Sense {
    match (sense, flipped) {
        (Sense::Same, false) | (Sense::Reversed, true) => Sense::Same,
        _ => Sense::Reversed,
    }
}

/// Overall summary across all classified vertices.
fn side_summary(class: &HashMap<Id<Vertex>, Sign3>, keep: KeepSide) -> (bool, bool, bool) {
    let (mut kept, mut dropped, mut on) = (false, false, false);
    for &s in class.values() {
        match s {
            Sign3::On => on = true,
            s if is_strictly_kept(s, keep) => kept = true,
            _ => dropped = true,
        }
    }
    (kept, dropped, on)
}

// ── conic split parameter ─────────────────────────────────────────────────

// ── boundary-walk data ────────────────────────────────────────────────────

/// One node in a face loop's augmented boundary walk.
#[derive(Debug, Clone, Copy)]
struct BoundaryNode {
    vertex_out: Id<Vertex>,
    point: Point3,
    sign: Sign3,
    /// The source curve id of the edge *leaving* this node.
    edge_curve: CurveId,
    /// The source curve geometry of the edge leaving this node.
    edge_geom: CurveGeom,
    /// Boundary interval of the edge leaving this node, possibly trimmed at a
    /// split.
    edge_b0: f64,
    edge_b1: f64,
    /// Whether the edge *leaving* this node is on the kept side, decided by its
    /// midpoint (so a sub-arc with both endpoints `On` is classified correctly).
    /// `None` until filled in by [`Cutter::fill_edge_sides`].
    edge_kept: Option<bool>,
    /// Whether this node was inserted by an edge split.
    #[allow(dead_code)]
    split: bool,
}

/// A section segment connecting two On portals through the face interior.
#[derive(Debug, Clone, Copy)]
struct SectionSeg {
    a: Id<Vertex>,
    b: Id<Vertex>,
    pa: Point3,
    pb: Point3,
}

/// A built loop with its 2-D-projectable point ring (for nesting).
#[derive(Debug, Clone)]
struct BuiltLoop {
    loop_id: Id<Loop>,
    points: Vec<Point3>,
}

/// A nesting group: one outer loop and its hole loops.
struct LoopGroup {
    outer: BuiltLoop,
    inners: Vec<BuiltLoop>,
}

/// `true` if the two endpoint signs are on strictly opposite sides.
fn straddles(a: Sign3, b: Sign3) -> bool {
    matches!(
        (a, b),
        (Sign3::Below, Sign3::Above) | (Sign3::Above, Sign3::Below)
    )
}

/// Pick the directed angular interval between `t0` and `t1` that passes through
/// `rep` (an interior representative angle), returning `[from, to]`.
///
/// Angles are on a circle; the two endpoints split it into two arcs. We choose
/// the arc that contains `rep` and orient the interval so a parameter sweep from
/// `from` to `to` stays on that arc (allowing one endpoint to exceed `2π` so the
/// half-edge boundary is monotone, matching the curve parameterisation).
fn arc_interval_containing(t0: f64, t1: f64, rep: f64) -> (f64, f64) {
    let two_pi = std::f64::consts::TAU;
    let norm = |a: f64| {
        let mut x = a % two_pi;
        if x < 0.0 {
            x += two_pi;
        }
        x
    };
    let (n0, n1, nr) = (norm(t0), norm(t1), norm(rep));
    // Arc going from n0 increasing to n1 (wrapping): does it contain nr?
    let span_inc = norm(n1 - n0);
    let rep_off = norm(nr - n0);
    if rep_off <= span_inc {
        // Increasing arc n0 → n0 + span_inc contains rep.
        (n0, n0 + span_inc)
    } else {
        // The other arc: n1 → n1 + (2π − span_inc) = n0 + 2π.
        (n1, n1 + (two_pi - span_inc))
    }
}

/// Midpoint of two points.
fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
}

/// Bring an angle `root` into the directed interval `[lo, hi]` by adding /
/// subtracting `period`.
fn align_into_interval(root: f64, lo: f64, hi: f64, period: f64) -> f64 {
    let mut t = root;
    if lo <= hi {
        while t < lo {
            t += period;
        }
        while t > hi {
            t -= period;
        }
    } else {
        // Decreasing interval [lo, hi] with lo > hi.
        while t > lo {
            t -= period;
        }
        while t < hi {
            t += period;
        }
    }
    t
}

/// Extract maximal kept boundary fragments from an augmented loop.
///
/// A fragment is a maximal run of consecutive nodes such that each edge between
/// them has its midpoint on the kept side (so the edge is part of the kept
/// region's boundary). Fragments begin and end at On vertices (the cut entry /
/// exit) or wrap the whole loop when nothing is dropped (handled earlier).
fn extract_kept_fragments(nodes: &[BoundaryNode], keep: KeepSide) -> Vec<Vec<KeptNode>> {
    let n = nodes.len();
    if n < 2 {
        return Vec::new();
    }
    // For each directed edge i → i+1, decide if it is a kept boundary edge.
    // An edge is kept if both endpoints are kept-or-on AND not both On-dropped.
    // Use the endpoint signs: keep the edge if neither endpoint is on the
    // dropped side, i.e. both are kept-or-on, AND at least one is strictly kept
    // OR the edge midpoint is on the kept side. For straddle-free augmented
    // loops, an edge between two non-dropped vertices is kept.
    // Use the midpoint-classified `edge_kept` when available (it handles a
    // sub-arc whose endpoints are both `On` but whose interior is dropped, e.g.
    // a chord cut of a cylinder rim); fall back to the endpoint-sign rule.
    let kept_edge: Vec<bool> = (0..n)
        .map(|i| {
            nodes[i].edge_kept.unwrap_or_else(|| {
                let a = nodes[i].sign;
                let b = nodes[(i + 1) % n].sign;
                edge_is_kept(a, b, keep)
            })
        })
        .collect();

    // Find a starting index where the previous edge is NOT kept (a fragment
    // boundary). If all edges are kept, the whole loop is one cycle.
    let all_kept = kept_edge.iter().all(|&k| k);
    if all_kept {
        // Entire loop kept: single closed fragment (start == end).
        let mut frag: Vec<KeptNode> = nodes.iter().map(KeptNode::from).collect();
        // Close it by repeating the first as the end marker.
        frag.push(KeptNode::from(&nodes[0]));
        return vec![frag];
    }

    let start = (0..n).find(|&i| !kept_edge[(i + n - 1) % n] && kept_edge[i]);
    let Some(start) = start else {
        return Vec::new();
    };

    let mut frags: Vec<Vec<KeptNode>> = Vec::new();
    let mut i = start;
    let mut visited = 0;
    while visited < n {
        if kept_edge[i] {
            // Begin a fragment at node i.
            let mut frag = vec![KeptNode::from(&nodes[i])];
            let mut j = i;
            while kept_edge[j] {
                let next = (j + 1) % n;
                frag.push(KeptNode::from(&nodes[next]));
                j = next;
                visited += 1;
                if visited >= n {
                    break;
                }
            }
            frags.push(frag);
            i = j;
        } else {
            i = (i + 1) % n;
            visited += 1;
        }
    }
    frags
}

/// Decide whether the directed edge between endpoint signs is a kept boundary
/// edge. Dropped if either endpoint is strictly dropped; kept otherwise (both
/// kept-or-on). An On–On edge is kept (it is a coincident boundary).
fn edge_is_kept(a: Sign3, b: Sign3, keep: KeepSide) -> bool {
    let dropped = |s: Sign3| match keep {
        KeepSide::Below => s == Sign3::Above,
        KeepSide::Above => s == Sign3::Below,
    };
    !dropped(a) && !dropped(b)
}

/// A node of a kept fragment exposing the portal API.
#[derive(Debug, Clone, Copy)]
struct KeptNode {
    vertex_out: Id<Vertex>,
    point: Point3,
    sign: Sign3,
    edge_curve: CurveId,
    edge_geom: CurveGeom,
    edge_b0: f64,
    edge_b1: f64,
}

impl From<&BoundaryNode> for KeptNode {
    fn from(n: &BoundaryNode) -> Self {
        Self {
            vertex_out: n.vertex_out,
            point: n.point,
            sign: n.sign,
            edge_curve: n.edge_curve,
            edge_geom: n.edge_geom,
            edge_b0: n.edge_b0,
            edge_b1: n.edge_b1,
        }
    }
}

/// `true` for a sign strictly on the dropped side.
fn is_dropped(sign: Sign3, keep: KeepSide) -> bool {
    matches!(
        (keep, sign),
        (KeepSide::Below, Sign3::Above) | (KeepSide::Above, Sign3::Below)
    )
}

/// Nest built loops into outer + holes by exact 2-D containment.
///
/// Loops are projected by their own point ring; the loop with the largest
/// absolute area that contains a smaller loop is its outer. For the simple
/// caps produced here the common case is one outer per face (with possibly one
/// hole when the face had an inner loop straddling the cut).
fn nest_loops(loops: &[BuiltLoop]) -> Vec<LoopGroup> {
    // Project each loop to an arbitrary but consistent 2-D frame derived from
    // the first three non-collinear points of the loop. For trimmed planar
    // sub-faces the loops are co-planar (same face plane), so we can use a
    // shared frame from the first loop.
    if loops.is_empty() {
        return Vec::new();
    }
    // Build a frame from the first loop's plane.
    let frame = match frame_from_points(&loops[0].points) {
        Some(f) => f,
        None => {
            // Degenerate: treat every loop as its own outer.
            return loops
                .iter()
                .cloned()
                .map(|l| LoopGroup {
                    outer: l,
                    inners: Vec::new(),
                })
                .collect();
        }
    };
    let projected: Vec<Vec<[f64; 2]>> = loops
        .iter()
        .map(|l| l.points.iter().map(|&p| project_with(&frame, p)).collect())
        .collect();
    let areas: Vec<f64> = projected.iter().map(|p| signed_area_2d(p).abs()).collect();

    // For each loop, find the smallest loop that strictly contains it.
    let n = loops.len();
    let mut parent: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let rep = representative_point(&projected[i]);
        let mut best: Option<usize> = None;
        for j in 0..n {
            if i == j {
                continue;
            }
            if areas[j] > areas[i] && point_in_polygon(rep, &projected[j]) {
                match best {
                    Some(b) if areas[j] < areas[b] => best = Some(j),
                    None => best = Some(j),
                    _ => {}
                }
            }
        }
        parent[i] = best;
    }

    let mut groups: Vec<LoopGroup> = Vec::new();
    let mut group_index: HashMap<usize, usize> = HashMap::new();
    for i in 0..n {
        if parent[i].is_none() {
            group_index.insert(i, groups.len());
            groups.push(LoopGroup {
                outer: loops[i].clone(),
                inners: Vec::new(),
            });
        }
    }
    for i in 0..n {
        if let Some(p) = parent[i] {
            if let Some(&gi) = group_index.get(&p) {
                groups[gi].inners.push(loops[i].clone());
            } else {
                // Nested deeper than one level: treat as its own outer.
                groups.push(LoopGroup {
                    outer: loops[i].clone(),
                    inners: Vec::new(),
                });
            }
        }
    }
    groups
}

/// A simple 2-D frame: origin + two orthonormal in-plane axes.
#[derive(Clone, Copy)]
struct SimpleFrame {
    origin: Point3,
    u: Vec3,
    v: Vec3,
}

fn frame_from_points(points: &[Point3]) -> Option<SimpleFrame> {
    if points.len() < 3 {
        return None;
    }
    let origin = points[0];
    let mut e1 = Vec3::ZERO;
    for &p in &points[1..] {
        let d = p - origin;
        if d.norm() > 1e-9_f64 {
            e1 = d;
            break;
        }
    }
    let u = e1.try_unit()?.as_vec();
    // Find a second vector not parallel to u.
    let mut normal = Vec3::ZERO;
    for &p in &points[1..] {
        let d = p - origin;
        let c = u.cross(d);
        if c.norm() > 1e-9_f64 {
            normal = c;
            break;
        }
    }
    let n = normal.try_unit()?.as_vec();
    let v = n.cross(u);
    Some(SimpleFrame { origin, u, v })
}

fn project_with(f: &SimpleFrame, p: Point3) -> [f64; 2] {
    let d = p - f.origin;
    [d.dot(f.u), d.dot(f.v)]
}

/// An interior representative point of a projected polygon (vertex-average is
/// inside for convex; for safety use the centroid of the first triangle).
fn representative_point(poly: &[[f64; 2]]) -> [f64; 2] {
    if poly.len() < 3 {
        return poly.first().copied().unwrap_or([0.0, 0.0]);
    }
    [
        (poly[0][0] + poly[1][0] + poly[2][0]) / 3.0,
        (poly[0][1] + poly[1][1] + poly[2][1]) / 3.0,
    ]
}

// ── cap cycle assembly ────────────────────────────────────────────────────

/// A closed cap cycle: an ordered ring of (vertex, point) pairs.
#[derive(Debug, Clone)]
struct CapCycle {
    ring: Vec<(Id<Vertex>, Point3)>,
    /// Projected 2-D ring for nesting.
    proj: Vec<[f64; 2]>,
}

/// A nesting group of cap cycles: outer + holes.
struct CapGroup {
    outer: CapCycle,
    inners: Vec<CapCycle>,
}

/// Nest cap cycles into outer + holes by exact 2-D containment.
fn nest_cap_cycles(cycles: &[CapCycle]) -> Vec<CapGroup> {
    let n = cycles.len();
    let areas: Vec<f64> = cycles
        .iter()
        .map(|c| signed_area_2d(&c.proj).abs())
        .collect();
    let mut parent: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let rep = representative_point(&cycles[i].proj);
        let mut best: Option<usize> = None;
        for j in 0..n {
            if i == j {
                continue;
            }
            if areas[j] > areas[i] && point_in_polygon(rep, &cycles[j].proj) {
                match best {
                    Some(b) if areas[j] < areas[b] => best = Some(j),
                    None => best = Some(j),
                    _ => {}
                }
            }
        }
        parent[i] = best;
    }
    let mut groups: Vec<CapGroup> = Vec::new();
    let mut group_index: HashMap<usize, usize> = HashMap::new();
    for i in 0..n {
        if parent[i].is_none() {
            group_index.insert(i, groups.len());
            groups.push(CapGroup {
                outer: cycles[i].clone(),
                inners: Vec::new(),
            });
        }
    }
    for i in 0..n {
        if let Some(p) = parent[i] {
            if let Some(&gi) = group_index.get(&p) {
                groups[gi].inners.push(cycles[i].clone());
            } else {
                groups.push(CapGroup {
                    outer: cycles[i].clone(),
                    inners: Vec::new(),
                });
            }
        }
    }
    groups
}

/// Solve `n·(X(t) − plane.point) = 0` for a circle or ellipse, returning the
/// roots (angles) in `[0, 2π)`.
///
/// Substituting the conic parameterisation `X(t) = c + p·cos t + q·sin t` (with
/// `p`, `q` the two in-plane semi-axis vectors) into the plane equation gives
/// `A cos t + B sin t = C` where `A = n·p`, `B = n·q`, `C = −n·(c − plane.point)`.
/// With `R = hypot(A, B)`, `φ = atan2(B, A)` the equation is `R cos(t − φ) = C`,
/// whose roots are `t = φ ± acos(C / R)` when `|C| ≤ R` (a tangent when `|C| = R`,
/// none when `|C| > R`).
fn split_conic_param(c: Point3, p_vec: Vec3, q_vec: Vec3, plane: &Plane) -> Vec<f64> {
    let n = plane.normal().as_vec();
    let a = n.dot(p_vec);
    let b = n.dot(q_vec);
    let cc = -n.dot(c - plane.point());
    let r = (a * a + b * b).sqrt();
    if r == 0.0 {
        return Vec::new();
    }
    let ratio = cc / r;
    // Outside [−1, 1] beyond a hair (tangent tolerance) means no intersection.
    if !(-1.0..=1.0).contains(&ratio) && ratio.abs() - 1.0 > 1e-9_f64 {
        return Vec::new();
    }
    let clamped = ratio.clamp(-1.0, 1.0);
    let phi = b.atan2(a);
    let delta = clamped.acos();
    let two_pi = std::f64::consts::TAU;
    let norm = |t: f64| {
        let mut x = t % two_pi;
        if x < 0.0 {
            x += two_pi;
        }
        x
    };
    if delta < 1e-12_f64 {
        vec![norm(phi)]
    } else {
        vec![norm(phi + delta), norm(phi - delta)]
    }
}
