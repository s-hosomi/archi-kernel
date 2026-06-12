//! Extrude a 2-D profile along an axis into a closed solid.
//!
//! The cross-section lies on the plane through `axis.origin` perpendicular to
//! `axis.dir`; the swept solid runs from there to `axis.origin + dir·length`.
//! The 2-D → 3-D lift uses the deterministic orthonormal basis `(u, v)` derived
//! from `dir` by [`plane_basis`](crate::primitives::plane_basis) — the same
//! seed rule [`Circle3::point_at`](crate::primitives::Circle3::point_at) uses —
//! so circular cross-sections stay parameter-consistent with their edge curves.
//! With `u × v = dir`, a CCW outline swept along `+dir` produces outward-facing
//! side normals.
//!
//! Every plane (bottom cap, top cap, each side) is registered through
//! [`GeomStore::insert_plane`](crate::geom::GeomStore::insert_plane), so the
//! canonicalisation de-duplicates coincident planes and the returned `flipped`
//! flag is folded into each face's [`Sense`].
//!
//! Round sections are handled per `DESIGN.md` §6-1: the side surface is split at
//! a two-point seam into two half-cylinder faces, and each cap is a single disk
//! face whose boundary is two semicircular arcs. This avoids a single
//! closed-edge loop, which is degenerate for the later splitting phases.

use std::collections::HashMap;

use crate::brep::Brep;
use crate::csg::Profile2d;
use crate::error::KernelError;
use crate::geom::{CurveGeom, CurveId, SurfaceGeom, VertexGeom};
use crate::math::{Point3, Unit3, Vec3};
use crate::primitives::{plane_basis, Circle3, Cylinder, Line3, Plane};
use crate::profile::ProfileGeom;
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::{Face, HalfEdge, Loop, Sense, Shell, Solid, Vertex};

/// Extrude `profile` along `axis` for `length` metres into a closed [`Brep`].
///
/// The section sits on the plane through `axis.origin()` perpendicular to
/// `axis.dir()`; the solid is swept to `axis.origin() + dir·length`. The result
/// validates at [`ValidateLevel::Full`](crate::topo::ValidateLevel::Full).
///
/// # Errors
///
/// * [`KernelError::NonPositiveDimension`] if `length` is not strictly positive,
///   or if the profile outline is degenerate (see
///   [`Profile2d::outline`](crate::csg::Profile2d#method.outline)).
pub fn extrude(
    profile: &Profile2d,
    axis: &Line3,
    length: f64,
    tol: &Tol,
) -> Result<Brep, KernelError> {
    if length <= 0.0 {
        return Err(KernelError::NonPositiveDimension {
            name: "length",
            value: length,
        });
    }
    let outline = profile.outline()?;

    let dir = axis.dir();
    let (u, v) = plane_basis(dir);
    let mut builder = ExtrudeBuilder::new(*tol, *axis, u, v, length);

    match outline {
        ProfileGeom::Polygon(ring) => builder.build_polygon(&ring),
        ProfileGeom::Circle { radius } => builder.build_circle(radius),
    }

    Ok(builder.finish())
}

/// Lift a 2-D profile point to 3-D at height `s` along the axis.
struct ExtrudeBuilder {
    brep: Brep,
    tol: Tol,
    axis: Line3,
    u: Vec3,
    v: Vec3,
    length: f64,
    /// Deduplicated vertices, keyed on quantised coordinates.
    vertices: HashMap<CoordKey, Id<Vertex>>,
    coords: HashMap<CoordKey, Point3>,
    /// Deduplicated straight edges, keyed on the unordered vertex-coordinate
    /// pair, mapping to the shared line curve and its parameter-origin key.
    lines: HashMap<(CoordKey, CoordKey), (CurveId, CoordKey)>,
    faces: Vec<Id<Face>>,
}

/// Integer-quantised coordinate so identical corners deduplicate despite `f64`
/// not being `Hash`/`Eq`. The scale (1e9) resolves building dimensions
/// (1e-3..1e2 m) comfortably without overflowing `i64`.
type CoordKey = (i64, i64, i64);

fn key(p: Point3) -> CoordKey {
    let q = |x: f64| (x * 1.0e9_f64).round() as i64;
    (q(p.x), q(p.y), q(p.z))
}

impl ExtrudeBuilder {
    fn new(tol: Tol, axis: Line3, u: Vec3, v: Vec3, length: f64) -> Self {
        Self {
            brep: Brep::new(),
            tol,
            axis,
            u,
            v,
            length,
            vertices: HashMap::new(),
            coords: HashMap::new(),
            lines: HashMap::new(),
            faces: Vec::new(),
        }
    }

    /// The extrusion direction.
    #[inline]
    fn dir(&self) -> Unit3 {
        self.axis.dir()
    }

    /// Lift a local 2-D profile point to 3-D at axial height `s`.
    fn lift(&self, p: [f64; 2], s: f64) -> Point3 {
        self.axis.origin() + self.u * p[0] + self.v * p[1] + self.dir().as_vec() * s
    }

    /// Get or create the shared vertex at a 3-D point.
    fn vertex(&mut self, p: Point3) -> Id<Vertex> {
        let k = key(p);
        if let Some(&v) = self.vertices.get(&k) {
            return v;
        }
        let pid = self.brep.geom.insert_point(VertexGeom::Explicit(p));
        let v = self.brep.topo.add_vertex(Vertex { point: pid });
        self.vertices.insert(k, v);
        self.coords.insert(k, p);
        v
    }

    /// Get or create the shared straight-line curve between two points,
    /// returning the curve and the coord key chosen as its parameter origin.
    fn line_curve(&mut self, a: Point3, b: Point3) -> (CurveId, CoordKey) {
        let (ka, kb) = (key(a), key(b));
        let unordered = if ka <= kb { (ka, kb) } else { (kb, ka) };
        if let Some(&entry) = self.lines.get(&unordered) {
            return entry;
        }
        let (origin_key, origin_pt, other_pt) = if ka <= kb { (ka, a, b) } else { (kb, b, a) };
        let line = Line3::new(origin_pt, other_pt - origin_pt).expect("non-degenerate edge");
        let cid = self.brep.geom.insert_curve(CurveGeom::Line(line));
        let entry = (cid, origin_key);
        self.lines.insert(unordered, entry);
        entry
    }

    /// Add a straight half-edge from `a` to `b` on its shared line, with the
    /// boundary parameterised by signed arc length from the line origin.
    fn line_half_edge(&mut self, a: Point3, b: Point3) -> Id<HalfEdge> {
        let start = self.vertex(a);
        let _ = self.vertex(b);
        let (curve, origin_key) = self.line_curve(a, b);
        let origin_pt = self.coords[&origin_key];
        let line = match self.brep.geom.curve(curve).expect("line curve") {
            CurveGeom::Line(l) => *l,
            _ => unreachable!("line_curve stores a line"),
        };
        let ta = (a - origin_pt).dot(line.dir().as_vec());
        let tb = (b - origin_pt).dot(line.dir().as_vec());
        self.brep.topo.add_half_edge(HalfEdge {
            start,
            curve,
            boundary: [ta, tb],
        })
    }

    /// Add an arc half-edge on a shared circle curve, with an angular boundary.
    fn arc_half_edge(
        &mut self,
        circle: CurveId,
        start: Id<Vertex>,
        a: f64,
        b: f64,
    ) -> Id<HalfEdge> {
        self.brep.topo.add_half_edge(HalfEdge {
            start,
            curve: circle,
            boundary: [a, b],
        })
    }

    /// Build a loop from a list of half-edges.
    fn loop_of(&mut self, hes: Vec<Id<HalfEdge>>) -> Id<Loop> {
        self.brep.topo.add_loop(Loop { half_edges: hes })
    }

    /// Register a plane (canonicalised) and turn the `flipped` flag into a
    /// face sense.
    fn plane_face(&mut self, plane: Plane, outer: Id<Loop>) -> Id<Face> {
        let (surface, flipped) = self.brep.geom.insert_plane(plane, &self.tol);
        let sense = if flipped {
            Sense::Reversed
        } else {
            Sense::Same
        };
        self.brep.topo.add_face(Face {
            surface,
            sense,
            outer,
            inners: Vec::new(),
        })
    }

    /// Register a cylinder surface verbatim and build a side face.
    fn cylinder_face(&mut self, cyl: Cylinder, outer: Id<Loop>) -> Id<Face> {
        let surface = self.brep.geom.insert_surface(SurfaceGeom::Cylinder(cyl));
        self.brep.topo.add_face(Face {
            surface,
            sense: Sense::Same,
            outer,
            inners: Vec::new(),
        })
    }

    fn push_face(&mut self, f: Id<Face>) {
        self.faces.push(f);
    }

    // ── polygon extrusion ──────────────────────────────────────────────────

    /// Build the solid for a CCW polygonal profile: bottom + top caps and one
    /// side face per edge.
    fn build_polygon(&mut self, ring: &[[f64; 2]]) {
        let n = ring.len();
        let bottom: Vec<Point3> = ring.iter().map(|&p| self.lift(p, 0.0)).collect();
        let top: Vec<Point3> = ring.iter().map(|&p| self.lift(p, self.length)).collect();

        // Bottom cap: outward normal −dir, so traverse the ring reversed
        // (CCW as seen from below).
        let bottom_ring: Vec<Point3> = bottom.iter().rev().copied().collect();
        let hes: Vec<Id<HalfEdge>> = (0..n)
            .map(|i| {
                let a = bottom_ring[i];
                let b = bottom_ring[(i + 1) % n];
                self.line_half_edge(a, b)
            })
            .collect();
        let lp = self.loop_of(hes);
        let f = self.plane_face(Plane::new_unchecked(bottom[0], flip_unit(self.dir())), lp);
        self.push_face(f);

        // Top cap: outward normal +dir, ring order.
        let hes: Vec<Id<HalfEdge>> = (0..n)
            .map(|i| {
                let a = top[i];
                let b = top[(i + 1) % n];
                self.line_half_edge(a, b)
            })
            .collect();
        let lp = self.loop_of(hes);
        let f = self.plane_face(Plane::new_unchecked(top[0], self.dir()), lp);
        self.push_face(f);

        // Side faces, one per edge: B_i → B_{i+1} → T_{i+1} → T_i.
        for i in 0..n {
            let j = (i + 1) % n;
            let (bi, bj, ti, tj) = (bottom[i], bottom[j], top[j], top[i]);
            let h0 = self.line_half_edge(bi, bj);
            let h1 = self.line_half_edge(bj, ti);
            let h2 = self.line_half_edge(ti, tj);
            let h3 = self.line_half_edge(tj, bi);
            let lp = self.loop_of(vec![h0, h1, h2, h3]);
            // Outward normal: edge direction (in-plane) crossed with the axis.
            let edge = bj - bi;
            let n_out = edge.cross(self.dir().as_vec());
            let plane = Plane::new(bi, n_out).expect("non-degenerate side normal");
            let f = self.plane_face(plane, lp);
            self.push_face(f);
        }
    }

    // ── circular extrusion ─────────────────────────────────────────────────

    /// Build the solid for a circular profile: two half-cylinder side faces and
    /// two disk caps split at a two-point seam (`DESIGN.md` §6-1).
    fn build_circle(&mut self, radius: f64) {
        use std::f64::consts::PI;
        let two_pi = 2.0 * PI;
        let dir = self.dir();

        // Bottom and top circles, sharing the axis basis (so the angle
        // parameters match `Circle3::point_at`).
        let bottom_centre = self.lift([0.0, 0.0], 0.0);
        let top_centre = self.lift([0.0, 0.0], self.length);
        let bottom_circle = Circle3::new_unchecked(bottom_centre, dir, radius);
        let top_circle = Circle3::new_unchecked(top_centre, dir, radius);

        // Seam points A (φ = 0) and B (φ = π) on each rim.
        let a_bot = bottom_circle.point_at(0.0);
        let b_bot = bottom_circle.point_at(PI);
        let a_top = top_circle.point_at(0.0);
        let b_top = top_circle.point_at(PI);

        // Pre-create the shared vertices so the arcs reference identical ids.
        let va_bot = self.vertex(a_bot);
        let vb_bot = self.vertex(b_bot);
        let va_top = self.vertex(a_top);
        let vb_top = self.vertex(b_top);

        let cb = self
            .brep
            .geom
            .insert_curve(CurveGeom::Circle(bottom_circle));
        let ct = self.brep.geom.insert_curve(CurveGeom::Circle(top_circle));

        // Cylinder surface for the sides.
        let cyl = Cylinder::new_unchecked(self.axis, radius);

        // Bottom cap (outward −dir): two arcs 2π→π→0 (CCW as seen from below).
        let h0 = self.arc_half_edge(cb, va_bot, two_pi, PI);
        let h1 = self.arc_half_edge(cb, vb_bot, PI, 0.0);
        let lp = self.loop_of(vec![h0, h1]);
        let f = self.plane_face(Plane::new_unchecked(bottom_centre, flip_unit(dir)), lp);
        self.push_face(f);

        // Top cap (outward +dir): two arcs 0→π→2π.
        let h0 = self.arc_half_edge(ct, va_top, 0.0, PI);
        let h1 = self.arc_half_edge(ct, vb_top, PI, two_pi);
        let lp = self.loop_of(vec![h0, h1]);
        let f = self.plane_face(Plane::new_unchecked(top_centre, dir), lp);
        self.push_face(f);

        // Side 1: φ ∈ [0, π].  bottom arc → seam(B) up → top arc → seam(A) down.
        let s_bottom = self.arc_half_edge(cb, va_bot, 0.0, PI);
        let s_up = self.line_half_edge(b_bot, b_top);
        let s_top = self.arc_half_edge(ct, vb_top, PI, 0.0);
        let s_down = self.line_half_edge(a_top, a_bot);
        let lp = self.loop_of(vec![s_bottom, s_up, s_top, s_down]);
        let f = self.cylinder_face(cyl, lp);
        self.push_face(f);

        // Side 2: φ ∈ [π, 2π]. bottom arc → seam(A) up → top arc → seam(B) down.
        let s_bottom = self.arc_half_edge(cb, vb_bot, PI, two_pi);
        let s_up = self.line_half_edge(a_bot, a_top);
        let s_top = self.arc_half_edge(ct, va_top, two_pi, PI);
        let s_down = self.line_half_edge(b_top, b_bot);
        let lp = self.loop_of(vec![s_bottom, s_up, s_top, s_down]);
        let f = self.cylinder_face(cyl, lp);
        self.push_face(f);
    }

    /// Wrap the accumulated faces into a single-shell solid.
    fn finish(mut self) -> Brep {
        let shell = self.brep.topo.add_shell(Shell {
            faces: self.faces.clone(),
        });
        let solid = self.brep.topo.add_solid(Solid {
            shells: vec![shell],
        });
        self.brep.solids = vec![solid];
        self.brep
    }
}

/// Flip a unit vector (kept as a `Unit3` since negation preserves the norm).
fn flip_unit(u: Unit3) -> Unit3 {
    Unit3::new_unchecked(-u.as_vec())
}
