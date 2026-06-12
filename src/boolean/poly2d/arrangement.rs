//! Planar arrangement: snap → split → dedup → DCEL → faces.
//!
//! This module turns the union of both operands' edges into a planar
//! subdivision (an arrangement) and extracts its bounded faces. It is the
//! structural core of the engine; the public boolean ops are a thin layer of
//! winding classification on top.
//!
//! # Pipeline
//!
//! 1. **Snap / merge** every input vertex through [`crate::boolean::poly2d::snap::VertexStore`]
//!    so coincident points become one vertex (degeneracy collapse).
//! 2. **Intersect** every pair of edges and snap the crossing points, so each
//!    edge knows all the vertices that split it.
//! 3. **Split** each input edge at its interior vertices into atomic edges
//!    between consecutive vertices. Straight inputs split into sub-segments;
//!    arc inputs split into sub-arcs ordered by angle along the sweep.
//! 4. **Dedup** coincident atomic edges: two atomic edges between the same
//!    vertex pair *with the same geometry* (both straight, or both arcs of the
//!    same circle) are the same arrangement edge — this is where shared /
//!    overlapping boundary edges in buildings collapse. Each surviving
//!    arrangement edge becomes a pair of opposite **half-edges**.
//! 5. **Build the DCEL**: sort half-edges leaving each vertex by **tangent
//!    angle** (so arcs order correctly against straight edges), link
//!    `next`/`prev` by the "next clockwise" rule, and trace face loops.
//! 6. **Extract faces**: every loop is kept; the unbounded outer wrap always
//!    classifies as outside both operands and is therefore never selected.
//!
//! Robustness rests on the exact [`orient2d`] (combinatorial decisions for
//! segments) plus the snap (coincidence collapse). Float coordinates are used
//! only for *where*, never for *whether*; arc sidedness uses an exact radial
//! test (inside / outside the circle) which is unambiguous away from the curve.

use std::collections::{BTreeMap, HashMap};

use crate::boolean::poly2d::error::Poly2Error;
use crate::boolean::poly2d::geom::{
    directed_sweep, eps_sq, orient2d, Arc, Edge2, Orient, Point2, Vec2,
};
use crate::boolean::poly2d::intersect::intersect;
use crate::boolean::poly2d::region::{Contour, Region};
use crate::boolean::poly2d::snap::{VertexId, VertexStore};
use crate::boolean::support::quantize;
use crate::tolerance::Tol;

/// Which operand an input edge came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operand {
    /// The first operand (`A`).
    A,
    /// The second operand (`B`).
    B,
}

/// The geometry kind of an input or arrangement edge.
#[derive(Debug, Clone, Copy)]
enum EdgeGeom {
    /// A straight segment.
    Seg,
    /// An arc on a circle: centre, radius. The traversal direction and the
    /// angular span are recovered from the edge's two endpoints plus the
    /// `ccw` flag (whether the input arc swept counter-clockwise).
    Arc {
        center: Point2,
        radius: f64,
        ccw: bool,
    },
}

/// An input edge after snapping its endpoints to vertex ids, carrying its
/// operand, original traversal direction, and geometry.
#[derive(Debug, Clone, Copy)]
struct InputEdge {
    a: VertexId,
    b: VertexId,
    operand: Operand,
    geom: EdgeGeom,
}

/// An arrangement edge (an undirected edge between two distinct vertices),
/// annotated with its geometry and the net directed winding contribution of each
/// operand.
///
/// `wind_a` / `wind_b` are the sum over all coincident input edges of `+1` if
/// the input ran `a → b` and `-1` if it ran `b → a`. After dedup an arrangement
/// edge with `wind == 0` for *both* operands carries no boundary and is dropped.
#[derive(Debug, Clone, Copy)]
struct ArrEdge {
    a: VertexId,
    b: VertexId,
    wind_a: i32,
    wind_b: i32,
    /// Geometry of the edge in its **canonical** `(a ≤ b)` direction.
    geom: EdgeGeom,
}

/// A directed half-edge in the DCEL.
#[derive(Debug, Clone, Copy)]
struct HalfEdge {
    origin: VertexId,
    dest: VertexId,
    twin: usize,
    next: usize,
    /// Geometry as traversed `origin → dest`.
    geom: EdgeGeom,
    /// Index of the face loop this half-edge belongs to (assigned during trace).
    face: Option<usize>,
}

/// One directed boundary edge of a traced face loop, with the geometry needed to
/// reconstruct an arc.
#[derive(Debug, Clone, Copy)]
pub struct LoopEdge {
    /// Start vertex id.
    pub a: VertexId,
    /// End vertex id.
    pub b: VertexId,
    /// Geometry as traversed `a → b`.
    geom: EdgeGeom,
}

/// A traced face loop: the ordered boundary edges of the face on its left.
#[derive(Debug, Clone)]
pub struct FaceLoop {
    /// Directed boundary edges in traversal order.
    pub edges: Vec<LoopEdge>,
    /// Start-point coordinates of each edge, parallel to `edges`.
    pub vertices: Vec<Point2>,
}

impl FaceLoop {
    /// Signed area of the loop (positive = CCW bounded face / outer boundary,
    /// negative = CW hole boundary or the unbounded outer wrap), arc-aware.
    #[cfg(test)]
    pub fn signed_area(&self) -> f64 {
        let mut acc = 0.0_f64;
        let n = self.vertices.len();
        for i in 0..n {
            let a = self.vertices[i];
            let b = self.vertices[(i + 1) % n];
            acc += a.x * b.y - b.x * a.y;
        }
        let mut area = 0.5 * acc;
        for (i, e) in self.edges.iter().enumerate() {
            if let EdgeGeom::Arc {
                center,
                ccw,
                radius,
            } = e.geom
            {
                let pa = self.vertices[i];
                let pb = self.vertices[(i + 1) % n];
                let dtheta = directed_sweep(center, pa, pb, ccw);
                area += 0.5 * radius * radius * (dtheta - dtheta.sin());
            }
        }
        area
    }
}

/// The built arrangement.
pub struct Arrangement {
    store: VertexStore,
    inputs: Vec<InputEdge>,
    /// Bounded faces, in no particular order.
    pub faces: Vec<FaceLoop>,
}

impl Arrangement {
    /// Build the arrangement from two operand regions.
    pub fn build(a: &Region, b: &Region, tol: &Tol) -> Result<Self, Poly2Error> {
        validate_region(a, tol)?;
        validate_region(b, tol)?;

        let mut store = VertexStore::new(*tol);
        let mut inputs: Vec<InputEdge> = Vec::new();

        // Vertex-on-edge pre-snap is segment-only (the grazing case the building
        // domain needs); arc inputs skip it (their endpoints are seam points that
        // already snap vertex-to-vertex).
        let a = project_grazing_vertices(a, b, tol);
        let b = project_grazing_vertices(b, &a, tol);

        ingest(&a, Operand::A, &mut store, &mut inputs, tol)?;
        ingest(&b, Operand::B, &mut store, &mut inputs, tol)?;

        // Intersect every pair, snap crossings, collect split points.
        let n = inputs.len();
        let mut split_points: Vec<Vec<VertexId>> = vec![Vec::new(); n];
        for i in 0..n {
            let ei = input_edge2(&store, &inputs[i]);
            for j in (i + 1)..n {
                let ej = input_edge2(&store, &inputs[j]);
                let cr = intersect(&ei, &ej, tol)?;
                for p in cr.points {
                    let v = store.insert(p);
                    push_split(&mut split_points[i], v, &inputs[i]);
                    push_split(&mut split_points[j], v, &inputs[j]);
                }
            }
        }

        // Vertex-on-edge snapping for straight edges (segment grazing only).
        let vcount = store.len();
        for i in 0..inputs.len() {
            if !matches!(inputs[i].geom, EdgeGeom::Seg) {
                continue;
            }
            let (sa, sb) = (inputs[i].a, inputs[i].b);
            let a = store.point(sa);
            let b = store.point(sb);
            for vid in 0..vcount {
                let v = VertexId(vid);
                if v == sa || v == sb {
                    continue;
                }
                let p = store.point(v);
                if point_on_segment_interior(a, b, p, tol) {
                    push_split(&mut split_points[i], v, &inputs[i]);
                }
            }
        }

        // Split each input into atomic edges and dedup by (vertex pair, geom).
        let mut arr_map: BTreeMap<(VertexId, VertexId, GeomKey), ArrEdge> = BTreeMap::new();
        for (i, edge) in inputs.iter().enumerate() {
            let chain = ordered_chain(&store, edge, &split_points[i], tol);
            for w in chain.windows(2) {
                let (u, v) = (w[0], w[1]);
                if u == v {
                    continue;
                }
                accumulate_edge(&mut store, &mut arr_map, u, v, edge);
            }
        }

        let arr_edges: Vec<ArrEdge> = arr_map
            .into_values()
            .filter(|e| e.wind_a != 0 || e.wind_b != 0)
            .collect();

        let halfs = build_dcel(&store, &arr_edges, tol)?;
        let faces = trace_faces(&store, halfs);

        Ok(Self {
            store,
            inputs,
            faces,
        })
    }

    /// Winding number of `p` with respect to the given operand's *original*
    /// input edges (before arrangement), counting both straight and arc edges.
    pub fn winding(&self, p: Point2, operand: Operand) -> i32 {
        let mut wind = 0_i32;
        for edge in self.inputs.iter().filter(|s| s.operand == operand) {
            let a = self.store.point(edge.a);
            let b = self.store.point(edge.b);
            match edge.geom {
                EdgeGeom::Seg => wind += ray_cross_seg(p, a, b),
                EdgeGeom::Arc {
                    center,
                    radius,
                    ccw,
                } => {
                    wind += ray_cross_arc(p, a, b, center, radius, ccw);
                }
            }
        }
        wind
    }

    /// Borrow the snapped vertex store (for reconstruction).
    #[inline]
    pub fn store(&self) -> &VertexStore {
        &self.store
    }

    /// The `(in_a, in_b)` winding classification of the **face to the left** of a
    /// directed loop, robust against thin neighbouring features and curved-edge
    /// tangencies.
    ///
    /// The sample is offset from an edge midpoint along its left normal, then
    /// shrunk until the **loop's own self-winding** at the sample matches the
    /// target for "the face on the left": `+1` for a CCW loop (its interior),
    /// `0` for a CW loop (the surrounding region). Anchoring on this topological
    /// invariant — rather than on the raw step landing in the right place — makes
    /// the classification robust where a curved edge's midpoint sits at a
    /// tangent contact with another circle (the step would otherwise slip into
    /// the neighbouring disc).
    pub fn loop_sample_point(&self, face: &FaceLoop) -> (i32, i32) {
        let classify = |p: Point2| (self.winding(p, Operand::A), self.winding(p, Operand::B));
        let n = face.edges.len();
        if n == 0 {
            return (0, 0);
        }
        // Pick the longest edge (by chord) for a stable normal.
        let mut best_i = 0usize;
        let mut best_len = -1.0_f64;
        for i in 0..n {
            let e = &face.edges[i];
            let l = self.store.point(e.a).dist(self.store.point(e.b));
            if l > best_len {
                best_len = l;
                best_i = i;
            }
        }
        let e = &face.edges[best_i];
        let edge2 = arr_edge2(&self.store, e.a, e.b, e.geom);
        // Sample at a **generic fraction** of the edge (not the midpoint): an arc
        // midpoint can sit exactly at the circle's diametral height, where a
        // horizontal winding ray grazes the circle and double-counts. A fraction
        // like 0.37 keeps the base point off every circle's centre height and off
        // the seam, so the winding ray-cast stays unambiguous.
        let base = edge2.point_at(0.37);
        let tan = edge_tangent_at_fraction(&edge2, 0.37);
        if tan.len_sq() <= 0.0 {
            return classify(self.centroid(face));
        }
        // Left normal of the tangent.
        let nx = -tan.y;
        let ny = tan.x;
        let sample = |s: f64| Point2::new(base.x + nx * s, base.y + ny * s);

        // Target self-winding for "the face on the left".
        let target = if loop_signed_area(&self.store, face) > 0.0 {
            1
        } else {
            0
        };

        // Shrink the step until the sample lands in the loop's own left face
        // (self-winding == target) *and* the global classification is stable.
        let mut step = 1.0e-6_f64;
        let mut last = classify(sample(step));
        for _ in 0..40 {
            let p = sample(step);
            if self.loop_self_winding(face, p) == target {
                return classify(p);
            }
            last = classify(p);
            step *= 0.5;
        }
        last
    }

    /// Winding number of `p` with respect to **one face loop's own edges**, used
    /// to confirm a sample lies in the loop's left face.
    fn loop_self_winding(&self, face: &FaceLoop, p: Point2) -> i32 {
        let mut w = 0_i32;
        for e in &face.edges {
            let a = self.store.point(e.a);
            let b = self.store.point(e.b);
            match e.geom {
                EdgeGeom::Seg => w += ray_cross_seg(p, a, b),
                EdgeGeom::Arc {
                    center,
                    radius,
                    ccw,
                } => {
                    w += ray_cross_arc(p, a, b, center, radius, ccw);
                }
            }
        }
        w
    }

    fn centroid(&self, face: &FaceLoop) -> Point2 {
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for v in &face.vertices {
            cx += v.x;
            cy += v.y;
        }
        let k = face.vertices.len().max(1) as f64;
        Point2::new(cx / k, cy / k)
    }
}

impl LoopEdge {
    /// Build the [`Edge2`] this loop edge represents, given the snapped store.
    pub fn to_edge2(self, store: &VertexStore) -> Edge2 {
        arr_edge2(store, self.a, self.b, self.geom)
    }
}

/// Unit tangent of an edge at fraction `t ∈ [0, 1]`, in traversal direction.
fn edge_tangent_at_fraction(edge: &Edge2, t: f64) -> Vec2 {
    match edge {
        Edge2::Seg { start, end } => {
            let d = start.to(*end);
            let l = d.len();
            if l > 0.0 {
                Vec2::new(d.x / l, d.y / l)
            } else {
                Vec2::new(0.0, 0.0)
            }
        }
        Edge2::Arc(a) => {
            let theta = a.start_angle + a.sweep * t;
            let s = if a.sweep >= 0.0 { 1.0 } else { -1.0 };
            Vec2::new(-theta.sin() * s, theta.cos() * s)
        }
    }
}

/// Arc-aware signed area of a traced face loop (positive = CCW).
fn loop_signed_area(store: &VertexStore, face: &FaceLoop) -> f64 {
    let n = face.edges.len();
    if n == 0 {
        return 0.0;
    }
    let pts: Vec<Point2> = face.edges.iter().map(|e| store.point(e.a)).collect();
    let mut acc = 0.0_f64;
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        acc += a.x * b.y - b.x * a.y;
    }
    let mut area = 0.5 * acc;
    for (i, e) in face.edges.iter().enumerate() {
        if let EdgeGeom::Arc {
            center,
            radius,
            ccw,
        } = e.geom
        {
            let pa = pts[i];
            let pb = pts[(i + 1) % n];
            let dtheta = directed_sweep(center, pa, pb, ccw);
            area += 0.5 * radius * radius * (dtheta - dtheta.sin());
        }
    }
    area
}

/// A coarse geometry discriminator for the dedup key: straight edges share one
/// key; arcs are keyed on quantised circle centre + radius so two arc fragments
/// of the *same* circle dedup while distinct circles stay separate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum GeomKey {
    Seg,
    /// An arc keyed on its circle (centre + radius) **and** a quantised midpoint
    /// of the arc as traversed in canonical `(a ≤ b)` order. The midpoint is
    /// essential: two semicircles of the *same* circle between the *same* seam
    /// vertices (upper vs lower) share the circle key but have opposite midpoints,
    /// so they must not dedup into one (which would cancel both and erase the
    /// circle). Distinct fragments of one input arc still share circle + midpoint
    /// only when they are the same atomic edge.
    Arc {
        cx: i64,
        cy: i64,
        r: i64,
        mx: i64,
        my: i64,
    },
}

/// Compute the dedup key for an atomic edge whose **canonical** geometry runs
/// `pa → pb`.
fn geom_key_for(pa: Point2, pb: Point2, geom: &EdgeGeom) -> GeomKey {
    match geom {
        EdgeGeom::Seg => GeomKey::Seg,
        EdgeGeom::Arc {
            center,
            radius,
            ccw,
        } => {
            let sa = (pa.y - center.y).atan2(pa.x - center.x);
            let sweep = directed_sweep(*center, pa, pb, *ccw);
            let mid = Point2::new(
                center.x + radius * (sa + 0.5 * sweep).cos(),
                center.y + radius * (sa + 0.5 * sweep).sin(),
            );
            GeomKey::Arc {
                cx: quantize(center.x),
                cy: quantize(center.y),
                r: quantize(*radius),
                mx: quantize(mid.x),
                my: quantize(mid.y),
            }
        }
    }
}

/// Build the [`Edge2`] for an atomic edge from `a → b` with the given geometry.
fn arr_edge2(store: &VertexStore, a: VertexId, b: VertexId, geom: EdgeGeom) -> Edge2 {
    let pa = store.point(a);
    let pb = store.point(b);
    match geom {
        EdgeGeom::Seg => Edge2::seg(pa, pb),
        EdgeGeom::Arc {
            center,
            radius,
            ccw,
        } => {
            let start_angle = (pa.y - center.y).atan2(pa.x - center.x);
            let sweep = directed_sweep(center, pa, pb, ccw);
            Edge2::Arc(Arc::new(center, radius, start_angle, sweep))
        }
    }
}

/// Winding contribution of a directed segment `a → b`.
fn ray_cross_seg(p: Point2, a: Point2, b: Point2) -> i32 {
    if a.y <= p.y {
        if b.y > p.y && orient2d(a, b, p) == Orient::Left {
            return 1;
        }
    } else if b.y <= p.y && orient2d(a, b, p) == Orient::Right {
        return -1;
    }
    0
}

/// Winding contribution of a directed **arc** `a → b` (on its circle) to the
/// winding number of `p`, by counting signed crossings of the horizontal ray to
/// the right of `p` with the actual arc geometry.
///
/// The horizontal line `y = p.y` meets the circle in at most two points; for
/// each that lies strictly to the right of `p` *and* within the arc's angular
/// sweep, the arc contributes `+1` if it is locally rising (`dy/dθ` positive in
/// the traversal direction) or `−1` if falling — the standard non-zero winding
/// rule applied to a curved edge. This is exact away from the curve (the only
/// place it matters, since the sample point is always taken in a face interior).
fn ray_cross_arc(p: Point2, a: Point2, b: Point2, center: Point2, radius: f64, ccw: bool) -> i32 {
    let dy = p.y - center.y;
    // A ray whose height grazes the circle (|dy| ≈ radius) is **tangent**: by the
    // standard ray-casting convention a tangent touch is not a crossing (it does
    // not change inside/outside). Treating it as two coincident crossings would
    // double-count (±2). A small absolute slack catches the case where a sample's
    // y lands exactly on the circle's top/bottom — which happens when it sits at
    // another circle's centre height. The sample is always in a face interior, so
    // this never suppresses a real crossing.
    let graze = 1e-9_f64;
    if dy.abs() >= radius - graze {
        return 0;
    }
    let dx = (radius * radius - dy * dy).max(0.0).sqrt();
    // The two circle points at height p.y: x = center.x ± dx.
    let arc = {
        let start_angle = (a.y - center.y).atan2(a.x - center.x);
        let sweep = directed_sweep(center, a, b, ccw);
        Arc::new(center, radius, start_angle, sweep)
    };
    let mut acc = 0_i32;
    let tol = Tol::default();
    for &xx in &[center.x + dx, center.x - dx] {
        let pt = Point2::new(xx, p.y);
        if xx <= p.x {
            continue; // not to the right of p
        }
        // Must lie within the arc's sweep.
        if arc.angle_of_point(pt, &tol).is_none() {
            continue;
        }
        // Half-open rule: exclude the crossing when it coincides with the arc's
        // **end** endpoint, so a seam vertex shared by two arcs is counted by
        // exactly one of them (the one for which it is the start). Without this a
        // ray grazing a shared diametral vertex double-counts.
        if pt.coincident(arc.end(), &tol) {
            continue;
        }
        // Local vertical direction of the arc at this crossing. The CCW tangent at
        // angle θ is (−sinθ, cosθ); for CW negate. We need d(y)/d(travel) sign.
        let theta = (pt.y - center.y).atan2(pt.x - center.x);
        let tan_y = if ccw { theta.cos() } else { -theta.cos() };
        if tan_y > 0.0 {
            acc += 1; // rising crossing
        } else if tan_y < 0.0 {
            acc -= 1; // falling crossing
        }
    }
    acc
}

/// Reject malformed input: degenerate contours and self-intersections.
fn validate_region(r: &Region, tol: &Tol) -> Result<(), Poly2Error> {
    for (ci, c) in r.contours.iter().enumerate() {
        // An all-segment contour needs ≥3 distinct vertices; a contour with arcs
        // can enclose area with only the seam vertices, so only reject it when it
        // bounds no area.
        if c.has_arc() {
            if c.signed_area().abs() <= eps_sq(tol) {
                return Err(Poly2Error::DegenerateContour { contour_index: ci });
            }
        } else if c.distinct_vertex_count(tol) < 3 {
            return Err(Poly2Error::DegenerateContour { contour_index: ci });
        }
    }
    Ok(())
}

/// Snap every contour edge's endpoints and emit [`InputEdge`]s.
fn ingest(
    r: &Region,
    operand: Operand,
    store: &mut VertexStore,
    inputs: &mut Vec<InputEdge>,
    _tol: &Tol,
) -> Result<(), Poly2Error> {
    for c in &r.contours {
        for e in &c.edges {
            match e {
                Edge2::Seg { start, end } => {
                    let a = store.insert(*start);
                    let b = store.insert(*end);
                    if a != b {
                        inputs.push(InputEdge {
                            a,
                            b,
                            operand,
                            geom: EdgeGeom::Seg,
                        });
                    }
                }
                Edge2::Arc(arc) => {
                    let a = store.insert(arc.start());
                    let b = store.insert(arc.end());
                    inputs.push(InputEdge {
                        a,
                        b,
                        operand,
                        geom: EdgeGeom::Arc {
                            center: arc.center,
                            radius: arc.radius,
                            ccw: arc.sweep >= 0.0,
                        },
                    });
                }
            }
        }
    }
    Ok(())
}

/// Project every (segment-input) region vertex within `tol` of another operand's
/// straight edge interior onto that edge, collapsing a sub-tolerance grazing gap.
/// Arc edges are passed through unchanged.
fn project_grazing_vertices(region: &Region, other: &Region, tol: &Tol) -> Region {
    let edges: Vec<(Point2, Point2)> = other
        .contours
        .iter()
        .flat_map(|c| c.edges.iter())
        .filter_map(|e| match e {
            Edge2::Seg { start, end } => Some((*start, *end)),
            Edge2::Arc(_) => None,
        })
        .collect();

    let project = |p: Point2| -> Point2 {
        for &(a, b) in &edges {
            if let Some(proj) = project_if_on_segment_interior(a, b, p, tol) {
                return proj;
            }
        }
        p
    };

    let contours = region
        .contours
        .iter()
        .map(|c| {
            if c.has_arc() {
                // Don't reshape arc contours; pass them through verbatim.
                c.clone()
            } else {
                let pts: Vec<Point2> = c.vertices().into_iter().map(project).collect();
                Contour::from_points(&pts)
            }
        })
        .collect();
    Region::new(contours)
}

fn project_if_on_segment_interior(a: Point2, b: Point2, p: Point2, tol: &Tol) -> Option<Point2> {
    let d = a.to(b);
    let len_sq = d.len_sq();
    if len_sq <= eps_sq(tol) {
        return None;
    }
    let t_raw = a.to(p).dot(d) / len_sq;
    let t = t_raw.clamp(0.0, 1.0);
    let foot = Point2::new(a.x + d.x * t, a.y + d.y * t);
    if foot.dist_sq(p) > eps_sq(tol) {
        return None;
    }
    Some(foot)
}

fn point_on_segment_interior(a: Point2, b: Point2, p: Point2, tol: &Tol) -> bool {
    let d = a.to(b);
    let len_sq = d.len_sq();
    if len_sq <= eps_sq(tol) {
        return false;
    }
    let t = a.to(p).dot(d) / len_sq;
    let margin = tol.length / len_sq.sqrt();
    if t <= margin || t >= 1.0 - margin {
        return false;
    }
    let cross = a.to(p).cross(d);
    (cross * cross) <= eps_sq(tol) * len_sq
}

/// Reconstruct the [`Edge2`] of an input edge from snapped coordinates.
#[inline]
fn input_edge2(store: &VertexStore, e: &InputEdge) -> Edge2 {
    arr_edge2(store, e.a, e.b, e.geom)
}

/// Record `v` as a split point of input edge `s` if it is interior.
fn push_split(splits: &mut Vec<VertexId>, v: VertexId, s: &InputEdge) {
    if v != s.a && v != s.b && !splits.contains(&v) {
        splits.push(v);
    }
}

/// Order an input edge's vertices `a, (splits…), b` along the edge: by chord
/// projection for a segment, by angle along the sweep for an arc.
fn ordered_chain(
    store: &VertexStore,
    s: &InputEdge,
    splits: &[VertexId],
    tol: &Tol,
) -> Vec<VertexId> {
    let a = store.point(s.a);
    let b = store.point(s.b);
    match s.geom {
        EdgeGeom::Seg => {
            let d = a.to(b);
            let len_sq = d.len_sq().max(f64::MIN_POSITIVE);
            let mut mids: Vec<(f64, VertexId)> = splits
                .iter()
                .map(|&v| (a.to(store.point(v)).dot(d) / len_sq, v))
                .collect();
            mids.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
            let mut chain = Vec::with_capacity(mids.len() + 2);
            chain.push(s.a);
            chain.extend(mids.into_iter().map(|(_, v)| v));
            chain.push(s.b);
            chain
        }
        EdgeGeom::Arc {
            center,
            radius,
            ccw,
        } => {
            let arc = {
                let start_angle = (a.y - center.y).atan2(a.x - center.x);
                let sweep = directed_sweep(center, a, b, ccw);
                Arc::new(center, radius, start_angle, sweep)
            };
            // Parameter along the sweep, in [0, 1], for each split.
            let span = arc.sweep.abs().max(f64::MIN_POSITIVE);
            let mut mids: Vec<(f64, VertexId)> = splits
                .iter()
                .filter_map(|&v| {
                    let p = store.point(v);
                    arc.angle_of_point(p, tol).map(|theta| {
                        let off = (theta - arc.start_angle).abs() / span;
                        (off, v)
                    })
                })
                .collect();
            mids.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
            let mut chain = Vec::with_capacity(mids.len() + 2);
            chain.push(s.a);
            chain.extend(mids.into_iter().map(|(_, v)| v));
            chain.push(s.b);
            chain
        }
    }
}

/// Add a directed atomic edge `u → v` to the arrangement-edge map.
fn accumulate_edge(
    store: &mut VertexStore,
    map: &mut BTreeMap<(VertexId, VertexId, GeomKey), ArrEdge>,
    u: VertexId,
    v: VertexId,
    edge: &InputEdge,
) {
    let (key_uv, dir) = if u <= v { ((u, v), 1) } else { ((v, u), -1) };
    // Canonical geometry is the geometry as traversed in the (a ≤ b) direction.
    let geom = match edge.geom {
        EdgeGeom::Seg => EdgeGeom::Seg,
        EdgeGeom::Arc {
            center,
            radius,
            ccw,
        } => EdgeGeom::Arc {
            center,
            radius,
            // If we flipped to canonical order, the traversal direction flips too.
            ccw: if dir == 1 { ccw } else { !ccw },
        },
    };
    let pa = store.point(key_uv.0);
    let pb = store.point(key_uv.1);
    let gk = geom_key_for(pa, pb, &geom);
    let entry = map.entry((key_uv.0, key_uv.1, gk)).or_insert(ArrEdge {
        a: key_uv.0,
        b: key_uv.1,
        wind_a: 0,
        wind_b: 0,
        geom,
    });
    match edge.operand {
        Operand::A => entry.wind_a += dir,
        Operand::B => entry.wind_b += dir,
    }
}

/// Build the DCEL half-edge graph from undirected arrangement edges.
fn build_dcel(
    store: &VertexStore,
    edges: &[ArrEdge],
    tol: &Tol,
) -> Result<Vec<HalfEdge>, Poly2Error> {
    let mut halfs: Vec<HalfEdge> = Vec::with_capacity(edges.len() * 2);
    let mut outgoing: HashMap<VertexId, Vec<usize>> = HashMap::new();

    for e in edges {
        let h0 = halfs.len();
        let h1 = h0 + 1;
        // h0: a → b with canonical geometry; h1: b → a with reversed geometry.
        let rev_geom = match e.geom {
            EdgeGeom::Seg => EdgeGeom::Seg,
            EdgeGeom::Arc {
                center,
                radius,
                ccw,
            } => EdgeGeom::Arc {
                center,
                radius,
                ccw: !ccw,
            },
        };
        halfs.push(HalfEdge {
            origin: e.a,
            dest: e.b,
            twin: h1,
            next: usize::MAX,
            geom: e.geom,
            face: None,
        });
        halfs.push(HalfEdge {
            origin: e.b,
            dest: e.a,
            twin: h0,
            next: usize::MAX,
            geom: rev_geom,
            face: None,
        });
        outgoing.entry(e.a).or_default().push(h0);
        outgoing.entry(e.b).or_default().push(h1);
    }

    // Sort outgoing half-edges by *tangent* angle at the vertex (so an arc orders
    // by the direction it leaves the vertex, not by its far endpoint), breaking
    // ties between collinear-tangent edges by signed curvature so a tangent pinch
    // (an arc that leaves a shared vertex collinear with a segment or another arc)
    // gets a well-defined ring order instead of an arbitrary one.
    for (&v, outs) in &mut outgoing {
        let vp = store.point(v);
        outs.sort_by(|&x, &y| {
            let kx = leave_key(store, &halfs[x], vp);
            let ky = leave_key(store, &halfs[y], vp);
            kx.cmp_ring(&ky, tol)
        });
    }

    let n = halfs.len();
    #[allow(clippy::needless_range_loop)]
    for h_in in 0..n {
        let v = halfs[h_in].dest;
        let twin = halfs[h_in].twin;
        let outs = outgoing.get(&v).ok_or(Poly2Error::Internal {
            what: "vertex missing from outgoing map",
        })?;
        let pos = outs
            .iter()
            .position(|&e| e == twin)
            .ok_or(Poly2Error::Internal {
                what: "twin not found in vertex ring",
            })?;
        let k = outs.len();
        let next_out = outs[(pos + k - 1) % k];
        halfs[h_in].next = next_out;
    }

    Ok(halfs)
}

/// Ring-ordering key for a half-edge leaving its origin vertex: the **tangent
/// angle** (the direction it departs) plus a **signed curvature** tie-break.
///
/// Two half-edges that depart a shared vertex with the *same* tangent direction
/// (a tangent pinch — e.g. an arc tangent to a segment, or two circles tangent at
/// a vertex) are indistinguishable by tangent angle alone, which leaves the DCEL
/// ring order undefined and makes face tracing drop a bounded face. The signed
/// curvature `kappa` resolves the tie: an arc bends toward its centre, so just
/// past the vertex it deviates from the common tangent line to the side the
/// centre lies on. `kappa > 0` means it veers to the **left** of the tangent
/// (toward higher angle), `kappa < 0` to the right; a segment has `kappa = 0`.
/// Ordering ties by ascending `kappa` therefore matches the ascending-angle
/// (CCW) order of the actual departing curves, and larger `|kappa|` (sharper
/// curve) orders further from the straight edge when several curves share both
/// tangent and bend side.
#[derive(Debug, Clone, Copy)]
struct LeaveKey {
    angle: f64,
    kappa: f64,
}

impl LeaveKey {
    /// Compare two leave keys for the outgoing ring ordering. Primary key is the
    /// tangent angle; when the two tangents are parallel and same-sense (within
    /// the angular tolerance) the signed curvature breaks the tie.
    fn cmp_ring(&self, other: &Self, tol: &Tol) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        // Tangent directions that differ by more than the angular tolerance are
        // ordered by angle alone. Genuine distinct directions in this (building-
        // scale, mostly axis-aligned) domain are O(1) radians apart, far above
        // `tol.angular`; only a true tangent pinch falls into the curvature branch.
        let da = self.angle - other.angle;
        if da.abs() > tol.angular {
            return self
                .angle
                .partial_cmp(&other.angle)
                .unwrap_or(Ordering::Equal);
        }
        self.kappa
            .partial_cmp(&other.kappa)
            .unwrap_or(Ordering::Equal)
    }
}

/// The [`LeaveKey`] (tangent angle + signed curvature) of a half-edge departing
/// `vp`, so arcs and straight edges order consistently in the DCEL ring even at a
/// tangent pinch.
fn leave_key(store: &VertexStore, h: &HalfEdge, vp: Point2) -> LeaveKey {
    match h.geom {
        EdgeGeom::Seg => {
            let dir = vp.to(store.point(h.dest));
            LeaveKey {
                angle: dir.y.atan2(dir.x),
                kappa: 0.0,
            }
        }
        EdgeGeom::Arc {
            center,
            radius,
            ccw,
        } => {
            let rad = center.to(vp); // outward radial (centre → vp)
                                     // CCW tangent is the radial rotated +90°; CW is −90°.
            let dir = if ccw {
                Vec2::new(-rad.y, rad.x)
            } else {
                Vec2::new(rad.y, -rad.x)
            };
            // Unit inward normal (vp → centre). The arc curves toward the centre,
            // so the sign of cross(tangent, inward_normal) gives the bend side:
            // `+1` when the centre lies to the left of the tangent. Scaling by the
            // curvature `1/radius` lets sharper curves order further from straight.
            let inward = vp.to(center);
            let inward_len = inward.len();
            let kappa = if radius > 0.0 && inward_len > 0.0 {
                let dir_len = dir.len().max(f64::MIN_POSITIVE);
                let side = (dir.x * inward.y - dir.y * inward.x) / (dir_len * inward_len);
                side / radius
            } else {
                0.0
            };
            LeaveKey {
                angle: dir.y.atan2(dir.x),
                kappa,
            }
        }
    }
}

/// Trace all face loops by following `next` pointers.
fn trace_faces(store: &VertexStore, mut halfs: Vec<HalfEdge>) -> Vec<FaceLoop> {
    let mut faces: Vec<FaceLoop> = Vec::new();
    let n = halfs.len();
    for start in 0..n {
        if halfs[start].face.is_some() || halfs[start].next == usize::MAX {
            continue;
        }
        let face_id = faces.len();
        let mut edges: Vec<LoopEdge> = Vec::new();
        let mut verts: Vec<Point2> = Vec::new();
        let mut cur = start;
        let cap = n + 1;
        let mut steps = 0;
        loop {
            halfs[cur].face = Some(face_id);
            let h = halfs[cur];
            verts.push(store.point(h.origin));
            edges.push(LoopEdge {
                a: h.origin,
                b: h.dest,
                geom: h.geom,
            });
            cur = halfs[cur].next;
            steps += 1;
            if cur == start || cur == usize::MAX || steps > cap {
                break;
            }
        }
        faces.push(FaceLoop {
            edges,
            vertices: verts,
        });
    }
    faces
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean::poly2d::geom::Point2;

    fn square(cx: f64, cy: f64, s: f64) -> Region {
        Region::from_points(&[
            Point2::new(cx, cy),
            Point2::new(cx + s, cy),
            Point2::new(cx + s, cy + s),
            Point2::new(cx, cy + s),
        ])
    }

    #[test]
    fn single_square_has_inner_and_outer_loop() {
        let tol = Tol::default();
        let a = square(0.0_f64, 0.0_f64, 1.0_f64);
        let b = Region::empty();
        let arr = Arrangement::build(&a, &b, &tol).unwrap();
        let pos = arr.faces.iter().filter(|f| f.signed_area() > 0.0).count();
        let neg = arr.faces.iter().filter(|f| f.signed_area() < 0.0).count();
        assert_eq!((pos, neg), (1, 1));
    }

    #[test]
    fn overlapping_squares_have_three_bounded_faces() {
        let tol = Tol::default();
        let a = square(0.0_f64, 0.0_f64, 1.0_f64);
        let b = square(0.5_f64, 0.5_f64, 1.0_f64);
        let arr = Arrangement::build(&a, &b, &tol).unwrap();
        let bounded = arr.faces.iter().filter(|f| f.signed_area() > 0.0).count();
        assert_eq!(bounded, 3);
    }

    #[test]
    fn winding_inside_and_outside() {
        let tol = Tol::default();
        let a = square(0.0_f64, 0.0_f64, 2.0_f64);
        let b = Region::empty();
        let arr = Arrangement::build(&a, &b, &tol).unwrap();
        assert_eq!(arr.winding(Point2::new(1.0_f64, 1.0_f64), Operand::A), 1);
        assert_eq!(arr.winding(Point2::new(5.0_f64, 5.0_f64), Operand::A), 0);
    }

    #[test]
    fn circle_winding_inside_outside() {
        let tol = Tol::default();
        let a = Region::circle(Point2::new(0.0_f64, 0.0_f64), 1.0_f64);
        let b = Region::empty();
        let arr = Arrangement::build(&a, &b, &tol).unwrap();
        assert_eq!(arr.winding(Point2::new(0.0_f64, 0.0_f64), Operand::A), 1);
        assert_eq!(arr.winding(Point2::new(5.0_f64, 0.0_f64), Operand::A), 0);
    }
}
