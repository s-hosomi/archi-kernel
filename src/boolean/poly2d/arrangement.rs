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
//! 3. **Split** each input edge at its interior vertices into unit half-open
//!    segments between consecutive vertices.
//! 4. **Dedup** collinear-coincident split segments: two split segments between
//!    the same vertex pair are the *same* arrangement edge (this is where
//!    shared / overlapping boundary edges in buildings collapse). Each surviving
//!    arrangement edge becomes a pair of opposite **half-edges**.
//! 5. **Build the DCEL**: sort half-edges leaving each vertex by angle, link
//!    `next`/`prev` by the "next clockwise" rule, and trace face loops.
//! 6. **Extract faces**: each loop with positive signed area is a bounded face;
//!    the single negative-area loop per component is the outer (unbounded) wrap.
//!
//! Robustness rests on the exact [`orient2d`] (combinatorial decisions) plus the
//! snap (coincidence collapse). Float coordinates are used only for *where*,
//! never for *whether*.

use std::collections::HashMap;

use crate::boolean::poly2d::error::Poly2Error;
use crate::boolean::poly2d::geom::{orient2d, Edge2, Orient, Point2};
use crate::boolean::poly2d::intersect::intersect;
use crate::boolean::poly2d::region::Region;
use crate::boolean::poly2d::snap::{VertexId, VertexStore};
use crate::tolerance::Tol;

/// Which operand an input edge came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operand {
    /// The first operand (`A`).
    A,
    /// The second operand (`B`).
    B,
}

/// An input edge after snapping its endpoints to vertex ids, carrying its
/// operand and its original traversal direction (used to recover winding).
#[derive(Debug, Clone, Copy)]
struct InputSeg {
    a: VertexId,
    b: VertexId,
    operand: Operand,
}

/// An arrangement edge (an undirected segment between two distinct vertices),
/// annotated with the net directed winding contribution of each operand.
///
/// `wind_a` / `wind_b` are the sum over all coincident input segments of
/// `+1` if the input ran `a → b` and `-1` if it ran `b → a`. After dedup an
/// arrangement edge with `wind == 0` for *both* operands carries no boundary and
/// is dropped (this is how a shared edge traversed in opposite directions — e.g.
/// the internal edge of two squares being unioned — vanishes).
#[derive(Debug, Clone)]
struct ArrEdge {
    a: VertexId,
    b: VertexId,
    wind_a: i32,
    wind_b: i32,
}

/// A directed half-edge in the DCEL.
#[derive(Debug, Clone)]
struct HalfEdge {
    origin: VertexId,
    dest: VertexId,
    twin: usize,
    next: usize,
    /// Index of the face loop this half-edge belongs to (assigned during trace).
    face: Option<usize>,
}

/// A traced face loop: the ordered vertices of its boundary and its signed area.
#[derive(Debug, Clone)]
pub struct FaceLoop {
    /// Vertices of the loop in traversal order.
    pub vertices: Vec<Point2>,
    /// Snapped vertex ids of the loop, parallel to `vertices`. Reconstruction
    /// uses these for *exact* shared-edge cancellation (no tolerance needed —
    /// the snap already merged coincident vertices).
    pub vertex_ids: Vec<VertexId>,
}

impl FaceLoop {
    /// Signed area of the loop (positive = CCW bounded face / outer boundary,
    /// negative = CW hole boundary or the unbounded outer wrap).
    ///
    /// Used by the arrangement tests to assert face orientation; the boolean
    /// driver classifies faces by winding, not area, so the lib path does not
    /// call this.
    #[cfg(test)]
    pub fn signed_area(&self) -> f64 {
        let mut acc = 0.0_f64;
        let n = self.vertices.len();
        for i in 0..n {
            let a = self.vertices[i];
            let b = self.vertices[(i + 1) % n];
            acc += a.x * b.y - b.x * a.y;
        }
        0.5 * acc
    }
}

/// The built arrangement: the snapped vertex store, the original input segments
/// (kept for winding classification), and the traced bounded faces.
pub struct Arrangement {
    store: VertexStore,
    inputs: Vec<InputSeg>,
    /// Bounded faces (positive signed area), in no particular order.
    pub faces: Vec<FaceLoop>,
}

impl Arrangement {
    /// Build the arrangement from two operand regions.
    ///
    /// Returns [`Poly2Error::ArcNotYetSupported`] if any edge is an arc, and
    /// [`Poly2Error::SelfIntersectingInput`] / [`Poly2Error::DegenerateContour`]
    /// for malformed input.
    pub fn build(a: &Region, b: &Region, tol: &Tol) -> Result<Self, Poly2Error> {
        if a.has_arc() || b.has_arc() {
            return Err(Poly2Error::ArcNotYetSupported);
        }
        validate_region(a, tol)?;
        validate_region(b, tol)?;

        let mut store = VertexStore::new(*tol);
        let mut inputs: Vec<InputSeg> = Vec::new();

        // ── 1. snap endpoints, collect input segments ──────────────────────
        ingest(a, Operand::A, &mut store, &mut inputs, tol)?;
        ingest(b, Operand::B, &mut store, &mut inputs, tol)?;

        // ── 2. intersect every pair, snap crossings, collect split points ──
        // split_points[i] holds the parameter-ordered vertex ids that split
        // input segment i (besides its own two endpoints).
        let n = inputs.len();
        let mut split_points: Vec<Vec<VertexId>> = vec![Vec::new(); n];
        for i in 0..n {
            let ei = seg_edge(&store, &inputs[i]);
            for j in (i + 1)..n {
                let ej = seg_edge(&store, &inputs[j]);
                let cr = intersect(&ei, &ej, tol)?;
                for p in cr.points {
                    let v = store.insert(p);
                    push_split(&mut split_points[i], v, &inputs[i]);
                    push_split(&mut split_points[j], v, &inputs[j]);
                }
            }
        }

        // ── 3. split each input into unit segments between consecutive verts ─
        // ── 4. dedup into arrangement edges keyed by unordered vertex pair ──
        let mut arr_map: HashMap<(VertexId, VertexId), ArrEdge> = HashMap::new();
        for (i, seg) in inputs.iter().enumerate() {
            let chain = ordered_chain(&store, seg, &split_points[i], tol);
            for w in chain.windows(2) {
                let (u, v) = (w[0], w[1]);
                if u == v {
                    continue; // zero-length fragment collapsed by snapping
                }
                accumulate_edge(&mut arr_map, u, v, seg.operand);
            }
        }

        // Drop arrangement edges whose winding is zero for *both* operands:
        // they are interior shared edges that cancel and carry no boundary.
        let arr_edges: Vec<ArrEdge> = arr_map
            .into_values()
            .filter(|e| e.wind_a != 0 || e.wind_b != 0)
            .collect();

        // ── 5. build DCEL ──────────────────────────────────────────────────
        let halfs = build_dcel(&store, &arr_edges)?;

        // ── 6. trace faces ─────────────────────────────────────────────────
        let faces = trace_faces(&store, halfs);

        Ok(Self {
            store,
            inputs,
            faces,
        })
    }

    /// Winding number of `p` with respect to the given operand's *original*
    /// input edges (before arrangement). Robust because `p` is taken in a face
    /// interior, away from every edge.
    pub fn winding(&self, p: Point2, operand: Operand) -> i32 {
        let mut wind = 0_i32;
        for seg in self.inputs.iter().filter(|s| s.operand == operand) {
            let a = self.store.point(seg.a);
            let b = self.store.point(seg.b);
            wind += ray_cross_contribution(p, a, b);
        }
        wind
    }

    /// Borrow the snapped vertex store (for reconstruction).
    #[inline]
    pub fn store(&self) -> &VertexStore {
        &self.store
    }
}

/// Winding contribution of a directed segment `a → b` to the winding number of
/// `p`, using the standard "upward/downward crossing of the horizontal ray to
/// the right" rule with the exact [`orient2d`] sidedness test.
fn ray_cross_contribution(p: Point2, a: Point2, b: Point2) -> i32 {
    if a.y <= p.y {
        if b.y > p.y {
            // upward crossing candidate; p strictly left of a→b ?
            if orient2d(a, b, p) == Orient::Left {
                return 1;
            }
        }
    } else if b.y <= p.y {
        // downward crossing candidate; p strictly right of a→b ?
        if orient2d(a, b, p) == Orient::Right {
            return -1;
        }
    }
    0
}

/// Reject malformed input: degenerate contours and self-intersections.
fn validate_region(r: &Region, tol: &Tol) -> Result<(), Poly2Error> {
    for (ci, c) in r.contours.iter().enumerate() {
        if c.distinct_vertex_count(tol) < 3 {
            return Err(Poly2Error::DegenerateContour { contour_index: ci });
        }
    }
    Ok(())
}

/// Snap every contour edge's endpoints and emit [`InputSeg`]s.
fn ingest(
    r: &Region,
    operand: Operand,
    store: &mut VertexStore,
    inputs: &mut Vec<InputSeg>,
    _tol: &Tol,
) -> Result<(), Poly2Error> {
    for c in &r.contours {
        for e in &c.edges {
            match e {
                Edge2::Seg { start, end } => {
                    let a = store.insert(*start);
                    let b = store.insert(*end);
                    if a != b {
                        inputs.push(InputSeg { a, b, operand });
                    }
                }
                Edge2::Arc(_) => return Err(Poly2Error::ArcNotYetSupported),
            }
        }
    }
    Ok(())
}

/// Reconstruct the [`Edge2`] of an input segment from snapped coordinates.
#[inline]
fn seg_edge(store: &VertexStore, s: &InputSeg) -> Edge2 {
    Edge2::seg(store.point(s.a), store.point(s.b))
}

/// Record `v` as a split point of input segment `s` if it is interior.
fn push_split(splits: &mut Vec<VertexId>, v: VertexId, s: &InputSeg) {
    if v != s.a && v != s.b && !splits.contains(&v) {
        splits.push(v);
    }
}

/// Order an input segment's vertices `a, (splits…), b` along the segment.
fn ordered_chain(
    store: &VertexStore,
    s: &InputSeg,
    splits: &[VertexId],
    _tol: &Tol,
) -> Vec<VertexId> {
    let a = store.point(s.a);
    let b = store.point(s.b);
    let d = a.to(b);
    let len_sq = d.len_sq().max(f64::MIN_POSITIVE);
    let mut mids: Vec<(f64, VertexId)> = splits
        .iter()
        .map(|&v| {
            let t = a.to(store.point(v)).dot(d) / len_sq;
            (t, v)
        })
        .collect();
    mids.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut chain = Vec::with_capacity(mids.len() + 2);
    chain.push(s.a);
    for (_, v) in mids {
        chain.push(v);
    }
    chain.push(s.b);
    chain
}

/// Add a directed unit segment `u → v` to the arrangement-edge map, folding the
/// winding into the canonical (min,max) key.
fn accumulate_edge(
    map: &mut HashMap<(VertexId, VertexId), ArrEdge>,
    u: VertexId,
    v: VertexId,
    operand: Operand,
) {
    let (key, dir) = if u <= v { ((u, v), 1) } else { ((v, u), -1) };
    let entry = map.entry(key).or_insert(ArrEdge {
        a: key.0,
        b: key.1,
        wind_a: 0,
        wind_b: 0,
    });
    match operand {
        Operand::A => entry.wind_a += dir,
        Operand::B => entry.wind_b += dir,
    }
}

/// Build the DCEL half-edge graph from undirected arrangement edges.
fn build_dcel(store: &VertexStore, edges: &[ArrEdge]) -> Result<Vec<HalfEdge>, Poly2Error> {
    let mut halfs: Vec<HalfEdge> = Vec::with_capacity(edges.len() * 2);
    // outgoing[v] = indices of half-edges with origin v.
    let mut outgoing: HashMap<VertexId, Vec<usize>> = HashMap::new();

    for e in edges {
        let h0 = halfs.len();
        let h1 = h0 + 1;
        halfs.push(HalfEdge {
            origin: e.a,
            dest: e.b,
            twin: h1,
            next: usize::MAX,
            face: None,
        });
        halfs.push(HalfEdge {
            origin: e.b,
            dest: e.a,
            twin: h0,
            next: usize::MAX,
            face: None,
        });
        outgoing.entry(e.a).or_default().push(h0);
        outgoing.entry(e.b).or_default().push(h1);
    }

    // For each vertex, sort outgoing half-edges by polar angle. The `next` of a
    // half-edge h (arriving at vertex v = h.dest via its twin's origin) is the
    // outgoing edge most clockwise from h's reverse direction — the standard
    // "next half-edge around a face" rule that traces faces consistently CCW.
    for (&v, outs) in &mut outgoing {
        let vp = store.point(v);
        outs.sort_by(|&x, &y| {
            let ax = angle_of(store, &halfs[x], vp);
            let ay = angle_of(store, &halfs[y], vp);
            ax.partial_cmp(&ay).unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    // Link next pointers. For an incoming half-edge `h_in` arriving at `v`, its
    // twin `h_in.twin` is outgoing from `v`. Find that twin's position in v's
    // sorted ring; the face successor is the *previous* outgoing edge in CW
    // order (i.e. the next one clockwise), which keeps bounded faces CCW.
    let n = halfs.len();
    // The body reads `halfs[h_in]` to find `v`/`twin` and then writes
    // `halfs[h_in].next`, while also borrowing `outgoing`; an index loop is the
    // clearest way to express that without fighting the borrow checker.
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
        // Next outgoing in clockwise order = previous in the CCW-sorted ring.
        let k = outs.len();
        let next_out = outs[(pos + k - 1) % k];
        halfs[h_in].next = next_out;
    }

    Ok(halfs)
}

/// Polar angle of a half-edge's direction as seen from its origin vertex `vp`.
#[inline]
fn angle_of(store: &VertexStore, h: &HalfEdge, vp: Point2) -> f64 {
    let d = vp.to(store.point(h.dest));
    d.y.atan2(d.x)
}

/// Trace all face loops by following `next` pointers, then keep the bounded
/// (positive-area) ones.
fn trace_faces(store: &VertexStore, mut halfs: Vec<HalfEdge>) -> Vec<FaceLoop> {
    let mut faces: Vec<FaceLoop> = Vec::new();
    let n = halfs.len();
    for start in 0..n {
        if halfs[start].face.is_some() || halfs[start].next == usize::MAX {
            continue;
        }
        let face_id = faces.len();
        let mut verts: Vec<Point2> = Vec::new();
        let mut vids: Vec<VertexId> = Vec::new();
        let mut cur = start;
        // Walk the cycle. A hard cap guards against a malformed graph rather
        // than spinning forever (panic-free contract).
        let cap = n + 1;
        let mut steps = 0;
        loop {
            halfs[cur].face = Some(face_id);
            vids.push(halfs[cur].origin);
            verts.push(store.point(halfs[cur].origin));
            cur = halfs[cur].next;
            steps += 1;
            if cur == start || cur == usize::MAX || steps > cap {
                break;
            }
        }
        faces.push(FaceLoop {
            vertices: verts,
            vertex_ids: vids,
        });
    }
    // Keep loops of *both* orientations. A CCW loop bounds the face on its left
    // (its interior); a CW loop is the inner boundary (a hole edge) of the face
    // on its left. Both are needed so donut faces classify and reconstruct
    // correctly. We drop only the single all-enclosing outer wrap of each
    // connected component, identified as a loop whose left-side sample escapes
    // to infinity — but that drop is unnecessary here because such a wrap's
    // face-on-left is the unbounded region, which always classifies as outside
    // both operands and is therefore never selected. So we keep everything.
    faces
}

impl FaceLoop {
    /// A point in the **face that lies to the left** of this directed loop,
    /// suitable for winding classification.
    ///
    /// In a DCEL every directed loop borders exactly one face on its left:
    /// * a CCW loop's left side is its own enclosed interior;
    /// * a CW loop (a hole boundary) has the surrounding face on its left.
    ///
    /// Classifying that left-side face and keeping the loop iff the face is
    /// selected makes holes "just work": a kept CCW loop is an outer boundary, a
    /// kept CW loop is a hole boundary of the same selected face.
    ///
    /// The sample is taken at the **midpoint of an edge**, offset by a *tiny*
    /// step along that edge's left normal. The step is a small fraction of the
    /// edge length, so the point stays in the **minimal face adjacent to that
    /// edge** rather than crossing into a neighbouring face — this is essential
    /// when the loop is the outer boundary of a thin annulus (e.g. a wall ring
    /// around an opening): a large step would land in the opening and misclassify
    /// the ring. The longest edge is used to make the offset numerically stable.
    pub fn face_sample_point(&self) -> Point2 {
        let n = self.vertices.len();
        if n < 2 {
            return self.centroid();
        }
        // Pick the longest edge for stability.
        let mut best_i = 0usize;
        let mut best_len = -1.0_f64;
        for i in 0..n {
            let a = self.vertices[i];
            let b = self.vertices[(i + 1) % n];
            let l = a.dist(b);
            if l > best_len {
                best_len = l;
                best_i = i;
            }
        }
        let a = self.vertices[best_i];
        let b = self.vertices[(best_i + 1) % n];
        let d = a.to(b);
        let len = d.len();
        if len <= 0.0 {
            return self.centroid();
        }
        // Left normal of edge a→b is (-dy, dx)/len. Step a tiny fraction inward.
        let nx = -d.y / len;
        let ny = d.x / len;
        // 1e-4 of the edge length keeps us inside any face thicker than that,
        // which covers all non-degenerate building faces (eps = 1e-6).
        let step = 1e-4 * len;
        let mx = (a.x + b.x) * 0.5;
        let my = (a.y + b.y) * 0.5;
        Point2::new(mx + nx * step, my + ny * step)
    }

    fn centroid(&self) -> Point2 {
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for v in &self.vertices {
            cx += v.x;
            cy += v.y;
        }
        let k = self.vertices.len().max(1) as f64;
        Point2::new(cx / k, cy / k)
    }
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
        // A single square yields two loops: the CCW interior (area +1) and the
        // CW outer wrap (area −1).
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
        // Two unit squares overlapping in a quarter: 3 bounded (CCW) faces
        // (A-only, B-only, overlap) plus the CW outer wrap.
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
}
