//! The shared global 2-D arrangement over both operands' cross-sections.
//!
//! This is the robustness keystone of the prismatic build (`DESIGN.md` §4.2,
//! "区間×2D アレンジメント一括構築"). Both prism operands are constant along the
//! shared direction `d`, so every band's 2-D region boundary is a subset of the
//! two operand boundaries `∂R_a ∪ ∂R_b`. We therefore overlay **both** operand
//! boundaries into **one** planar subdivision, split it **once** at all mutual
//! intersections, and reuse that single segmentation for every band and every
//! axial level.
//!
//! The arrangement yields:
//!
//! * a set of **atomic cells** (the minimal 2-D faces of the subdivision), each
//!   tagged with its residency `(in_a, in_b)` taken at an interior sample point
//!   via exact winding;
//! * a set of **atomic edges** (the minimal segments of the subdivision), each
//!   carrying the two cells it borders (left / right).
//!
//! Walls are generated per atomic edge per band and interface caps per atomic
//! cell per level. Because both read the *same* vertices, a wall's vertical edge
//! and an interface's boundary edge meet at identical coordinates, so the
//! 3-D vertex / curve interning in [`super::build`] pairs every sibling exactly
//! — watertightness by construction, with over-splitting tolerated.
//!
//! The snapping, segment intersection and exact `orient2d` are reused verbatim
//! from the proven 2-D engine ([`crate::boolean::poly2d`]), so this module only
//! adds the DCEL trace and the cell adjacency it needs.

use std::collections::{BTreeSet, HashMap};

use crate::boolean::poly2d::geom::{orient2d, Edge2, Orient, Point2};
use crate::boolean::poly2d::intersect::intersect;
use crate::boolean::poly2d::snap::{VertexId, VertexStore};
use crate::boolean::poly2d::{Poly2Error, Region};
use crate::tolerance::Tol;

use super::error::PrismError;

/// A directed input boundary segment after snapping its endpoints, tagged with
/// the index of the operand region it came from.
#[derive(Debug, Clone, Copy)]
struct InputSeg {
    a: VertexId,
    b: VertexId,
    operand: usize,
}

/// An atomic cell of the arrangement: a traced bounded face (with any inner hole
/// rings) and its residency.
#[derive(Debug, Clone)]
pub(super) struct Cell {
    /// Vertex ids of the cell's outer boundary loop, in CCW order.
    pub vertex_ids: Vec<VertexId>,
    /// Inner hole rings (CW), each a loop of vertex ids. A cell that strictly
    /// contains another operand's cross-section (a slab with an interior shaft,
    /// a column piercing a slab) is an *annulus*: its outer ring is the
    /// enclosing boundary and each inner ring is a contained void. Without this,
    /// the contained tool was silently ignored and a phantom solid emitted
    /// (`DESIGN.md` §4.2).
    pub inner_rings: Vec<Vec<VertexId>>,
    /// `inside[i]` is `true` if an interior point of the cell (in the annulus,
    /// not in any hole) lies inside the `i`-th operand region.
    pub inside: Vec<bool>,
}

/// An atomic undirected edge of the arrangement, with its two bordering cells.
///
/// `left` is the cell on the left of the directed edge `a → b`; `right` is the
/// cell on its right. Either may be `None` for the unbounded exterior.
#[derive(Debug, Clone, Copy)]
pub(super) struct ArrEdge {
    pub a: VertexId,
    pub b: VertexId,
    pub left: Option<usize>,
    pub right: Option<usize>,
}

/// The built arrangement.
pub(super) struct Arrangement {
    store: VertexStore,
    /// The atomic cells (bounded faces).
    pub cells: Vec<Cell>,
    /// The atomic edges with cell adjacency.
    pub edges: Vec<ArrEdge>,
}

impl Arrangement {
    /// 2-D coordinate of a vertex id.
    #[inline]
    pub(super) fn point(&self, v: VertexId) -> Point2 {
        self.store.point(v)
    }

    /// Build the arrangement overlaying every operand region.
    pub(super) fn build(regions: &[&Region], tol: &Tol) -> Result<Self, PrismError> {
        let mut store = VertexStore::new(*tol);
        let mut inputs: Vec<InputSeg> = Vec::new();
        for (i, r) in regions.iter().enumerate() {
            ingest(r, i, &mut store, &mut inputs)?;
        }

        // Split every input segment at every mutual / self intersection.
        let n = inputs.len();
        let mut split_points: Vec<Vec<VertexId>> = vec![Vec::new(); n];
        for i in 0..n {
            let ei = seg_edge(&store, &inputs[i]);
            for j in (i + 1)..n {
                let ej = seg_edge(&store, &inputs[j]);
                let cr = intersect(&ei, &ej, tol).map_err(map_poly2)?;
                for p in cr.points {
                    let v = store.insert(p);
                    push_split(&mut split_points[i], v, &inputs[i]);
                    push_split(&mut split_points[j], v, &inputs[j]);
                }
            }
        }

        // Dedup the split fragments into undirected arrangement edges. A
        // `BTreeSet` keeps the edge order deterministic (sorted by vertex-id
        // pair) so the DCEL trace, cell order, and every downstream tie-break are
        // reproducible for identical input — a `HashMap` here would leak its
        // RandomState seed into the result.
        let mut edge_set: BTreeSet<(VertexId, VertexId)> = BTreeSet::new();
        for (i, seg) in inputs.iter().enumerate() {
            let chain = ordered_chain(&store, seg, &split_points[i]);
            for w in chain.windows(2) {
                let (u, v) = (w[0], w[1]);
                if u == v {
                    continue;
                }
                let key = if u <= v { (u, v) } else { (v, u) };
                edge_set.insert(key);
            }
        }
        let undirected: Vec<(VertexId, VertexId)> = edge_set.into_iter().collect();

        // Build the DCEL and trace its cells, recording cell adjacency per edge.
        let (cells_raw, edges) = trace(&store, &undirected)?;

        // Tag each cell with its residency in every operand by exact winding at
        // an interior point. The sample is taken in the annulus (the outer ring
        // minus the holes); `face_sample_point` offsets just inside the outer
        // boundary, which is always annulus material, so a contained hole never
        // contaminates the residency of its parent cell.
        let n_ops = regions.len();
        let mut cells = Vec::with_capacity(cells_raw.len());
        for raw in cells_raw {
            let sample = face_sample_point(&store, &raw.outer);
            let inside: Vec<bool> = (0..n_ops)
                .map(|op| winding(&store, &inputs, op, sample) != 0)
                .collect();
            cells.push(Cell {
                vertex_ids: raw.outer,
                inner_rings: raw.inners,
                inside,
            });
        }

        Ok(Self {
            store,
            cells,
            edges,
        })
    }
}

/// Map a 2-D engine error into a prismatic error (arc → Phase 3c).
fn map_poly2(e: Poly2Error) -> PrismError {
    PrismError::from(e)
}

/// Snap every contour edge endpoint and emit directed input segments.
fn ingest(
    r: &Region,
    operand: usize,
    store: &mut VertexStore,
    inputs: &mut Vec<InputSeg>,
) -> Result<(), PrismError> {
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
                Edge2::Arc(_) => return Err(PrismError::ArcNotYetSupported),
            }
        }
    }
    Ok(())
}

#[inline]
fn seg_edge(store: &VertexStore, s: &InputSeg) -> Edge2 {
    Edge2::seg(store.point(s.a), store.point(s.b))
}

fn push_split(splits: &mut Vec<VertexId>, v: VertexId, s: &InputSeg) {
    if v != s.a && v != s.b && !splits.contains(&v) {
        splits.push(v);
    }
}

/// Order a segment's interior split vertices along the segment.
fn ordered_chain(store: &VertexStore, s: &InputSeg, splits: &[VertexId]) -> Vec<VertexId> {
    let a = store.point(s.a);
    let b = store.point(s.b);
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

/// A directed half-edge in the trace DCEL.
struct Half {
    origin: VertexId,
    dest: VertexId,
    twin: usize,
    next: usize,
    cell: Option<usize>,
}

/// A traced cell: its CCW outer ring plus any CW inner hole rings.
struct RawCell {
    outer: Vec<VertexId>,
    inners: Vec<Vec<VertexId>>,
}

/// Build the DCEL from undirected edges, trace cells, and return the cells'
/// outer/inner vertex loops plus the edge adjacency.
fn trace(
    store: &VertexStore,
    undirected: &[(VertexId, VertexId)],
) -> Result<(Vec<RawCell>, Vec<ArrEdge>), PrismError> {
    let mut halfs: Vec<Half> = Vec::with_capacity(undirected.len() * 2);
    let mut outgoing: HashMap<VertexId, Vec<usize>> = HashMap::new();
    for &(a, b) in undirected.iter() {
        let h0 = halfs.len();
        let h1 = h0 + 1;
        halfs.push(Half {
            origin: a,
            dest: b,
            twin: h1,
            next: usize::MAX,
            cell: None,
        });
        halfs.push(Half {
            origin: b,
            dest: a,
            twin: h0,
            next: usize::MAX,
            cell: None,
        });
        outgoing.entry(a).or_default().push(h0);
        outgoing.entry(b).or_default().push(h1);
    }

    for (&v, outs) in &mut outgoing {
        let vp = store.point(v);
        outs.sort_by(|&x, &y| {
            let ax = angle(store, &halfs[x], vp);
            let ay = angle(store, &halfs[y], vp);
            ax.partial_cmp(&ay).unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    let n = halfs.len();
    #[allow(clippy::needless_range_loop)]
    for h_in in 0..n {
        let v = halfs[h_in].dest;
        let twin = halfs[h_in].twin;
        let outs = outgoing
            .get(&v)
            .ok_or(PrismError::Poly2(Poly2Error::Internal {
                what: "vertex missing from arrangement ring",
            }))?;
        let pos = outs
            .iter()
            .position(|&e| e == twin)
            .ok_or(PrismError::Poly2(Poly2Error::Internal {
                what: "twin missing from arrangement ring",
            }))?;
        let k = outs.len();
        halfs[h_in].next = outs[(pos + k - 1) % k];
    }

    // Trace every directed loop. A positive-area (CCW) loop is the outer
    // boundary of a cell; a negative-area (CW) loop is an inner hole boundary
    // that must be attached to the cell whose interior contains it. Both kinds
    // of loop's half-edges must resolve to the same cell so the build's wall and
    // interface generation sees the annulus on the hole side (not `None`).
    struct Loop {
        ring: Vec<VertexId>,
        pts: Vec<Point2>,
        chain: Vec<usize>,
        area: f64,
    }
    let mut loops: Vec<Loop> = Vec::new();
    for start in 0..n {
        if halfs[start].cell.is_some() || halfs[start].next == usize::MAX {
            continue;
        }
        let mut ring: Vec<VertexId> = Vec::new();
        let mut pts: Vec<Point2> = Vec::new();
        let mut chain: Vec<usize> = Vec::new();
        let mut cur = start;
        let cap = n + 1;
        let mut steps = 0;
        loop {
            halfs[cur].cell = Some(usize::MAX); // mark visited
            ring.push(halfs[cur].origin);
            pts.push(store.point(halfs[cur].origin));
            chain.push(cur);
            cur = halfs[cur].next;
            steps += 1;
            if cur == start || cur == usize::MAX || steps > cap {
                break;
            }
        }
        let area = signed_area(&pts);
        loops.push(Loop {
            ring,
            pts,
            chain,
            area,
        });
    }

    // half-edge index -> cell index for the cell on its left, if bounded.
    let mut half_cell: Vec<Option<usize>> = vec![None; n];

    // Positive loops define cells (their index is their position in `cells`).
    let mut cells: Vec<RawCell> = Vec::new();
    let mut pos_loops: Vec<usize> = Vec::new();
    for (li, lp) in loops.iter().enumerate() {
        if lp.area > 0.0 {
            let cid = cells.len();
            for &h in &lp.chain {
                half_cell[h] = Some(cid);
            }
            cells.push(RawCell {
                outer: lp.ring.clone(),
                inners: Vec::new(),
            });
            pos_loops.push(li);
        }
    }

    // A sorted vertex-id key identifies a loop's boundary regardless of start /
    // orientation.
    let ring_key = |ring: &[VertexId]| {
        let mut k: Vec<usize> = ring.iter().map(|v| v.0).collect();
        k.sort_unstable();
        k
    };
    let pos_keys: Vec<Vec<usize>> = pos_loops
        .iter()
        .map(|&li| ring_key(&loops[li].ring))
        .collect();

    // Determine which negative loops are genuine **holes** and attach each to the
    // cell that encloses it.
    //
    // A genuine hole's boundary is exactly the reverse of some *positive* cell's
    // ring — the cell that fills the void (the shaft / inner column is itself an
    // atomic cell, traced CCW). So a negative loop is a hole iff its vertex key
    // matches a positive cell `V` (the void cell); the hole then belongs to the
    // smallest *other* cell whose outer ring encloses `V`'s interior (the
    // enclosing annulus). A negative loop matching *no* cell is a complex
    // unbounded wrap (the outline around protruding / multi-component material),
    // never a hole — skip it so its half-edges stay `None` (the exterior). This
    // is robust where geometric containment alone is not: a wrap's non-convex
    // outline could otherwise probe into an unrelated cell.
    for lp in &loops {
        if lp.area >= 0.0 {
            continue;
        }
        let lp_key = ring_key(&lp.ring);
        // The void cell `V`: a positive cell with this exact boundary.
        let Some(void_ci) = pos_keys.iter().position(|k| *k == lp_key) else {
            continue; // not any cell's boundary ⇒ unbounded wrap, not a hole
        };
        // Probe a point strictly inside the void (i.e. inside `V`), then find the
        // smallest enclosing cell other than `V`.
        let probe = hole_probe_point(&lp.pts);
        let mut parent: Option<usize> = None;
        let mut parent_area = f64::INFINITY;
        for (ci, &li) in pos_loops.iter().enumerate() {
            if ci == void_ci {
                continue;
            }
            let outer = &loops[li];
            if outer.area < parent_area && point_in_ring(&outer.pts, probe) {
                parent = Some(ci);
                parent_area = outer.area;
            }
        }
        if let Some(cid) = parent {
            cells[cid].inners.push(lp.ring.clone());
            for &h in &lp.chain {
                half_cell[h] = Some(cid);
            }
        }
    }

    // For each undirected edge, the two half-edges give its left / right cells.
    // Half-edge `a→b` has its traced cell on the left; the twin's cell is on
    // the right of `a→b`.
    let mut edges: Vec<ArrEdge> = Vec::with_capacity(undirected.len());
    for (ei, &(a, b)) in undirected.iter().enumerate() {
        // Find the half-edge oriented a→b and its twin.
        let h_ab = 2 * ei; // matches construction order (origin a, dest b)
        let h_ba = 2 * ei + 1;
        debug_assert_eq!(halfs[h_ab].origin, a);
        debug_assert_eq!(halfs[h_ab].dest, b);
        edges.push(ArrEdge {
            a,
            b,
            left: half_cell[h_ab],
            right: half_cell[h_ba],
        });
    }

    Ok((cells, edges))
}

/// A point strictly **inside** a CW loop (its enclosed void), used to find which
/// cell encloses the hole: the parent is the smallest cell whose outer ring
/// contains this interior point (excluding the cell whose own ring *is* this
/// loop). For a CW loop the interior is to the *right* of each directed edge, so
/// the right-normal offset of an edge midpoint points inward; the step is shrunk
/// until the point is verified inside the loop, so a thin or non-convex loop is
/// still sampled correctly.
fn hole_probe_point(pts: &[Point2]) -> Point2 {
    let n = pts.len();
    if n < 3 {
        return pts.first().copied().unwrap_or(Point2::new(0.0, 0.0));
    }
    // Longest edge for a stable normal.
    let mut best_i = 0usize;
    let mut best_len = -1.0_f64;
    for i in 0..n {
        let l = pts[i].dist(pts[(i + 1) % n]);
        if l > best_len {
            best_len = l;
            best_i = i;
        }
    }
    let a = pts[best_i];
    let b = pts[(best_i + 1) % n];
    let d = a.to(b);
    let len = d.len().max(f64::MIN_POSITIVE);
    // Right normal of a→b is (dy, -dx)/len; for a CW loop this points inward
    // (into the enclosed void).
    let nx = d.y / len;
    let ny = -d.x / len;
    let mx = (a.x + b.x) * 0.5;
    let my = (a.y + b.y) * 0.5;
    let mut step = 1e-4_f64 * len;
    for _ in 0..60 {
        let p = Point2::new(mx + nx * step, my + ny * step);
        if point_in_ring(pts, p) {
            return p;
        }
        step *= 0.5;
    }
    // Fallback: centroid (correct for convex loops).
    let (mut cx, mut cy) = (0.0_f64, 0.0_f64);
    for p in pts {
        cx += p.x;
        cy += p.y;
    }
    Point2::new(cx / n as f64, cy / n as f64)
}

/// `true` if `p` is strictly inside the polygon ring `pts`, by the exact
/// [`orient2d`] crossing rule.
fn point_in_ring(pts: &[Point2], p: Point2) -> bool {
    let n = pts.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let a = pts[i];
        let b = pts[j];
        if (a.y > p.y) != (b.y > p.y) {
            let upward = b.y > a.y;
            let left = orient2d(a, b, p) == Orient::Left;
            if left == upward {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

#[inline]
fn angle(store: &VertexStore, h: &Half, vp: Point2) -> f64 {
    let d = vp.to(store.point(h.dest));
    d.y.atan2(d.x)
}

/// Signed area of a 2-D loop (positive = CCW).
fn signed_area(pts: &[Point2]) -> f64 {
    let n = pts.len();
    if n < 3 {
        return 0.0;
    }
    let mut acc = 0.0_f64;
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        acc += a.x * b.y - b.x * a.y;
    }
    0.5 * acc
}

/// A point strictly inside the cell bounded by the CCW loop `ring`, taken just
/// inside the longest edge (the same stable rule the 2-D engine uses).
fn face_sample_point(store: &VertexStore, ring: &[VertexId]) -> Point2 {
    let n = ring.len();
    if n < 2 {
        return centroid(store, ring);
    }
    let mut best_i = 0usize;
    let mut best_len = -1.0_f64;
    for i in 0..n {
        let a = store.point(ring[i]);
        let b = store.point(ring[(i + 1) % n]);
        let l = a.dist(b);
        if l > best_len {
            best_len = l;
            best_i = i;
        }
    }
    let a = store.point(ring[best_i]);
    let b = store.point(ring[(best_i + 1) % n]);
    let d = a.to(b);
    let len = d.len();
    if len <= 0.0 {
        return centroid(store, ring);
    }
    let nx = -d.y / len;
    let ny = d.x / len;
    let mx = (a.x + b.x) * 0.5;
    let my = (a.y + b.y) * 0.5;
    // Every traced cell here is CCW (only positive-area loops become cells), so
    // the cell's own polygon *is* its interior. A fixed `1e-4 · edge` step
    // overshoots a thin cell (a 0.5 mm residual strip over a 10 m wall, step
    // 1 mm) and lands in the neighbouring cell, dropping legitimate geometry.
    // Verify the candidate is inside the cell with the exact predicate and
    // shrink the step until it is.
    let pts: Vec<Point2> = ring.iter().map(|&v| store.point(v)).collect();
    let mut step = 1e-4_f64 * len;
    for _ in 0..60 {
        let p = Point2::new(mx + nx * step, my + ny * step);
        if point_in_cell(&pts, p) {
            return p;
        }
        step *= 0.5;
    }
    centroid(store, ring)
}

/// `true` if `p` is strictly inside the simple polygon `ring`, decided with the
/// exact [`orient2d`] predicate for every edge crossing.
fn point_in_cell(ring: &[Point2], p: Point2) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let a = ring[i];
        let b = ring[j];
        if (a.y > p.y) != (b.y > p.y) {
            let upward = b.y > a.y;
            let left = orient2d(a, b, p) == Orient::Left;
            if left == upward {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

fn centroid(store: &VertexStore, ring: &[VertexId]) -> Point2 {
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    for &v in ring {
        let p = store.point(v);
        cx += p.x;
        cy += p.y;
    }
    let k = ring.len().max(1) as f64;
    Point2::new(cx / k, cy / k)
}

/// Winding number of `p` w.r.t. the given operand's original directed segments.
fn winding(store: &VertexStore, inputs: &[InputSeg], operand: usize, p: Point2) -> i32 {
    let mut w = 0_i32;
    for s in inputs.iter().filter(|s| s.operand == operand) {
        let a = store.point(s.a);
        let b = store.point(s.b);
        if a.y <= p.y {
            if b.y > p.y && orient2d(a, b, p) == Orient::Left {
                w += 1;
            }
        } else if b.y <= p.y && orient2d(a, b, p) == Orient::Right {
            w -= 1;
        }
    }
    w
}
