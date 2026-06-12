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

use std::collections::HashMap;

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

/// An atomic cell of the arrangement: a traced bounded face and its residency.
#[derive(Debug, Clone)]
pub(super) struct Cell {
    /// Vertex ids of the cell's boundary loop, in CCW order.
    pub vertex_ids: Vec<VertexId>,
    /// `inside[i]` is `true` if an interior point of the cell lies inside the
    /// `i`-th operand region.
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

        // Dedup the split fragments into undirected arrangement edges.
        let mut edge_set: HashMap<(VertexId, VertexId), ()> = HashMap::new();
        for (i, seg) in inputs.iter().enumerate() {
            let chain = ordered_chain(&store, seg, &split_points[i]);
            for w in chain.windows(2) {
                let (u, v) = (w[0], w[1]);
                if u == v {
                    continue;
                }
                let key = if u <= v { (u, v) } else { (v, u) };
                edge_set.insert(key, ());
            }
        }
        let undirected: Vec<(VertexId, VertexId)> = edge_set.into_keys().collect();

        // Build the DCEL and trace its cells, recording cell adjacency per edge.
        let (cells_raw, edges) = trace(&store, &undirected)?;

        // Tag each cell with its residency in every operand by exact winding at
        // an interior point.
        let n_ops = regions.len();
        let mut cells = Vec::with_capacity(cells_raw.len());
        for verts in cells_raw {
            let sample = face_sample_point(&store, &verts);
            let inside: Vec<bool> = (0..n_ops)
                .map(|op| winding(&store, &inputs, op, sample) != 0)
                .collect();
            cells.push(Cell {
                vertex_ids: verts,
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

/// Build the DCEL from undirected edges, trace cells, and return the cells'
/// CCW vertex loops plus the edge adjacency.
fn trace(
    store: &VertexStore,
    undirected: &[(VertexId, VertexId)],
) -> Result<(Vec<Vec<VertexId>>, Vec<ArrEdge>), PrismError> {
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

    // Trace every directed loop; a positive-area loop bounds an atomic cell.
    let mut cells: Vec<Vec<VertexId>> = Vec::new();
    // half-edge index -> cell index for the cell on its left, if bounded.
    let mut half_cell: Vec<Option<usize>> = vec![None; n];
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
        if signed_area(&pts) > 0.0 {
            let cid = cells.len();
            for &h in &chain {
                half_cell[h] = Some(cid);
            }
            cells.push(ring);
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
    let step = 1e-4_f64 * len;
    Point2::new((a.x + b.x) * 0.5 + nx * step, (a.y + b.y) * 0.5 + ny * step)
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
