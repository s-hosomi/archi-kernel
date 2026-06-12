//! Reconstruct a [`Region`] from the set of selected arrangement faces.
//!
//! After classification we have a set of bounded faces (each a CCW loop of
//! vertex ids) that should be kept. Their union is the result. Reconstruction:
//!
//! 1. **Cancel internal edges.** Each selected face contributes its CCW boundary
//!    as directed edges. An edge shared by two selected faces appears once in
//!    each direction (`u→v` and `v→u`) and cancels — it is interior to the
//!    union, not part of the result boundary. This is what makes the internal
//!    edge of two squares being unioned *disappear*.
//! 2. **Trace boundary loops.** The surviving directed edges form closed loops.
//!    At each vertex we leave along the edge that makes the tightest left turn
//!    from the reverse of the incoming direction, which keeps the trace hugging
//!    the boundary so outer loops come out CCW and holes CW.
//! 3. **Normalize orientation.** Outers are emitted CCW, holes CW (the trace
//!    sign already encodes this); zero-area slivers are dropped.
//!
//! Because step 1 cancels shared edges by **exact vertex id** (the snap already
//! merged coincident vertices), adjacency needs no tolerance here.

use std::collections::HashMap;

use crate::boolean::poly2d::geom::{Point2, Vec2};
use crate::boolean::poly2d::region::{Contour, Region};
use crate::boolean::poly2d::snap::{VertexId, VertexStore};
use crate::tolerance::Tol;

/// A directed edge between two snapped vertices.
type DEdge = (VertexId, VertexId);

/// Rebuild a region from selected faces, each given as a CCW loop of vertex ids.
pub fn reconstruct(store: &VertexStore, selected_loops: &[Vec<VertexId>], tol: &Tol) -> Region {
    // ── 1. cancel internal edges ───────────────────────────────────────────
    let mut count: HashMap<DEdge, i32> = HashMap::new();
    for loop_v in selected_loops {
        let n = loop_v.len();
        for i in 0..n {
            let u = loop_v[i];
            let v = loop_v[(i + 1) % n];
            if u == v {
                continue;
            }
            *count.entry((u, v)).or_insert(0) += 1;
        }
    }

    // Surviving directed edges: those whose count exceeds their reverse count.
    let mut adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();
    let keys: Vec<DEdge> = count.keys().copied().collect();
    for (u, v) in keys {
        let fwd = *count.get(&(u, v)).unwrap_or(&0);
        let rev = *count.get(&(v, u)).unwrap_or(&0);
        let net = fwd - rev;
        if net > 0 {
            for _ in 0..net {
                adj.entry(u).or_default().push(v);
            }
        }
    }
    if adj.is_empty() {
        return Region::empty();
    }

    // ── 2. trace boundary loops ─────────────────────────────────────────────
    let mut contours: Vec<Contour> = Vec::new();
    let mut starts: Vec<VertexId> = adj.keys().copied().collect();
    starts.sort();

    for start in starts {
        while adj.get(&start).map(|v| !v.is_empty()).unwrap_or(false) {
            if let Some(loop_ids) = trace_one_loop(store, &mut adj, start) {
                if loop_ids.len() >= 3 {
                    let pts: Vec<Point2> = loop_ids.iter().map(|&id| store.point(id)).collect();
                    contours.push(Contour::from_points(&pts));
                }
            } else {
                break;
            }
        }
    }

    // Regularization: drop sliver contours whose mean thickness is at or below
    // `eps`. A sliver's area is at most `eps × perimeter` (mean thickness =
    // 2·area / perimeter ≤ eps), which is the residue of near-coincident
    // vertices that landed just outside the pairwise merge radius. Genuine thin
    // features sit well above this bound (their vertices are spread far beyond
    // `eps`), so this removes snap artifacts without erasing real geometry.
    contours.retain(|c| !is_sliver(c, tol));
    Region::new(contours)
}

/// Multiple of `eps` below which a whole contour is treated as a snap-cluster
/// artifact. `eps`-merging is not transitively closed (A≈B, B≈C, A≉C), so a
/// handful of chained near-merges can leave a micro-contour whose every vertex
/// is within a few `eps` of the others — a meaningless residue, not geometry.
/// 64 is a safety factor over the worst-case chained slack; at the default
/// `eps = 1e-6` m this drops sub-0.1-mm clusters, far below any real building
/// feature.
const CLUSTER_EPS_FACTOR: f64 = 64.0;

/// `true` if the contour is a degenerate sliver to be regularized away.
///
/// Two complementary criteria, both rooted in the single tolerance:
/// * **Thin sliver**: mean thickness (`2·area / perimeter`) at or below `eps` —
///   a long, hair-thin band from a near-tangent overlap.
/// * **Tiny cluster**: the whole contour fits within `CLUSTER_EPS_FACTOR · eps`
///   — a micro-triangle left by non-transitive snapping.
fn is_sliver(c: &Contour, tol: &Tol) -> bool {
    let area = c.signed_area().abs();
    if area <= 0.0 {
        return true;
    }
    let verts = c.vertices();
    let n = verts.len();

    // Bounding-box extent of the contour.
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut perim = 0.0_f64;
    for i in 0..n {
        let p = verts[i];
        min_x = min_x.min(p.x);
        min_y = min_y.min(p.y);
        max_x = max_x.max(p.x);
        max_y = max_y.max(p.y);
        perim += p.dist(verts[(i + 1) % n]);
    }
    let extent = (max_x - min_x).max(max_y - min_y);
    if extent <= CLUSTER_EPS_FACTOR * tol.length {
        return true;
    }
    if perim <= 0.0 {
        return true;
    }
    // mean thickness = 2*area/perimeter; thin sliver iff <= eps.
    2.0 * area <= tol.length * perim
}

/// Trace one closed loop from `start`, consuming edges from `adj`.
fn trace_one_loop(
    store: &VertexStore,
    adj: &mut HashMap<VertexId, Vec<VertexId>>,
    start: VertexId,
) -> Option<Vec<VertexId>> {
    let mut loop_ids: Vec<VertexId> = Vec::new();
    let first = pop_any(adj, start)?;
    let mut prev = start;
    let mut cur = first;
    loop_ids.push(start);

    let cap = adj.values().map(|v| v.len()).sum::<usize>() + 4;
    let mut steps = 0usize;
    while cur != start {
        loop_ids.push(cur);
        let next = pick_next(store, adj, prev, cur)?;
        prev = cur;
        cur = next;
        steps += 1;
        if steps > cap {
            return None;
        }
    }
    Some(loop_ids)
}

/// Remove and return any outgoing neighbour of `v`.
fn pop_any(adj: &mut HashMap<VertexId, Vec<VertexId>>, v: VertexId) -> Option<VertexId> {
    adj.get_mut(&v)?.pop()
}

/// From `cur` (arrived from `prev`), leave along the edge making the tightest
/// left turn from the reverse of the incoming direction. Consume and return it.
fn pick_next(
    store: &VertexStore,
    adj: &mut HashMap<VertexId, Vec<VertexId>>,
    prev: VertexId,
    cur: VertexId,
) -> Option<VertexId> {
    let outs = adj.get(&cur)?;
    if outs.is_empty() {
        return None;
    }
    let pc = store.point(cur);
    // Reference: direction pointing back to prev.
    let refd = pc.to(store.point(prev));

    let mut best_idx = 0usize;
    let mut best_angle = f64::INFINITY;
    for (i, &nbr) in outs.iter().enumerate() {
        let out_dir = pc.to(store.point(nbr));
        let ang = left_turn_angle(refd, out_dir);
        if ang < best_angle {
            best_angle = ang;
            best_idx = i;
        }
    }
    adj.get_mut(&cur)?.swap_remove(best_idx).into()
}

/// CCW angle in `[0, 2π)` from `refd` to `out_dir`. The smallest such angle is
/// the tightest left turn, which hugs the boundary so outer loops trace CCW.
fn left_turn_angle(refd: Vec2, out_dir: Vec2) -> f64 {
    let cross = refd.cross(out_dir);
    let dot = refd.dot(out_dir);
    let ang = cross.atan2(dot);
    if ang <= 0.0 {
        ang + std::f64::consts::TAU
    } else {
        ang
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tol;

    fn store_with(pts: &[Point2]) -> (VertexStore, Vec<VertexId>) {
        let mut s = VertexStore::new(Tol::default());
        let ids = pts.iter().map(|&p| s.insert(p)).collect();
        (s, ids)
    }

    #[test]
    fn single_ccw_face_round_trips() {
        let (store, ids) = store_with(&[
            Point2::new(0.0_f64, 0.0_f64),
            Point2::new(1.0_f64, 0.0_f64),
            Point2::new(1.0_f64, 1.0_f64),
            Point2::new(0.0_f64, 1.0_f64),
        ]);
        let region = reconstruct(&store, &[ids], &Tol::default());
        assert_eq!(region.contours.len(), 1);
        assert!((region.area() - 1.0_f64).abs() <= 1e-9_f64);
        assert!(region.signed_area() > 0.0_f64);
    }

    #[test]
    fn two_adjacent_faces_merge_dropping_shared_edge() {
        let mut s = VertexStore::new(Tol::default());
        let p = |x: f64, y: f64| Point2::new(x, y);
        let f1: Vec<VertexId> = [p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)]
            .iter()
            .map(|&q| s.insert(q))
            .collect();
        let f2: Vec<VertexId> = [p(1.0, 0.0), p(2.0, 0.0), p(2.0, 1.0), p(1.0, 1.0)]
            .iter()
            .map(|&q| s.insert(q))
            .collect();
        let region = reconstruct(&s, &[f1, f2], &Tol::default());
        assert_eq!(region.contours.len(), 1);
        assert!((region.area() - 2.0_f64).abs() <= 1e-9_f64);
    }
}
