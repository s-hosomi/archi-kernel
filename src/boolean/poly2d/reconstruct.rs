//! Reconstruct a [`Region`] from the set of selected arrangement faces.
//!
//! After classification we have a set of bounded faces (each a directed loop of
//! [`LoopEdge`]s) that should be kept. Their union is the result.
//! Reconstruction:
//!
//! 1. **Cancel internal edges.** Each selected face contributes its boundary as
//!    directed edges. An edge shared by two selected faces appears once in each
//!    direction (`u→v` and `v→u`) and cancels — it is interior to the union, not
//!    part of the result boundary. Arc edges cancel against the reverse arc of
//!    the *same circle*.
//! 2. **Trace boundary loops.** The surviving directed edges form closed loops.
//!    At each vertex we leave along the edge that makes the tightest left turn
//!    from the reverse of the incoming **tangent** direction, which keeps the
//!    trace hugging the boundary so outer loops come out CCW and holes CW.
//! 3. **Normalize orientation.** Outers are emitted CCW, holes CW; zero-area
//!    slivers are dropped. Arc edges are emitted as [`Edge2::Arc`], so the result
//!    keeps its circular boundary (no polyline approximation).
//!
//! Because step 1 cancels shared edges by **exact vertex id + circle key** (the
//! snap already merged coincident vertices), adjacency needs no tolerance here.

use std::collections::HashMap;

use crate::boolean::poly2d::arrangement::FaceLoop;
use crate::boolean::poly2d::geom::{Arc, Edge2, Point2, Vec2};
use crate::boolean::poly2d::region::{Contour, Region};
use crate::boolean::poly2d::snap::{VertexId, VertexStore};
use crate::boolean::support::quantize;
use crate::tolerance::Tol;

/// A directed boundary edge of the result, carrying enough geometry to emit an
/// arc and to order around a vertex.
#[derive(Debug, Clone, Copy)]
struct DEdge {
    a: VertexId,
    b: VertexId,
    edge: Edge2,
}

/// Rebuild a region from the selected face loops.
pub fn reconstruct(store: &VertexStore, selected: &[&FaceLoop], tol: &Tol) -> Region {
    // ── 1. cancel internal edges ───────────────────────────────────────────
    // Key each directed edge by (a, b, geom-discriminator). A forward edge
    // cancels its exact reverse on the same curve.
    let mut count: HashMap<(VertexId, VertexId, GeomKey), (i32, Edge2)> = HashMap::new();
    for face in selected {
        for le in &face.edges {
            if le.a == le.b {
                continue;
            }
            let edge = le.to_edge2(store);
            let gk = geom_key(&edge);
            count.entry((le.a, le.b, gk)).or_insert((0, edge)).0 += 1;
        }
    }

    // Surviving directed edges: forward count minus reverse count (same curve).
    let mut adj: HashMap<VertexId, Vec<DEdge>> = HashMap::new();
    let mut keys: Vec<(VertexId, VertexId, GeomKey)> = count.keys().copied().collect();
    keys.sort_by_key(|x| (x.0, x.1, x.2));
    for k in keys {
        let (a, b, gk) = k;
        let fwd = count.get(&(a, b, gk)).map(|e| e.0).unwrap_or(0);
        let rev = count.get(&(b, a, gk.reversed())).map(|e| e.0).unwrap_or(0);
        let net = fwd - rev;
        if net > 0 {
            let edge = count.get(&(a, b, gk)).map(|e| e.1).unwrap();
            for _ in 0..net {
                adj.entry(a).or_default().push(DEdge { a, b, edge });
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
            if let Some(loop_edges) = trace_one_loop(store, &mut adj, start) {
                if let Some(contour) = build_contour(&loop_edges) {
                    contours.push(contour);
                }
            } else {
                break;
            }
        }
    }

    contours.retain(|c| !is_sliver(c, tol));
    Region::new(contours)
}

/// Geometry discriminator for cancellation and ordering.
///
/// An arc carries its circle (centre + radius) **and** a quantised midpoint of
/// the physical arc, so the two semicircles of one circle between the same seam
/// vertices stay distinct (they bulge to opposite sides → opposite midpoints) and
/// do not wrongly cancel, while a forward arc `a→b` and its reverse `b→a` share
/// the same midpoint and cancel correctly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum GeomKey {
    Seg,
    Arc {
        cx: i64,
        cy: i64,
        r: i64,
        mx: i64,
        my: i64,
    },
}

impl GeomKey {
    /// The key of the reverse edge: identical (same curve + same midpoint), so
    /// cancellation of `a→b` against `b→a` matches on the same `GeomKey`.
    fn reversed(self) -> Self {
        self
    }
}

fn geom_key(edge: &Edge2) -> GeomKey {
    match edge {
        Edge2::Seg { .. } => GeomKey::Seg,
        Edge2::Arc(a) => {
            let mid = a.mid_point();
            GeomKey::Arc {
                cx: quantize(a.center.x),
                cy: quantize(a.center.y),
                r: quantize(a.radius),
                mx: quantize(mid.x),
                my: quantize(mid.y),
            }
        }
    }
}

/// Multiple of `eps` below which a whole contour is treated as a snap-cluster
/// artifact.
const CLUSTER_EPS_FACTOR: f64 = 64.0;

/// `true` if the contour is a degenerate sliver to be regularized away.
fn is_sliver(c: &Contour, tol: &Tol) -> bool {
    let area = c.signed_area().abs();
    if area <= 0.0 {
        return true;
    }
    // An arc contour (e.g. a small circular hole or island) is real even though
    // its bounding box can be small; only the straight-edge sliver / cluster
    // heuristics apply to polygonal contours.
    if c.has_arc() {
        return area <= tol.length * tol.length;
    }
    let verts = c.vertices();
    let n = verts.len();
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
    2.0 * area <= tol.length * perim
}

/// Trace one closed loop from `start`, consuming directed edges from `adj`.
fn trace_one_loop(
    store: &VertexStore,
    adj: &mut HashMap<VertexId, Vec<DEdge>>,
    start: VertexId,
) -> Option<Vec<DEdge>> {
    let mut loop_edges: Vec<DEdge> = Vec::new();
    let first = pop_any(adj, start)?;
    let mut prev_in = first; // the edge we arrived on
    loop_edges.push(first);
    let mut cur = first.b;

    let cap = adj.values().map(|v| v.len()).sum::<usize>() + 4;
    let mut steps = 0usize;
    while cur != start {
        let next = pick_next(store, adj, prev_in, cur)?;
        loop_edges.push(next);
        prev_in = next;
        cur = next.b;
        steps += 1;
        if steps > cap {
            return None;
        }
    }
    Some(loop_edges)
}

fn pop_any(adj: &mut HashMap<VertexId, Vec<DEdge>>, v: VertexId) -> Option<DEdge> {
    adj.get_mut(&v)?.pop()
}

/// From `cur` (arrived via `prev_in`), leave along the edge making the tightest
/// left turn from the reverse of the incoming tangent direction. Consume it.
fn pick_next(
    store: &VertexStore,
    adj: &mut HashMap<VertexId, Vec<DEdge>>,
    prev_in: DEdge,
    cur: VertexId,
) -> Option<DEdge> {
    let outs = adj.get(&cur)?;
    if outs.is_empty() {
        return None;
    }
    // Incoming tangent (direction of travel arriving at `cur`); we turn from its
    // reverse.
    let in_tan = edge_tangent_at_end(store, &prev_in);
    let refd = Vec2::new(-in_tan.x, -in_tan.y);

    let mut best_idx = 0usize;
    let mut best_angle = f64::INFINITY;
    for (i, de) in outs.iter().enumerate() {
        let out_tan = edge_tangent_at_start(store, de);
        let ang = left_turn_angle(refd, out_tan);
        if ang < best_angle {
            best_angle = ang;
            best_idx = i;
        }
    }
    Some(adj.get_mut(&cur)?.swap_remove(best_idx))
}

/// The unit tangent of an edge at its **end** (arrival direction).
fn edge_tangent_at_end(store: &VertexStore, de: &DEdge) -> Vec2 {
    match de.edge {
        Edge2::Seg { .. } => {
            let d = store.point(de.a).to(store.point(de.b));
            unit(d)
        }
        Edge2::Arc(arc) => arc_tangent(&arc, store.point(de.b), arc.sweep >= 0.0),
    }
}

/// The unit tangent of an edge at its **start** (departure direction).
fn edge_tangent_at_start(store: &VertexStore, de: &DEdge) -> Vec2 {
    match de.edge {
        Edge2::Seg { .. } => {
            let d = store.point(de.a).to(store.point(de.b));
            unit(d)
        }
        Edge2::Arc(arc) => arc_tangent(&arc, store.point(de.a), arc.sweep >= 0.0),
    }
}

/// Tangent of an arc's circle at point `p`, in the traversal direction.
fn arc_tangent(arc: &Arc, p: Point2, ccw: bool) -> Vec2 {
    let rad = arc.center.to(p);
    let t = if ccw {
        Vec2::new(-rad.y, rad.x)
    } else {
        Vec2::new(rad.y, -rad.x)
    };
    unit(t)
}

fn unit(v: Vec2) -> Vec2 {
    let l = v.len();
    if l > 0.0 {
        Vec2::new(v.x / l, v.y / l)
    } else {
        Vec2::new(0.0, 0.0)
    }
}

/// CCW angle in `[0, 2π)` from `refd` to `out_dir`.
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

/// Build a [`Contour`] from a traced directed-edge loop, preserving arcs.
fn build_contour(loop_edges: &[DEdge]) -> Option<Contour> {
    if loop_edges.len() < 2 {
        return None;
    }
    let edges: Vec<Edge2> = loop_edges.iter().map(|de| de.edge).collect();
    Some(Contour::new(edges))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean::poly2d::arrangement::Arrangement;
    use crate::boolean::poly2d::Region;
    use crate::tolerance::Tol;

    #[test]
    fn single_ccw_face_round_trips() {
        let tol = Tol::default();
        let a = Region::from_points(&[
            Point2::new(0.0_f64, 0.0_f64),
            Point2::new(1.0_f64, 0.0_f64),
            Point2::new(1.0_f64, 1.0_f64),
            Point2::new(0.0_f64, 1.0_f64),
        ]);
        let arr = Arrangement::build(&a, &Region::empty(), &tol).unwrap();
        let faces: Vec<&FaceLoop> = arr.faces.iter().filter(|f| f.signed_area() > 0.0).collect();
        let region = reconstruct(arr.store(), &faces, &tol);
        assert_eq!(region.contours.len(), 1);
        assert!((region.area() - 1.0_f64).abs() <= 1e-9_f64);
        assert!(region.signed_area() > 0.0_f64);
    }
}
