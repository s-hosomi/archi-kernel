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
//! * a set of **atomic edges** (the minimal segments / arcs of the subdivision),
//!   each carrying the two cells it borders (left / right) and its geometry.
//!
//! Walls are generated per atomic edge per band and interface caps per atomic
//! cell per level. A straight atomic edge becomes a planar wall; an **arc**
//! atomic edge becomes a **cylinder wall** (Phase 3c). Because both walls and
//! caps read the *same* vertices and circle, a wall's vertical edge and an
//! interface's boundary edge meet at identical coordinates and the cylinder
//! surface is shared, so the 3-D vertex / curve / surface interning in
//! [`super::build`] pairs every sibling exactly — watertightness by construction.
//!
//! The snapping, edge intersection and exact `orient2d` are reused verbatim from
//! the proven 2-D engine ([`crate::boolean::poly2d`]).

use std::collections::{BTreeMap, HashMap};

use crate::boolean::poly2d::geom::{orient2d, Arc, Edge2, Orient, Point2};
use crate::boolean::poly2d::intersect::intersect;
use crate::boolean::poly2d::snap::{VertexId, VertexStore};
use crate::boolean::poly2d::{Poly2Error, Region};
use crate::tolerance::Tol;

use super::error::PrismError;

/// Geometry kind of an arrangement edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum EdgeGeom {
    /// A straight segment.
    Seg,
    /// An arc on a circle (centre + radius), traversed CCW when `ccw`.
    Arc {
        center: Point2,
        radius: f64,
        ccw: bool,
    },
}

/// A directed input boundary edge after snapping its endpoints, tagged with the
/// index of the operand region it came from and its geometry.
#[derive(Debug, Clone, Copy)]
struct InputEdge {
    a: VertexId,
    b: VertexId,
    operand: usize,
    geom: EdgeGeom,
}

/// One directed boundary edge of a cell ring, with geometry.
#[derive(Debug, Clone, Copy)]
pub(super) struct RingEdge {
    pub a: VertexId,
    pub b: VertexId,
    pub geom: EdgeGeom,
}

/// An atomic cell of the arrangement: a traced bounded face (with any inner hole
/// rings) and its residency.
#[derive(Debug, Clone)]
pub(super) struct Cell {
    /// Directed edges of the cell's outer boundary loop, in CCW order.
    pub outer: Vec<RingEdge>,
    /// Inner hole rings (CW), each a loop of directed edges.
    pub inner_rings: Vec<Vec<RingEdge>>,
    /// `inside[i]` is `true` if an interior point of the cell (in the annulus,
    /// not in any hole) lies inside the `i`-th operand region.
    pub inside: Vec<bool>,
}

/// An atomic undirected edge of the arrangement, with its two bordering cells.
///
/// `left` is the cell on the left of the directed edge `a → b`; `right` is the
/// cell on its right. Either may be `None` for the unbounded exterior. `geom` is
/// the geometry as traversed `a → b`.
#[derive(Debug, Clone, Copy)]
pub(super) struct ArrEdge {
    pub a: VertexId,
    pub b: VertexId,
    pub left: Option<usize>,
    pub right: Option<usize>,
    pub geom: EdgeGeom,
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
        let mut inputs: Vec<InputEdge> = Vec::new();
        for (i, r) in regions.iter().enumerate() {
            ingest(r, i, &mut store, &mut inputs)?;
        }

        // Split every input edge at every mutual / self intersection.
        let n = inputs.len();
        let mut split_points: Vec<Vec<VertexId>> = vec![Vec::new(); n];
        for i in 0..n {
            let ei = input_edge2(&store, &inputs[i]);
            for j in (i + 1)..n {
                let ej = input_edge2(&store, &inputs[j]);
                let cr = intersect(&ei, &ej, tol).map_err(map_poly2)?;
                for p in cr.points {
                    let v = store.insert(p);
                    push_split(&mut split_points[i], v, &inputs[i]);
                    push_split(&mut split_points[j], v, &inputs[j]);
                }
            }
        }

        // Dedup the split fragments into undirected arrangement edges, keyed on
        // the canonical vertex-id pair **plus a geometry discriminator** (so two
        // arcs of the same circle on opposite sides do not collide). A `BTreeMap`
        // keeps the order deterministic.
        let mut edge_map: BTreeMap<(VertexId, VertexId, GeomKey), EdgeGeom> = BTreeMap::new();
        for (i, edge) in inputs.iter().enumerate() {
            let chain = ordered_chain(&store, edge, &split_points[i], tol);
            for w in chain.windows(2) {
                let (u, v) = (w[0], w[1]);
                if u == v {
                    continue;
                }
                let (key_uv, dir) = if u <= v {
                    ((u, v), true)
                } else {
                    ((v, u), false)
                };
                let canon_geom = canonical_geom(edge.geom, dir);
                let pa = store.point(key_uv.0);
                let pb = store.point(key_uv.1);
                let gk = geom_key_for(pa, pb, &canon_geom);
                edge_map
                    .entry((key_uv.0, key_uv.1, gk))
                    .or_insert(canon_geom);
            }
        }
        let undirected: Vec<(VertexId, VertexId, EdgeGeom)> = edge_map
            .into_iter()
            .map(|((a, b, _), g)| (a, b, g))
            .collect();

        // Build the DCEL and trace its cells, recording cell adjacency per edge.
        let (cells_raw, edges) = trace(&store, &undirected)?;

        // Tag each cell with its residency in every operand by exact winding at an
        // interior point.
        let n_ops = regions.len();
        let mut cells = Vec::with_capacity(cells_raw.len());
        for raw in cells_raw {
            let sample = face_sample_point(&store, &raw.outer);
            let inside: Vec<bool> = (0..n_ops)
                .map(|op| winding(&store, &inputs, op, sample) != 0)
                .collect();
            cells.push(Cell {
                outer: raw.outer,
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

/// Map a 2-D engine error into a prismatic error (arc-degeneracy → Phase 3c).
fn map_poly2(e: Poly2Error) -> PrismError {
    PrismError::from(e)
}

/// Geometry discriminator for the dedup key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
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
                cx: (center.x * 1e9).round() as i64,
                cy: (center.y * 1e9).round() as i64,
                r: (radius * 1e9).round() as i64,
                mx: (mid.x * 1e9).round() as i64,
                my: (mid.y * 1e9).round() as i64,
            }
        }
    }
}

/// The geometry as traversed in canonical (a ≤ b) order: flip arc direction when
/// the edge was reversed to canonical order.
fn canonical_geom(geom: EdgeGeom, forward: bool) -> EdgeGeom {
    match geom {
        EdgeGeom::Seg => EdgeGeom::Seg,
        EdgeGeom::Arc {
            center,
            radius,
            ccw,
        } => EdgeGeom::Arc {
            center,
            radius,
            ccw: if forward { ccw } else { !ccw },
        },
    }
}

/// Signed sweep from `pa` to `pb` about `center` in the `ccw` direction.
fn directed_sweep(center: Point2, pa: Point2, pb: Point2, ccw: bool) -> f64 {
    let aa = (pa.y - center.y).atan2(pa.x - center.x);
    let ab = (pb.y - center.y).atan2(pb.x - center.x);
    if ccw {
        (ab - aa).rem_euclid(std::f64::consts::TAU)
    } else {
        -((aa - ab).rem_euclid(std::f64::consts::TAU))
    }
}

/// Build the [`Edge2`] for an `a → b` edge with the given geometry.
fn edge2_of(store: &VertexStore, a: VertexId, b: VertexId, geom: EdgeGeom) -> Edge2 {
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

/// Snap every contour edge endpoint and emit directed input edges.
fn ingest(
    r: &Region,
    operand: usize,
    store: &mut VertexStore,
    inputs: &mut Vec<InputEdge>,
) -> Result<(), PrismError> {
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

#[inline]
fn input_edge2(store: &VertexStore, s: &InputEdge) -> Edge2 {
    edge2_of(store, s.a, s.b, s.geom)
}

fn push_split(splits: &mut Vec<VertexId>, v: VertexId, s: &InputEdge) {
    if v != s.a && v != s.b && !splits.contains(&v) {
        splits.push(v);
    }
}

/// Order an edge's interior split vertices along the edge (chord param for a
/// segment, sweep angle for an arc).
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
            let start_angle = (a.y - center.y).atan2(a.x - center.x);
            let sweep = directed_sweep(center, a, b, ccw);
            let arc = Arc::new(center, radius, start_angle, sweep);
            let span = sweep.abs().max(f64::MIN_POSITIVE);
            let mut mids: Vec<(f64, VertexId)> = splits
                .iter()
                .filter_map(|&v| {
                    let p = store.point(v);
                    arc.angle_of_point(p, tol)
                        .map(|theta| ((theta - start_angle).abs() / span, v))
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

/// A directed half-edge in the trace DCEL.
struct Half {
    origin: VertexId,
    dest: VertexId,
    twin: usize,
    next: usize,
    geom: EdgeGeom,
    cell: Option<usize>,
}

/// A traced cell: its CCW outer ring plus any CW inner hole rings.
struct RawCell {
    outer: Vec<RingEdge>,
    inners: Vec<Vec<RingEdge>>,
}

/// Build the DCEL from undirected edges, trace cells, and return the cells'
/// outer/inner loops plus the edge adjacency.
fn trace(
    store: &VertexStore,
    undirected: &[(VertexId, VertexId, EdgeGeom)],
) -> Result<(Vec<RawCell>, Vec<ArrEdge>), PrismError> {
    let mut halfs: Vec<Half> = Vec::with_capacity(undirected.len() * 2);
    let mut outgoing: HashMap<VertexId, Vec<usize>> = HashMap::new();
    for &(a, b, geom) in undirected.iter() {
        let h0 = halfs.len();
        let h1 = h0 + 1;
        let rev = reverse_geom(geom);
        halfs.push(Half {
            origin: a,
            dest: b,
            twin: h1,
            next: usize::MAX,
            geom,
            cell: None,
        });
        halfs.push(Half {
            origin: b,
            dest: a,
            twin: h0,
            next: usize::MAX,
            geom: rev,
            cell: None,
        });
        outgoing.entry(a).or_default().push(h0);
        outgoing.entry(b).or_default().push(h1);
    }

    for (&v, outs) in &mut outgoing {
        let vp = store.point(v);
        outs.sort_by(|&x, &y| {
            let ax = leave_angle(store, &halfs[x], vp);
            let ay = leave_angle(store, &halfs[y], vp);
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

    struct Loop {
        ring: Vec<RingEdge>,
        pts: Vec<Point2>,
        chain: Vec<usize>,
        area: f64,
    }
    let mut loops: Vec<Loop> = Vec::new();
    for start in 0..n {
        if halfs[start].cell.is_some() || halfs[start].next == usize::MAX {
            continue;
        }
        let mut ring: Vec<RingEdge> = Vec::new();
        let mut pts: Vec<Point2> = Vec::new();
        let mut chain: Vec<usize> = Vec::new();
        let mut cur = start;
        let cap = n + 1;
        let mut steps = 0;
        loop {
            halfs[cur].cell = Some(usize::MAX);
            let h = &halfs[cur];
            ring.push(RingEdge {
                a: h.origin,
                b: h.dest,
                geom: h.geom,
            });
            pts.push(store.point(h.origin));
            chain.push(cur);
            cur = halfs[cur].next;
            steps += 1;
            if cur == start || cur == usize::MAX || steps > cap {
                break;
            }
        }
        let area = ring_signed_area(store, &ring, &pts);
        loops.push(Loop {
            ring,
            pts,
            chain,
            area,
        });
    }

    let mut half_cell: Vec<Option<usize>> = vec![None; n];

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

    let ring_key = |ring: &[RingEdge]| {
        let mut k: Vec<usize> = ring.iter().flat_map(|e| [e.a.0, e.b.0]).collect();
        k.sort_unstable();
        k
    };
    let pos_keys: Vec<Vec<usize>> = pos_loops
        .iter()
        .map(|&li| ring_key(&loops[li].ring))
        .collect();

    let _ = &pos_keys;
    for lp in &loops {
        if lp.area >= 0.0 {
            continue;
        }
        // Each negative (CW) loop bounds an enclosed void. It is an inner hole of
        // the **smallest positive cell whose outer ring contains the void**: that
        // cell's material surrounds the void. The void itself may be a single cell
        // (a simple round hole) or a *composite* region (several cells of a
        // merged/overlapping void — a fused circular-sleeve pair); either way the
        // attachment is by geometric containment, not by matching one cell's
        // boundary (which fails for the composite case and dropped the void). The
        // probe point is taken strictly inside the void.
        let void_ci = pos_keys.iter().position(|k| *k == ring_key(&lp.ring));
        // Probe the **material** side (just outside the CW void loop, to the left
        // of its directed edges): that point lies in the surrounding annulus, not
        // in any enclosed void cell, so the smallest positive cell containing it is
        // the correct parent even when the void is composite (overlapping sleeves).
        let probe = material_probe_point(&lp.pts);
        let mut parent: Option<usize> = None;
        let mut parent_area = f64::INFINITY;
        for (ci, &li) in pos_loops.iter().enumerate() {
            if Some(ci) == void_ci {
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

    let mut edges: Vec<ArrEdge> = Vec::with_capacity(undirected.len());
    for (ei, &(a, b, geom)) in undirected.iter().enumerate() {
        let h_ab = 2 * ei;
        let h_ba = 2 * ei + 1;
        debug_assert_eq!(halfs[h_ab].origin, a);
        debug_assert_eq!(halfs[h_ab].dest, b);
        edges.push(ArrEdge {
            a,
            b,
            left: half_cell[h_ab],
            right: half_cell[h_ba],
            geom,
        });
    }

    Ok((cells, edges))
}

/// Reverse an edge's geometry (flip arc direction).
fn reverse_geom(geom: EdgeGeom) -> EdgeGeom {
    match geom {
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
    }
}

/// A point just **outside** a CW loop (on the material side), used to find the
/// enclosing parent cell. For a CW loop the interior (void) is on the right of
/// each directed edge, so the **left** normal of an edge midpoint points outward
/// into the surrounding material. The step is shrunk until the point is verified
/// *outside* the loop, so a thin or composite (peanut) void is sampled in the
/// material, never in an enclosed void cell.
fn material_probe_point(pts: &[Point2]) -> Point2 {
    let n = pts.len();
    if n < 2 {
        return pts.first().copied().unwrap_or(Point2::new(0.0, 0.0));
    }
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
    // Left normal of a→b is (−dy, dx)/len; for a CW loop this points outward
    // (into the surrounding material).
    let nx = -d.y / len;
    let ny = d.x / len;
    let mx = (a.x + b.x) * 0.5;
    let my = (a.y + b.y) * 0.5;
    let mut step = 1e-4_f64 * len;
    for _ in 0..60 {
        let p = Point2::new(mx + nx * step, my + ny * step);
        if !point_in_ring(pts, p) {
            return p;
        }
        step *= 0.5;
    }
    Point2::new(mx + nx * 1e-7, my + ny * 1e-7)
}

/// Angle at which a half-edge leaves its origin (tangent direction).
fn leave_angle(store: &VertexStore, h: &Half, vp: Point2) -> f64 {
    let dir = match h.geom {
        EdgeGeom::Seg => vp.to(store.point(h.dest)),
        EdgeGeom::Arc { center, ccw, .. } => {
            let rad = center.to(vp);
            if ccw {
                crate::boolean::poly2d::Vec2::new(-rad.y, rad.x)
            } else {
                crate::boolean::poly2d::Vec2::new(rad.y, -rad.x)
            }
        }
    };
    dir.y.atan2(dir.x)
}

/// Arc-aware signed area of a ring.
fn ring_signed_area(_store: &VertexStore, ring: &[RingEdge], pts: &[Point2]) -> f64 {
    let n = pts.len();
    if n < 2 {
        return 0.0;
    }
    let mut acc = 0.0_f64;
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        acc += a.x * b.y - b.x * a.y;
    }
    let mut area = 0.5 * acc;
    for (i, e) in ring.iter().enumerate() {
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

/// `true` if `p` is strictly inside the polygon ring `pts` (chord polygon test).
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

/// A point strictly inside the cell bounded by `ring`, robust against arc-edge
/// tangencies (sampled at a generic fraction of the longest edge, offset along
/// its left normal until the cell's self-winding confirms it is inside).
fn face_sample_point(store: &VertexStore, ring: &[RingEdge]) -> Point2 {
    let n = ring.len();
    if n == 0 {
        return Point2::new(0.0, 0.0);
    }
    let mut best_i = 0usize;
    let mut best_len = -1.0_f64;
    for (i, e) in ring.iter().enumerate() {
        let l = store.point(e.a).dist(store.point(e.b));
        if l > best_len {
            best_len = l;
            best_i = i;
        }
    }
    let e = ring[best_i];
    let ed = edge2_of(store, e.a, e.b, e.geom);
    let base = ed.point_at(0.37);
    let tan = tangent_at(&ed, 0.37);
    if tan.len_sq() <= 0.0 {
        return centroid(store, ring);
    }
    let nx = -tan.y;
    let ny = tan.x;
    let mut step = 1e-6_f64;
    for _ in 0..40 {
        let p = Point2::new(base.x + nx * step, base.y + ny * step);
        if ring_self_winding(store, ring, p) == 1 {
            return p;
        }
        step *= 0.5;
    }
    centroid(store, ring)
}

/// Tangent of an edge at fraction `t`.
fn tangent_at(edge: &Edge2, t: f64) -> crate::boolean::poly2d::Vec2 {
    use crate::boolean::poly2d::Vec2;
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

/// Winding number of `p` with respect to one ring's own edges.
fn ring_self_winding(store: &VertexStore, ring: &[RingEdge], p: Point2) -> i32 {
    let mut w = 0;
    for e in ring {
        let a = store.point(e.a);
        let b = store.point(e.b);
        match e.geom {
            EdgeGeom::Seg => w += ray_seg(p, a, b),
            EdgeGeom::Arc {
                center,
                radius,
                ccw,
            } => w += ray_arc(p, a, b, center, radius, ccw),
        }
    }
    w
}

fn centroid(store: &VertexStore, ring: &[RingEdge]) -> Point2 {
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    for e in ring {
        let p = store.point(e.a);
        cx += p.x;
        cy += p.y;
    }
    let k = ring.len().max(1) as f64;
    Point2::new(cx / k, cy / k)
}

/// Winding number of `p` w.r.t. the given operand's original directed edges,
/// counting both straight and arc edges.
fn winding(store: &VertexStore, inputs: &[InputEdge], operand: usize, p: Point2) -> i32 {
    let mut w = 0_i32;
    for s in inputs.iter().filter(|s| s.operand == operand) {
        let a = store.point(s.a);
        let b = store.point(s.b);
        match s.geom {
            EdgeGeom::Seg => w += ray_seg(p, a, b),
            EdgeGeom::Arc {
                center,
                radius,
                ccw,
            } => w += ray_arc(p, a, b, center, radius, ccw),
        }
    }
    w
}

/// Segment ray-crossing contribution.
fn ray_seg(p: Point2, a: Point2, b: Point2) -> i32 {
    if a.y <= p.y {
        if b.y > p.y && orient2d(a, b, p) == Orient::Left {
            return 1;
        }
    } else if b.y <= p.y && orient2d(a, b, p) == Orient::Right {
        return -1;
    }
    0
}

/// Arc ray-crossing contribution (see the poly2d arrangement for the rationale).
fn ray_arc(p: Point2, a: Point2, b: Point2, center: Point2, radius: f64, ccw: bool) -> i32 {
    let dy = p.y - center.y;
    let graze = 1e-9_f64;
    if dy.abs() >= radius - graze {
        return 0;
    }
    let dx = (radius * radius - dy * dy).max(0.0).sqrt();
    let start_angle = (a.y - center.y).atan2(a.x - center.x);
    let sweep = directed_sweep(center, a, b, ccw);
    let arc = Arc::new(center, radius, start_angle, sweep);
    let tol = Tol::default();
    let mut acc = 0_i32;
    for &xx in &[center.x + dx, center.x - dx] {
        let pt = Point2::new(xx, p.y);
        if xx <= p.x {
            continue;
        }
        if arc.angle_of_point(pt, &tol).is_none() {
            continue;
        }
        if pt.coincident(arc.end(), &tol) {
            continue;
        }
        let theta = (pt.y - center.y).atan2(pt.x - center.x);
        let tan_y = if ccw { theta.cos() } else { -theta.cos() };
        if tan_y > 0.0 {
            acc += 1;
        } else if tan_y < 0.0 {
            acc -= 1;
        }
    }
    acc
}
