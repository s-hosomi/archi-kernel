//! Build the 3-D result B-rep from the shared arrangement and the axial bands.
//!
//! Given the global 2-D [`Arrangement`] (atomic cells + edges, each cell tagged
//! with `(in_a, in_b)`) and the two operands' axial intervals, this module
//! constructs the prismatic boolean result as a watertight [`Brep`]
//! (`DESIGN.md` §4.2).
//!
//! # Bands
//!
//! The axial breakpoints are the four interval endpoints `{t0a, t1a, t0b, t1b}`
//! merged within tolerance. Consecutive distinct breakpoints define **bands**
//! `[z_k, z_{k+1}]`. Within a band each operand is either wholly present (the
//! band lies inside its interval) or wholly absent, so a cell's residency in the
//! band is `op.keep(in_a ∧ a_present, in_b ∧ b_present)` — the keep table of
//! `DESIGN.md` §4.4 evaluated once per `(cell, band)`.
//!
//! # Faces (watertight by shared interning)
//!
//! The solid is the set of resident `cell × band` voxels. Its boundary has two
//! kinds of face, generated from the **one** shared segmentation so coincident
//! edges meet exactly:
//!
//! * **Walls** (vertical, parallel to `d`): for each arrangement edge `g` and
//!   band `k`, if the two cells across `g` differ in residency, emit the quad
//!   `g × [z_k, z_{k+1}]`, oriented outward (toward the void side).
//! * **Interfaces** (horizontal, ⟂ `d`): for each cell `c` and level `z_k`, if
//!   residency differs between band `k−1` (or empty below the bottom) and band
//!   `k` (or empty above the top), emit the cell polygon at `z_k`, facing `+d`
//!   when material is below and `−d` when material is above. The bottom-most and
//!   top-most caps fall out of this rule (`DESIGN.md` §4.2 items 4–5).
//!
//! Every 3-D vertex is interned by quantised coordinate and every straight edge
//! by its endpoint-key pair, exactly as the extruder and the cut do, so two
//! faces sharing an edge reach the same [`CurveId`](crate::geom::CurveId) and the
//! sibling pairing (same curve, reversed boundary) succeeds — the result is
//! watertight even though cells are over-split.
//!
//! # Multiple solids
//!
//! A boolean can split a member into disconnected pieces (a wall cut full-height
//! by an opening becomes two solids). After all faces are built they are grouped
//! into connected components by **shared edge** (sibling) adjacency, and each
//! component becomes its own [`Solid`] (`DESIGN.md` §4.2 item 6).

use std::collections::HashMap;

use crate::boolean::poly2d::snap::VertexId;
use crate::boolean::poly2d::Op;
use crate::boolean::support::{key, quantize, uf_find, uf_union, CoordKey};
use crate::brep::Brep;
use crate::geom::{CurveGeom, CurveId, SurfaceGeom, VertexGeom};
use crate::math::{Point3, Unit3, Vec3};
use crate::primitives::{Circle3, Cylinder, Line3, Plane};
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::validate::ValidateLevel;
use crate::topo::{Face, HalfEdge, Loop, Sense, Shell, Solid, Vertex};

use super::arrange::{Arrangement, EdgeGeom, RingEdge};
use super::detect::{Frame, PrismOperand};
use super::error::PrismError;

/// A 2-D directed boundary edge of a cap ring, in the shared frame's plane,
/// carrying the geometry needed to build the matching 3-D edge (a straight line
/// or a circular arc).
#[derive(Debug, Clone, Copy)]
pub(super) enum BoundaryEdge2 {
    /// A straight segment `a → b` (2-D frame coords), with the wedge tags of its
    /// endpoints (`wa` at `a`, `wb` at `b`; see [`Wedges`]). Both default to `0`
    /// except at a non-manifold corner pinch.
    Seg {
        a: [f64; 2],
        b: [f64; 2],
        wa: u32,
        wb: u32,
    },
    /// An arc `a → b` on a circle of `center` / `radius`, swept `ccw` or not.
    Arc {
        a: [f64; 2],
        b: [f64; 2],
        center: [f64; 2],
        radius: f64,
        ccw: bool,
    },
}

/// Convert a cell ring of [`RingEdge`]s into 2-D [`BoundaryEdge2`]s, tagging each
/// straight endpoint with its wedge via `tag` (mapping a 2-D vertex id to its
/// wedge tag on the cap's material side).
fn ring_to_bedges(
    arr: &Arrangement,
    ring: &[RingEdge],
    tag: &impl Fn(VertexId, VertexId, VertexId) -> u32,
) -> Vec<BoundaryEdge2> {
    ring.iter()
        .map(|e| {
            let pa = arr.point(e.a);
            let pb = arr.point(e.b);
            match e.geom {
                EdgeGeom::Seg => BoundaryEdge2::Seg {
                    a: [pa.x, pa.y],
                    b: [pb.x, pb.y],
                    // The corner at each endpoint is bounded by this ring edge.
                    wa: tag(e.a, e.a, e.b),
                    wb: tag(e.b, e.a, e.b),
                },
                EdgeGeom::Arc {
                    center,
                    radius,
                    ccw,
                } => BoundaryEdge2::Arc {
                    a: [pa.x, pa.y],
                    b: [pb.x, pb.y],
                    center: [center.x, center.y],
                    radius,
                    ccw,
                },
            }
        })
        .collect()
}

/// Build the prismatic boolean result for a binary `op` from two operands.
///
/// `budget` bounds the total boundary complexity (`DESIGN.md` §4.5); exceeding
/// it yields [`PrismError::ComplexityLimit`] rather than an unbounded build.
pub(super) fn build(
    frame: &Frame,
    a: &PrismOperand,
    b: &PrismOperand,
    op: Op,
    tol: &Tol,
    budget: usize,
) -> Result<Brep, PrismError> {
    // The keep table folds the two operands' (present ∧ inside) flags.
    build_combined(
        frame,
        &[a.clone(), b.clone()],
        move |f| op.keep(f[0], f[1]),
        tol,
        budget,
    )
}

/// Build a prismatic result over **any number** of operands, combined by
/// `keep`. `keep` receives, per cell per band, the slice of `present ∧ inside`
/// flags (one per operand, in input order) and returns whether the voxel is
/// material. This single path serves the binary ops (a 2-element fold) and
/// opening subtraction (`base ∧ ¬(opening₀ ∨ opening₁ ∨ …)`), so the multi-band
/// arrangement is shared (`DESIGN.md` §4.2, §4.5).
pub(super) fn build_combined(
    frame: &Frame,
    operands: &[PrismOperand],
    keep: impl Fn(&[bool]) -> bool,
    tol: &Tol,
    budget: usize,
) -> Result<Brep, PrismError> {
    let regions: Vec<&_> = operands.iter().map(|o| &o.region).collect();

    // First-line budget guard on the raw input size, *before* the arrangement is
    // built. A pathological operand set (e.g. thousands of openings) can blow up
    // the O(n²) pairwise split inside the arrangement; checking the input edge
    // total up front isolates it without first paying that quadratic cost. The
    // arrangement's own output measure is still checked below.
    let bands = bands(operands, tol);
    let input_edges: usize = regions
        .iter()
        .map(|r| r.contours.iter().map(|c| c.edges.len()).sum::<usize>())
        .sum();
    let input_measure = input_edges
        .saturating_mul(input_edges)
        .saturating_mul(bands.len().max(1));
    if input_measure > budget {
        return Err(PrismError::ComplexityLimit {
            measure: input_measure,
            budget,
        });
    }

    let arr = Arrangement::build(&regions, tol)?;

    let measure = arr.edges.len().saturating_mul(bands.len().max(1))
        + arr.cells.len().saturating_mul(bands.len() + 1);
    if measure > budget {
        return Err(PrismError::ComplexityLimit { measure, budget });
    }
    if bands.is_empty() || arr.cells.is_empty() {
        // Nothing survives (e.g. disjoint intervals, or an empty first operand).
        return Ok(Brep::new());
    }

    let mut builder = Builder::new(frame, *tol);

    // Residency of each cell in each band: feed `keep` the per-operand flags.
    let mut flags = vec![false; operands.len()];
    let resident: Vec<Vec<bool>> = arr
        .cells
        .iter()
        .map(|c| {
            bands
                .iter()
                .map(|bd| {
                    for (i, f) in flags.iter_mut().enumerate() {
                        *f = c.inside[i] && bd.present[i];
                    }
                    keep(&flags)
                })
                .collect()
        })
        .collect();

    // ── connected components over resident voxels (cell × band) ────────────
    // A voxel is `(cell, band)`. Two resident voxels belong to the same solid
    // only when they are genuinely **face-adjacent**: they share an arrangement
    // edge in the same band (a 2-D edge with both cells resident), or they are
    // the same cell in vertically-adjacent resident bands. Two voxels meeting at
    // only a 2-D *vertex* (the checkerboard corner touch) are NOT adjacent, so
    // they land in different components — and are then interned independently,
    // so their shared corner edge does not become a non-manifold 4-way edge.
    let ncells = arr.cells.len();
    let voxel_id = |ci: usize, k: usize| ci * nbands_for(&bands) + k;
    let nvox = ncells * nbands_for(&bands);
    let mut vparent: Vec<usize> = (0..nvox).collect();
    // Horizontal adjacency: cells sharing an arrangement edge, both resident.
    for e in &arr.edges {
        let (Some(lc), Some(rc)) = (e.left, e.right) else {
            continue;
        };
        // `k` indexes `resident` *and* drives `voxel_id`; an index loop is clear.
        #[allow(clippy::needless_range_loop)]
        for k in 0..bands.len() {
            if resident[lc][k] && resident[rc][k] {
                uf_union(&mut vparent, voxel_id(lc, k), voxel_id(rc, k));
            }
        }
    }
    // Vertical adjacency: same cell, consecutive resident bands. The index `ci`
    // also drives `voxel_id`, so an index loop is clearer than an enumerate.
    #[allow(clippy::needless_range_loop)]
    for ci in 0..ncells {
        for k in 1..bands.len() {
            if resident[ci][k] && resident[ci][k - 1] {
                uf_union(&mut vparent, voxel_id(ci, k), voxel_id(ci, k - 1));
            }
        }
    }
    let comp_of = |ci: usize, k: usize, vp: &mut Vec<usize>| uf_find(vp, voxel_id(ci, k));

    // ── wedge tags: resolve non-manifold corner pinches ────────────────────
    // Two material cells that meet a single 2-D *vertex* but are not adjacent
    // through any edge there (the void-void corner pinch: each is separated from
    // the other by a void on both sides) belong to the same global solid yet form
    // a non-manifold pinch — the vertical line through that corner would otherwise
    // carry four wall half-edges on one interned curve. We therefore split the
    // interning locally per corner: each maximal locally-connected bundle of
    // resident voxels around a 2-D vertex (a "wedge") gets its own tag, so its
    // walls and caps intern a distinct corner vertex and vertical edge — two
    // manifold edges instead of one 4-way edge. Ordinary (degree-≤2) vertices have
    // exactly one wedge and tag `0`, so this is a no-op everywhere but the pinch.
    let wedges = Wedges::compute(&arr, &resident, bands.len());

    // ── walls: per arrangement edge, per band ──────────────────────────────
    // A straight edge becomes a planar wall quad; an **arc** edge becomes a
    // cylinder wall (Phase 3c). Both face outward toward the void side.
    for e in &arr.edges {
        let pa = arr.point(e.a);
        let pb = arr.point(e.b);
        for (k, bd) in bands.iter().enumerate() {
            let left = e.left.map(|c| resident[c][k]).unwrap_or(false);
            let right = e.right.map(|c| resident[c][k]).unwrap_or(false);
            if left == right {
                continue;
            }
            let mat_cell = if left {
                e.left.unwrap()
            } else {
                e.right.unwrap()
            };
            let comp = comp_of(mat_cell, k, &mut vparent);
            match e.geom {
                EdgeGeom::Seg => {
                    // Wedge tags for the two corners over `e.a` / `e.b`, taken on
                    // the material side's cell in this band. The corner at each end
                    // is bounded by this very arrangement edge `(e.a, e.b)`.
                    let w0 = wedges.tag(e.a, mat_cell, e.a, e.b, k);
                    let w1 = wedges.tag(e.b, mat_cell, e.a, e.b, k);
                    builder.wall(
                        comp,
                        [pa.x, pa.y],
                        [pb.x, pb.y],
                        bd.z0,
                        bd.z1,
                        /* material_on_left = */ left,
                        w0,
                        w1,
                    );
                }
                EdgeGeom::Arc {
                    center,
                    radius,
                    ccw,
                } => {
                    builder.cylinder_wall(
                        comp,
                        [pa.x, pa.y],
                        [pb.x, pb.y],
                        [center.x, center.y],
                        radius,
                        ccw,
                        bd.z0,
                        bd.z1,
                        /* material_on_left = */ left,
                    );
                }
            }
        }
    }

    // ── interfaces: per cell, per level (band boundaries, incl. bottom/top) ─
    let nbands = bands.len();
    for (ci, cell) in arr.cells.iter().enumerate() {
        for k in 0..=nbands {
            let below = if k == 0 { false } else { resident[ci][k - 1] };
            let above = if k == nbands { false } else { resident[ci][k] };
            if below == above {
                continue;
            }
            let z = level_z(&bands, k);
            // The cap belongs to the resident (material) side; its wedge tags are
            // taken in that side's band so they match the abutting wall rims.
            let mat_band = if below { k - 1 } else { k };
            let comp = comp_of(ci, mat_band, &mut vparent);
            let tag = |v: VertexId, va: VertexId, vb: VertexId| wedges.tag(v, ci, va, vb, mat_band);
            let outer = ring_to_bedges(&arr, &cell.outer, &tag);
            let inners: Vec<Vec<BoundaryEdge2>> = cell
                .inner_rings
                .iter()
                .map(|ring| ring_to_bedges(&arr, ring, &tag))
                .collect();
            builder.interface(comp, &outer, &inners, z, /* up = */ below && !above);
        }
    }

    let brep = builder.finish(tol)?;
    Ok(brep)
}

/// One axial band `[z0, z1]` with each operand's presence flag.
#[derive(Debug, Clone)]
struct Band {
    z0: f64,
    z1: f64,
    /// `present[i]` is `true` when the band lies inside operand `i`'s interval.
    present: Vec<bool>,
}

/// Compute the axial bands.
///
/// Breakpoints are every operand interval endpoint merged within `tol`. Between
/// each pair of consecutive breakpoints a band `[z0, z1]` is created (zero-width
/// bands are skipped). Each band records a `present[i]` flag for every operand:
/// whether the band's midpoint lies inside that operand's axial interval. The
/// keep/discard decision is made by the caller (`build_combined`'s `keep`
/// closure) based on these flags; no pruning is done here.
fn bands(operands: &[PrismOperand], tol: &Tol) -> Vec<Band> {
    let mut bps: Vec<f64> = Vec::with_capacity(operands.len() * 2);
    for o in operands {
        bps.push(o.t0);
        bps.push(o.t1);
    }
    bps.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    bps.dedup_by(|x, y| (*x - *y).abs() <= tol.length);

    let mut out = Vec::new();
    for w in bps.windows(2) {
        let (z0, z1) = (w[0], w[1]);
        if z1 - z0 <= tol.length {
            continue;
        }
        let zm = 0.5 * (z0 + z1);
        let present: Vec<bool> = operands
            .iter()
            .map(|o| zm >= o.t0 - tol.length && zm <= o.t1 + tol.length)
            .collect();
        out.push(Band { z0, z1, present });
    }
    out
}

/// The axial coordinate of band boundary level `k` (0 = bottom of first band).
fn level_z(bands: &[Band], k: usize) -> f64 {
    if k == 0 {
        bands[0].z0
    } else {
        bands[k - 1].z1
    }
}

/// A directed 3-D boundary edge of a deferred face: a straight line or an arc of
/// a circle (Phase 3c). Arc edges carry the [`Circle3`] they lie on so that a
/// cap's boundary arc and the cylinder wall's rim arc reach the **same**
/// interned curve and sibling-pair.
#[derive(Debug, Clone, Copy)]
enum BEdge3 {
    /// A straight edge `a → b`. `wa` / `wb` are **wedge tags** for the two
    /// endpoints (see [`Wedges`]): two coincident endpoints with the *same* wedge
    /// tag intern to one shared vertex, while a non-manifold pinch corner gives the
    /// two material wedges distinct tags so the corner vertex (and the vertical
    /// edge through it) split into two manifold instances. `0` is the default
    /// "ordinary" tag carried by every non-pinch endpoint.
    Line {
        a: Point3,
        b: Point3,
        wa: u32,
        wb: u32,
    },
    /// An arc `a → b` on `circle`, swept from `a_ang` to `b_ang` (radians, the
    /// circle's own parameterisation). Arc endpoints are never pinch corners, so
    /// they carry no wedge tag (always the default `0`).
    Arc {
        circle: Circle3,
        a: Point3,
        b: Point3,
        a_ang: f64,
        b_ang: f64,
    },
}

impl BEdge3 {
    /// A straight edge with the default (ordinary) wedge tags on both endpoints.
    #[inline]
    fn line(a: Point3, b: Point3) -> Self {
        BEdge3::Line { a, b, wa: 0, wb: 0 }
    }

    fn start(&self) -> Point3 {
        match self {
            BEdge3::Line { a, .. } => *a,
            BEdge3::Arc { a, .. } => *a,
        }
    }
}

/// The surface a deferred face lies on.
#[derive(Debug, Clone, Copy)]
enum FaceSurf {
    /// A plane with outward normal `n_out` through a reference point.
    Plane { n_out: Vec3, point: Point3 },
    /// A cylinder surface (built verbatim, no canonicalisation), with the face
    /// sense relating its outward normal to the cylinder's radial-out normal:
    /// `Same` for a solid column (material inside), `Reversed` for a void wall.
    Cylinder { cyl: Cylinder, sense: Sense },
}

/// A deferred face: its outer ring, inner hole rings (each a list of typed
/// boundary edges), and the surface it lies on. Faces are collected first,
/// grouped into connected components, and only *then* turned into B-rep topology
/// — each component interning its own vertices, curves and surfaces. This keeps
/// two solids that merely touch at a corner geometrically **independent**.
struct FaceSpec {
    /// The connected component (solid) this face belongs to.
    comp: usize,
    outer: Vec<BEdge3>,
    holes: Vec<Vec<BEdge3>>,
    surf: FaceSurf,
}

/// Number of axial bands, clamped to at least one for voxel indexing.
#[inline]
fn nbands_for(bands: &[Band]) -> usize {
    bands.len().max(1)
}

/// Undirected key of an arrangement edge, by its endpoint vertex-id pair.
type UEdgeKey = (VertexId, VertexId);

#[inline]
fn uedge(a: VertexId, b: VertexId) -> UEdgeKey {
    if a.0 <= b.0 {
        (a, b)
    } else {
        (b, a)
    }
}

/// Per-corner wedge tags that resolve non-manifold pinch points.
///
/// A boundary **corner** is the meeting of two consecutive ring edges of one cell
/// at a 2-D vertex. A single cell can visit the same vertex through *two distinct
/// corners* (the void-void pinch: the surrounding material's one cell wraps around
/// the touch point on two opposite quadrants), which the validator would reject as
/// a non-manifold pinched edge. We therefore group corners into **wedges**: two
/// corners are in the same wedge when they sit on opposite sides of a shared
/// arrangement edge at the vertex whose two cells are *both resident* in the band
/// (material is continuous across — no wall), or when they are the same corner of
/// the same cell in vertically adjacent resident bands. Each wedge gets a distinct
/// per-vertex tag, so its walls and caps intern an independent corner vertex and
/// vertical edge — two manifold edges in place of one 4-way edge.
///
/// Ordinary vertices have a single wedge and tag `0`, so interning is unchanged
/// away from a pinch.
struct Wedges {
    /// `tags[(v, cell, edge_key, band)] → u32`, the tag of the corner of `cell`
    /// bounded by the arrangement edge `edge_key` at vertex `v` in `band`. Absent
    /// ⇒ tag `0`.
    tags: HashMap<(VertexId, usize, UEdgeKey, usize), u32>,
}

/// A boundary corner: cell, and the two arrangement edges meeting at the vertex.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Corner {
    cell: usize,
    e_prev: UEdgeKey,
    e_cur: UEdgeKey,
}

impl Wedges {
    fn compute(arr: &Arrangement, resident: &[Vec<bool>], nbands: usize) -> Self {
        let nbands = nbands.max(1);

        // 1. Enumerate every boundary corner and index it. A corner is keyed by
        //    (cell, e_prev, e_cur); we also record, per (vertex, cell, incident
        //    edge), which corner that edge belongs to — both of a corner's edges
        //    point back to it, so a wall on an edge (and a cap ring meeting there)
        //    can find its corner.
        let mut corner_index: HashMap<Corner, usize> = HashMap::new();
        let mut corners: Vec<(VertexId, usize)> = Vec::new(); // (vertex, cell) per corner
        let mut edge_to_corner: HashMap<(VertexId, usize, UEdgeKey), usize> = HashMap::new();

        let add_ring =
            |ci: usize,
             ring: &[RingEdge],
             corner_index: &mut HashMap<Corner, usize>,
             corners: &mut Vec<(VertexId, usize)>,
             edge_to_corner: &mut HashMap<(VertexId, usize, UEdgeKey), usize>| {
                let m = ring.len();
                if m == 0 {
                    return;
                }
                for i in 0..m {
                    let prev = ring[i];
                    let cur = ring[(i + 1) % m];
                    // Corner at the shared vertex prev.b == cur.a.
                    debug_assert_eq!(prev.b, cur.a);
                    let v = cur.a;
                    let kp = uedge(prev.a, prev.b);
                    let kc = uedge(cur.a, cur.b);
                    let corner = Corner {
                        cell: ci,
                        e_prev: kp,
                        e_cur: kc,
                    };
                    let idx = *corner_index.entry(corner).or_insert_with(|| {
                        corners.push((v, ci));
                        corners.len() - 1
                    });
                    edge_to_corner.insert((v, ci, kp), idx);
                    edge_to_corner.insert((v, ci, kc), idx);
                }
            };

        for (ci, cell) in arr.cells.iter().enumerate() {
            add_ring(
                ci,
                &cell.outer,
                &mut corner_index,
                &mut corners,
                &mut edge_to_corner,
            );
            for ring in &cell.inner_rings {
                add_ring(
                    ci,
                    ring,
                    &mut corner_index,
                    &mut corners,
                    &mut edge_to_corner,
                );
            }
        }

        let ncorners = corners.len();
        // 2. Union corners into wedges, per band. Node space: corner * nbands + k.
        let node = |corner: usize, k: usize| corner * nbands + k;
        let mut parent: Vec<usize> = (0..ncorners * nbands).collect();

        // Horizontal: across each arrangement edge whose two cells are both
        // resident in band k, the two cells' corners that own that edge are one
        // wedge (no wall between them).
        for e in &arr.edges {
            let (Some(l), Some(r)) = (e.left, e.right) else {
                continue;
            };
            let kedge = uedge(e.a, e.b);
            for &v in &[e.a, e.b] {
                let (Some(&cl), Some(&cr)) = (
                    edge_to_corner.get(&(v, l, kedge)),
                    edge_to_corner.get(&(v, r, kedge)),
                ) else {
                    continue;
                };
                // `k` indexes `resident` and drives `node`; an index loop is clear.
                #[allow(clippy::needless_range_loop)]
                for k in 0..nbands {
                    if resident[l][k] && resident[r][k] {
                        uf_union(&mut parent, node(cl, k), node(cr, k));
                    }
                }
            }
        }
        // Vertical: the same corner in adjacent resident bands.
        for (corner, &(_, ci)) in corners.iter().enumerate() {
            for k in 1..nbands {
                if resident[ci][k] && resident[ci][k - 1] {
                    uf_union(&mut parent, node(corner, k), node(corner, k - 1));
                }
            }
        }

        // 3. Assign per-vertex compact tags. Within one vertex, each distinct
        //    resident wedge-root gets a small id (0, 1, …); ≤1 wedge ⇒ all `0`.
        //    Iterate corners in index order (deterministic) so the assignment is
        //    stable across the wall and cap passes (which only ever *read* `tags`).
        let mut per_vertex_next: HashMap<VertexId, u32> = HashMap::new();
        let mut root_tag: HashMap<usize, u32> = HashMap::new();
        for (corner, &(v, ci)) in corners.iter().enumerate() {
            // `k` indexes `resident` and drives `node`; an index loop is clear.
            #[allow(clippy::needless_range_loop)]
            for k in 0..nbands {
                if !resident[ci][k] {
                    continue;
                }
                let root = uf_find(&mut parent, node(corner, k));
                root_tag.entry(root).or_insert_with(|| {
                    let n = per_vertex_next.entry(v).or_insert(0);
                    let t = *n;
                    *n += 1;
                    t
                });
            }
        }

        // 4. Build the (vertex, cell, edge_key, band) → tag map. Both bounding
        //    edges of a corner resolve to the corner's wedge tag, so a wall (keyed
        //    by its arrangement edge) and a cap ring (keyed by its ring edges) at
        //    the same corner agree. Only non-zero tags are stored.
        let mut tags: HashMap<(VertexId, usize, UEdgeKey, usize), u32> = HashMap::new();
        for (corner_key, &corner) in &corner_index {
            let (v, ci) = corners[corner];
            // `k` indexes `resident` and drives `node`; an index loop is clear.
            #[allow(clippy::needless_range_loop)]
            for k in 0..nbands {
                if !resident[ci][k] {
                    continue;
                }
                let root = uf_find(&mut parent, node(corner, k));
                if let Some(&tag) = root_tag.get(&root) {
                    if tag != 0 {
                        tags.insert((v, ci, corner_key.e_prev, k), tag);
                        tags.insert((v, ci, corner_key.e_cur, k), tag);
                    }
                }
            }
        }

        Self { tags }
    }

    /// The wedge tag of the corner of `cell` bounded by arrangement edge
    /// `(va, vb)` at vertex `v` in `band` (`0` if ordinary).
    #[inline]
    fn tag(&self, v: VertexId, cell: usize, va: VertexId, vb: VertexId, band: usize) -> u32 {
        self.tags
            .get(&(v, cell, uedge(va, vb), band))
            .copied()
            .unwrap_or(0)
    }
}

/// Accumulates deferred face specs from the wall / interface pass.
struct Builder {
    frame: Frame,
    tol: Tol,
    specs: Vec<FaceSpec>,
}

impl Builder {
    fn new(frame: &Frame, tol: Tol) -> Self {
        Self {
            frame: *frame,
            tol,
            specs: Vec::new(),
        }
    }

    /// Lift a frame point at axial height `t` to 3-D.
    #[inline]
    fn lift(&self, xy: [f64; 2], t: f64) -> Point3 {
        self.frame.lift(xy, t)
    }

    /// The circle angle (matching [`Circle3::point_at`]) of a 2-D frame point
    /// about a 2-D centre, normalised to `[0, 2π)`. Because the frame's
    /// `(e1, e2)` basis is exactly the `plane_basis(d)` that `Circle3` uses, the
    /// 2-D polar angle is the circle's own parameter, so the lifted endpoint and
    /// `circle.point_at(angle)` agree.
    #[inline]
    fn circle_angle(center: [f64; 2], p: [f64; 2]) -> f64 {
        (p[1] - center[1])
            .atan2(p[0] - center[0])
            .rem_euclid(std::f64::consts::TAU)
    }

    /// The `(start, end)` angular boundary of a directed arc `pa → pb` on a circle
    /// of `center`, swept in the `ccw` direction.
    ///
    /// # The `0`/`2π` seam normalisation rule (watertightness contract)
    ///
    /// A half-edge and its sibling (the same physical arc traversed the other way)
    /// must carry exactly reversed boundaries `[s, e]` vs `[e, s]`, compared with
    /// an **absolute** tolerance — the validator does *not* compare angles modulo
    /// `2π` (`topo::validate::is_sibling`). So the two arcs that border the same
    /// rim — a cylinder-wall rim arc and a cap-ring arc — must agree on the *same*
    /// numeric angle for each shared endpoint, including its `2π` offset. The rule
    /// that guarantees this:
    ///
    /// 1. The arc's **half-span** `half = ½·|sweep|` and its **midpoint angle**
    ///    `mid` are computed from the two endpoints and the direction. Both are
    ///    direction-independent: reversing `pa↔pb` and `ccw` leaves the physical
    ///    midpoint and the unsigned span unchanged.
    /// 2. `mid` is taken as the polar angle of the arc's actual **midpoint
    ///    point** (via [`circle_angle`], one `atan2(..).rem_euclid(2π)` on the
    ///    identical physical point), *not* as `start + ½·sweep`. The latter sits
    ///    right on the `0`/`2π` discontinuity for a seam-crossing arc, so the wall
    ///    and the cap — computing it from opposite endpoints — landed on opposite
    ///    sides of the seam (one read `≈0`, the other `≈2π`), yielding boundaries
    ///    that differed by `2π` and failed to sibling-pair. Anchoring on the
    ///    shared midpoint point removes that ambiguity entirely.
    /// 3. The boundary is then `mid ± half`, oriented by `ccw`, so a seam-crossing
    ///    arc reads e.g. `[-0.61, 0.61]` and its reverse `[0.61, -0.61]` — exact
    ///    reverses regardless of which endpoint each face started from. The two
    ///    semicircles of a full circle (which share the seam endpoints) still get
    ///    *distinct* pairs (their midpoints differ by π), resolving the otherwise
    ///    ambiguous `MultipleSiblings`.
    ///
    /// [`Circle3::point_at`](crate::primitives::Circle3::point_at) is
    /// `2π`-periodic, so a boundary that leaves `[0, 2π)` still names the correct
    /// point.
    fn arc_angles(center: [f64; 2], pa: [f64; 2], pb: [f64; 2], ccw: bool) -> (f64, f64) {
        let a = Self::circle_angle(center, pa);
        let b_raw = Self::circle_angle(center, pb);
        let s = directed_sweep_angles(a, b_raw, ccw); // signed, in (-2π, 2π)
        let half = 0.5 * s.abs();
        // Midpoint angle of the directed arc, reduced to a canonical `[0, 2π)`
        // representative that is the **same for the forward and reverse
        // traversals**. The two traversals' raw midpoints `a + ½·s` differ by
        // exactly `0` or `2π` (the reverse folds the far endpoint into `[0,2π)`),
        // so reducing modulo `2π` agrees — *except* right on the seam, where one
        // raw value is `≈0` and the other `≈2π`; `seam_canonical_angle` snaps both
        // to the same side so they cannot disagree by a full turn.
        let mid = seam_canonical_angle(a + 0.5 * s);
        let lo = mid - half;
        let hi = mid + half;
        if ccw {
            (lo, hi)
        } else {
            (hi, lo)
        }
    }

    /// A [`Circle3`] for a 2-D circle (centre + radius) lifted to axial height `t`,
    /// with normal `+d`.
    fn circle3(&self, center: [f64; 2], radius: f64, t: f64) -> Circle3 {
        let c = self.lift(center, t);
        Circle3::new_unchecked(c, Unit3::new_unchecked(self.frame.d), radius)
    }

    /// Record a planar face from typed boundary edges with the given outward
    /// normal (deferred; built per-component in [`finish`]).
    fn planar_face(
        &mut self,
        comp: usize,
        outer: Vec<BEdge3>,
        holes: Vec<Vec<BEdge3>>,
        n_out: Vec3,
    ) {
        if outer.len() < 2 {
            return;
        }
        let point = outer[0].start();
        self.specs.push(FaceSpec {
            comp,
            outer,
            holes: holes.into_iter().filter(|h| h.len() >= 2).collect(),
            surf: FaceSurf::Plane { n_out, point },
        });
    }

    /// Emit a vertical wall quad on the 2-D segment `p0→p1` spanning `[z0, z1]`.
    ///
    /// `w0` / `w1` are the **wedge tags** of the corners over `p0` / `p1` (see
    /// [`Wedges`]). The two vertical edges (over `p0`, over `p1`) and the four
    /// corner vertices carry the tag of the 2-D point they sit over, so a
    /// non-manifold pinch corner — where four walls meet on one vertical line —
    /// splits into two manifold edges (one per material wedge) instead of a 4-way
    /// non-manifold edge. The two horizontal rim edges run between `p0` and `p1`,
    /// so each carries `w0` at its `p0` end and `w1` at its `p1` end, exactly
    /// matching the cap-ring edge along the same 2-D segment (which is tagged the
    /// same way), preserving the wall↔cap sibling pairing.
    #[allow(clippy::too_many_arguments)]
    fn wall(
        &mut self,
        comp: usize,
        p0: [f64; 2],
        p1: [f64; 2],
        z0: f64,
        z1: f64,
        material_on_left: bool,
        w0: u32,
        w1: u32,
    ) {
        let b0 = self.lift(p0, z0);
        let b1 = self.lift(p1, z0);
        let t1 = self.lift(p1, z1);
        let t0 = self.lift(p0, z1);
        let edge = b1 - b0;
        let right_normal = edge.cross(self.frame.d);
        // Each lifted corner carries the wedge tag of its 2-D point: b0/t0 over p0
        // (tag w0), b1/t1 over p1 (tag w1).
        if material_on_left {
            // Ring b0 →(rim,z0)→ b1 →(seam,p1)→ t1 →(rim,z1)→ t0 →(seam,p0)→ b0.
            let lines = vec![
                BEdge3::Line {
                    a: b0,
                    b: b1,
                    wa: w0,
                    wb: w1,
                },
                BEdge3::Line {
                    a: b1,
                    b: t1,
                    wa: w1,
                    wb: w1,
                },
                BEdge3::Line {
                    a: t1,
                    b: t0,
                    wa: w1,
                    wb: w0,
                },
                BEdge3::Line {
                    a: t0,
                    b: b0,
                    wa: w0,
                    wb: w0,
                },
            ];
            self.planar_face(comp, lines, Vec::new(), right_normal);
        } else {
            // Ring b0 →(seam,p0)→ t0 →(rim,z1)→ t1 →(seam,p1)→ b1 →(rim,z0)→ b0.
            let lines = vec![
                BEdge3::Line {
                    a: b0,
                    b: t0,
                    wa: w0,
                    wb: w0,
                },
                BEdge3::Line {
                    a: t0,
                    b: t1,
                    wa: w0,
                    wb: w1,
                },
                BEdge3::Line {
                    a: t1,
                    b: b1,
                    wa: w1,
                    wb: w1,
                },
                BEdge3::Line {
                    a: b1,
                    b: b0,
                    wa: w1,
                    wb: w0,
                },
            ];
            self.planar_face(comp, lines, Vec::new(), -right_normal);
        }
    }

    /// Emit a **cylinder wall** patch over the 2-D arc `p0→p1` (on a circle of
    /// `center` / `radius`, swept `ccw`) spanning `[z0, z1]` (Phase 3c).
    ///
    /// The patch is bounded by the bottom rim arc (at `z0`), a vertical seam, the
    /// top rim arc (at `z1`), and a vertical seam — exactly the half-cylinder face
    /// shape the extruder produces, so [`mass::signed_volume`] integrates it in
    /// closed form. The rim arcs lie on [`Circle3`]s shared (by interning) with the
    /// cap interfaces, and the surface is a [`Cylinder`] shared across bands, so
    /// every sibling pairs.
    #[allow(clippy::too_many_arguments)]
    fn cylinder_wall(
        &mut self,
        comp: usize,
        p0: [f64; 2],
        p1: [f64; 2],
        center: [f64; 2],
        radius: f64,
        ccw: bool,
        z0: f64,
        z1: f64,
        material_on_left: bool,
    ) {
        let (a0, a1) = Self::arc_angles(center, p0, p1, ccw);
        let bottom = self.circle3(center, radius, z0);
        let top = self.circle3(center, radius, z1);
        let b0 = self.lift(p0, z0);
        let b1 = self.lift(p1, z0);
        let t1 = self.lift(p1, z1);
        let t0 = self.lift(p0, z1);
        // Cylinder axis line through the 2-D centre, direction +d.
        let axis = Line3::new(self.lift(center, 0.0), self.frame.d).expect("cylinder axis");
        let cyl = Cylinder::new_unchecked(axis, radius);

        // The four rim arcs use the canonical endpoint angles (a0 at p0, a1 at
        // p1), so each pairs with the matching cap-ring arc (same circle, reversed
        // boundary). Loop: b0 →(arc)→ b1 →(seam up)→ t1 →(arc back)→ t0 →(seam
        // down)→ b0 when material is on the left, else reversed.
        let arc_bottom_fwd = BEdge3::Arc {
            circle: bottom,
            a: b0,
            b: b1,
            a_ang: a0,
            b_ang: a1,
        };
        let arc_top_back = BEdge3::Arc {
            circle: top,
            a: t1,
            b: t0,
            a_ang: a1,
            b_ang: a0,
        };
        let arc_bottom_rev = BEdge3::Arc {
            circle: bottom,
            a: b1,
            b: b0,
            a_ang: a1,
            b_ang: a0,
        };
        let arc_top_fwd = BEdge3::Arc {
            circle: top,
            a: t0,
            b: t1,
            a_ang: a0,
            b_ang: a1,
        };
        let outer = if material_on_left {
            vec![
                arc_bottom_fwd,
                BEdge3::line(b1, t1),
                arc_top_back,
                BEdge3::line(t0, b0),
            ]
        } else {
            vec![
                arc_bottom_rev,
                BEdge3::line(b0, t0),
                arc_top_fwd,
                BEdge3::line(t1, b1),
            ]
        };
        // The cylinder surface's stored normal is radial-out. The face's *outward*
        // normal must point to the void. The circle interior is on the left of a
        // CCW arc, so the material occupies the interior iff `material_on_left ==
        // ccw`. A solid column (material inside) faces radial-out (`Same`); a void
        // (material outside) faces radial-in (`Reversed`). This is recorded in the
        // spec so the cap/wall both carry the correct outward sense, which the
        // volume integral relies on.
        let material_inside = material_on_left == ccw;
        let cyl_sense = if material_inside {
            Sense::Same
        } else {
            Sense::Reversed
        };
        self.specs.push(FaceSpec {
            comp,
            outer,
            holes: Vec::new(),
            surf: FaceSurf::Cylinder {
                cyl,
                sense: cyl_sense,
            },
        });
    }

    /// Emit a horizontal interface face for a cell ring at axial level `z`.
    ///
    /// `up` selects an outward normal of `+d` (material below) versus `−d`
    /// (material above). Arc ring edges become arc boundary edges on a shared
    /// [`Circle3`], so the cap's circular boundary pairs with the cylinder wall's
    /// rim arc.
    fn interface(
        &mut self,
        comp: usize,
        outer2: &[BoundaryEdge2],
        inner2: &[Vec<BoundaryEdge2>],
        z: f64,
        up: bool,
    ) {
        let outer = self.bedges_at(outer2, z);
        let holes: Vec<Vec<BEdge3>> = inner2.iter().map(|h| self.bedges_at(h, z)).collect();
        // outer2 is CCW in the (e1, e2) frame, normal +d. A down-facing floor
        // reverses every ring so the outer is CW and holes CCW, matching −d.
        if up {
            self.planar_face(comp, outer, holes, self.frame.d);
        } else {
            let rev_outer = reverse_bedges(&outer);
            let rev_holes: Vec<Vec<BEdge3>> = holes.iter().map(|h| reverse_bedges(h)).collect();
            self.planar_face(comp, rev_outer, rev_holes, -self.frame.d);
        }
    }

    /// Lift a list of 2-D boundary edges to 3-D at axial height `z`.
    fn bedges_at(&self, edges: &[BoundaryEdge2], z: f64) -> Vec<BEdge3> {
        edges
            .iter()
            .map(|e| match *e {
                BoundaryEdge2::Seg { a, b, wa, wb } => BEdge3::Line {
                    a: self.lift(a, z),
                    b: self.lift(b, z),
                    wa,
                    wb,
                },
                BoundaryEdge2::Arc {
                    a,
                    b,
                    center,
                    radius,
                    ccw,
                } => {
                    let circle = self.circle3(center, radius, z);
                    let (a_ang, b_ang) = Self::arc_angles(center, a, b, ccw);
                    BEdge3::Arc {
                        circle,
                        a: self.lift(a, z),
                        b: self.lift(b, z),
                        a_ang,
                        b_ang,
                    }
                }
            })
            .collect()
    }

    /// Build the B-rep from the deferred face specs, one independent solid per
    /// connected component, and validate the result at `Full`.
    fn finish(self, tol: &Tol) -> Result<Brep, PrismError> {
        if self.specs.is_empty() {
            return Ok(Brep::new());
        }
        let mut comp_ids: Vec<usize> = self.specs.iter().map(|s| s.comp).collect();
        comp_ids.sort_unstable();
        comp_ids.dedup();

        let mut brep = Brep::new();
        let mut solids: Vec<Id<Solid>> = Vec::new();
        for comp in comp_ids {
            let mut cb = ComponentBuilder::new(&mut brep, self.tol);
            for spec in self.specs.iter().filter(|s| s.comp == comp) {
                cb.add_face(&spec.outer, &spec.holes, spec.surf);
            }
            let faces = cb.faces;
            if faces.is_empty() {
                continue;
            }
            let shell = brep.topo.add_shell(Shell { faces });
            let solid = brep.topo.add_solid(Solid {
                shells: vec![shell],
            });
            solids.push(solid);
        }
        brep.solids = solids;
        brep.validate(tol, ValidateLevel::Full)
            .map_err(PrismError::InvalidResult)?;
        Ok(brep)
    }
}

/// Reduce an angle to a canonical `[0, 2π)` representative that is **seam-stable**:
/// any value within a hair of a full turn (the `0`/`2π` seam) folds to exactly
/// `0`, not to `2π − ε`.
///
/// This is the keystone of the arc-boundary normalisation
/// ([`Builder::arc_angles`]). The midpoint angle of a directed arc is computed as
/// `start + ½·sweep`; the same physical arc traversed the other way yields a raw
/// midpoint that differs by exactly `0` or `2π`. A plain `rem_euclid(2π)` maps
/// those to the same number *except* when the midpoint sits on the seam, where
/// rounding pushes one copy to `≈0` and the other to `≈2π` — and the resulting
/// boundary pairs then differ by `2π` and fail to sibling-pair (the watertightness
/// break this function fixes). Snapping the seam removes that single ambiguity, so
/// forward and reverse arcs always produce exactly reversed boundaries.
fn seam_canonical_angle(theta: f64) -> f64 {
    use std::f64::consts::TAU;
    // A sub-tolerance angular epsilon: far below any geometric angle the build can
    // produce, but comfortably above floating-point seam jitter (~1e-12).
    const SEAM_EPS: f64 = 1.0e-9;
    let r = theta.rem_euclid(TAU);
    if r >= TAU - SEAM_EPS || r <= SEAM_EPS {
        0.0
    } else {
        r
    }
}

/// Signed sweep (radians, in `(-2π, 2π)`) from angle `a` to angle `b_raw`
/// (both in `[0, 2π)`) in the `ccw` direction.
fn directed_sweep_angles(a: f64, b_raw: f64, ccw: bool) -> f64 {
    use std::f64::consts::TAU;
    if ccw {
        (b_raw - a).rem_euclid(TAU)
    } else {
        -((a - b_raw).rem_euclid(TAU))
    }
}

/// Reverse a list of directed boundary edges (for a down-facing cap).
fn reverse_bedges(edges: &[BEdge3]) -> Vec<BEdge3> {
    edges
        .iter()
        .rev()
        .map(|e| match *e {
            BEdge3::Line { a, b, wa, wb } => BEdge3::Line {
                a: b,
                b: a,
                wa: wb,
                wb: wa,
            },
            BEdge3::Arc {
                circle,
                a,
                b,
                a_ang,
                b_ang,
            } => BEdge3::Arc {
                circle,
                a: b,
                b: a,
                a_ang: b_ang,
                b_ang: a_ang,
            },
        })
        .collect()
}

/// Key identifying a shared circle curve: quantised centre + radius. Two arcs of
/// the same rim circle (a cap arc and a cylinder-wall rim arc) reach the same key
/// and share the curve, so their half-edges sibling-pair.
type CircleKey = (CoordKey, i64);

/// Builds one connected component's faces into a shared [`Brep`], interning its
/// own vertices, curves and surfaces (a fresh namespace per component).
struct ComponentBuilder<'a> {
    brep: &'a mut Brep,
    tol: Tol,
    vert_by_key: HashMap<VKey, Id<Vertex>>,
    coord_by_key: HashMap<VKey, Point3>,
    line_by_key: HashMap<(VKey, VKey), (CurveId, VKey)>,
    circle_by_key: HashMap<CircleKey, CurveId>,
    faces: Vec<Id<Face>>,
}

/// A vertex-identity key: the quantised coordinate **plus a wedge tag**. At a
/// non-manifold corner pinch the two material wedges share one coordinate but
/// carry distinct tags, so they intern as two independent vertices (and the
/// vertical edge through the corner as two independent curves) — splitting the
/// pinch into two manifold edges. Away from a pinch the tag is `0`, so identity
/// reduces to the coordinate as before.
type VKey = (CoordKey, u32);

impl<'a> ComponentBuilder<'a> {
    fn new(brep: &'a mut Brep, tol: Tol) -> Self {
        Self {
            brep,
            tol,
            vert_by_key: HashMap::new(),
            coord_by_key: HashMap::new(),
            line_by_key: HashMap::new(),
            circle_by_key: HashMap::new(),
            faces: Vec::new(),
        }
    }

    fn vertex(&mut self, p: Point3, wedge: u32) -> Id<Vertex> {
        let k: VKey = (key(p), wedge);
        if let Some(&v) = self.vert_by_key.get(&k) {
            return v;
        }
        let pid = self.brep.geom.insert_point(VertexGeom::Explicit(p));
        let v = self.brep.topo.add_vertex(Vertex { point: pid });
        self.vert_by_key.insert(k, v);
        self.coord_by_key.insert(k, p);
        v
    }

    fn line_curve(&mut self, a: Point3, b: Point3, wa: u32, wb: u32) -> (CurveId, VKey) {
        let (ka, kb): (VKey, VKey) = ((key(a), wa), (key(b), wb));
        let unordered = if ka <= kb { (ka, kb) } else { (kb, ka) };
        if let Some(&entry) = self.line_by_key.get(&unordered) {
            return entry;
        }
        let (origin_key, origin_pt, other_pt) = if ka <= kb { (ka, a, b) } else { (kb, b, a) };
        let line = Line3::new(origin_pt, other_pt - origin_pt).expect("non-degenerate edge");
        let cid = self.brep.geom.insert_curve(CurveGeom::Line(line));
        let entry = (cid, origin_key);
        self.line_by_key.insert(unordered, entry);
        entry
    }

    /// Intern a circle curve by its centre + radius key.
    fn circle_curve(&mut self, circle: Circle3) -> CurveId {
        let k: CircleKey = (key(circle.center()), quantize(circle.radius()));
        if let Some(&cid) = self.circle_by_key.get(&k) {
            return cid;
        }
        let cid = self.brep.geom.insert_curve(CurveGeom::Circle(circle));
        self.circle_by_key.insert(k, cid);
        cid
    }

    /// A straight half-edge from `a` to `b`, carrying the endpoints' wedge tags.
    fn line_half_edge(&mut self, a: Point3, b: Point3, wa: u32, wb: u32) -> Id<HalfEdge> {
        let start = self.vertex(a, wa);
        let _ = self.vertex(b, wb);
        let (curve, origin_key) = self.line_curve(a, b, wa, wb);
        let origin_pt = self.coord_by_key[&origin_key];
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

    /// An arc half-edge on a shared circle, with the given angular boundary.
    fn arc_half_edge(
        &mut self,
        circle: Circle3,
        a: Point3,
        b: Point3,
        a_ang: f64,
        b_ang: f64,
    ) -> Id<HalfEdge> {
        // Arc endpoints are never pinch corners; wedge tag `0`.
        let start = self.vertex(a, 0);
        let _ = self.vertex(b, 0);
        let curve = self.circle_curve(circle);
        self.brep.topo.add_half_edge(HalfEdge {
            start,
            curve,
            boundary: [a_ang, b_ang],
        })
    }

    /// Build the half-edges of one boundary ring.
    fn ring_half_edges(&mut self, ring: &[BEdge3]) -> Vec<Id<HalfEdge>> {
        ring.iter()
            .map(|e| match *e {
                BEdge3::Line { a, b, wa, wb } => self.line_half_edge(a, b, wa, wb),
                BEdge3::Arc {
                    circle,
                    a,
                    b,
                    a_ang,
                    b_ang,
                } => self.arc_half_edge(circle, a, b, a_ang, b_ang),
            })
            .collect()
    }

    fn add_face(&mut self, outer: &[BEdge3], holes: &[Vec<BEdge3>], surf: FaceSurf) {
        if outer.len() < 2 {
            return;
        }
        let hes = self.ring_half_edges(outer);
        let lp = self.brep.topo.add_loop(Loop { half_edges: hes });
        let mut inner_loops: Vec<Id<Loop>> = Vec::with_capacity(holes.len());
        for hole in holes {
            if hole.len() < 2 {
                continue;
            }
            let hhes = self.ring_half_edges(hole);
            inner_loops.push(self.brep.topo.add_loop(Loop { half_edges: hhes }));
        }
        let (surface, sense) = match surf {
            FaceSurf::Plane { n_out, point } => {
                let plane = match Plane::new(point, n_out) {
                    Ok(p) => p,
                    Err(_) => return,
                };
                let (surface, flipped) = self.brep.geom.insert_plane(plane, &self.tol);
                (
                    surface,
                    if flipped {
                        Sense::Reversed
                    } else {
                        Sense::Same
                    },
                )
            }
            FaceSurf::Cylinder { cyl, sense } => {
                let surface = self.brep.geom.insert_surface(SurfaceGeom::Cylinder(cyl));
                (surface, sense)
            }
        };
        let f = self.brep.topo.add_face(Face {
            surface,
            sense,
            outer: lp,
            inners: inner_loops,
        });
        self.faces.push(f);
    }
}
