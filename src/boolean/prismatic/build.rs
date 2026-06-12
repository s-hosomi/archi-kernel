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

use crate::boolean::poly2d::Op;
use crate::boolean::support::{key, CoordKey};
use crate::brep::Brep;
use crate::geom::{CurveGeom, CurveId, VertexGeom};
use crate::math::{Point3, Vec3};
use crate::primitives::{Line3, Plane};
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::validate::ValidateLevel;
use crate::topo::{Face, HalfEdge, Loop, Sense, Shell, Solid, Vertex};

use super::arrange::Arrangement;
use super::detect::{Frame, PrismOperand};
use super::error::PrismError;

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

    // ── walls: per arrangement edge, per band ──────────────────────────────
    for e in &arr.edges {
        let pa = arr.point(e.a);
        let pb = arr.point(e.b);
        for (k, bd) in bands.iter().enumerate() {
            let left = e.left.map(|c| resident[c][k]).unwrap_or(false);
            let right = e.right.map(|c| resident[c][k]).unwrap_or(false);
            if left == right {
                continue;
            }
            // The wall belongs to the resident (material) side's voxel.
            let comp = if left {
                comp_of(e.left.unwrap(), k, &mut vparent)
            } else {
                comp_of(e.right.unwrap(), k, &mut vparent)
            };
            // Material is on the side that is resident; the wall faces the void.
            // Edge a→b has its `left` cell on the left. Build the quad so its
            // outward normal points away from the material.
            builder.wall(
                comp,
                [pa.x, pa.y],
                [pb.x, pb.y],
                bd.z0,
                bd.z1,
                /* material_on_left = */ left,
            );
        }
    }

    // ── interfaces: per cell, per level (band boundaries, incl. bottom/top) ─
    let nbands = bands.len();
    let to_xy = |v| {
        let p = arr.point(v);
        [p.x, p.y]
    };
    for (ci, cell) in arr.cells.iter().enumerate() {
        let ring2: Vec<[f64; 2]> = cell.vertex_ids.iter().copied().map(to_xy).collect();
        // Inner hole rings of the cell (an annulus interface), so a contained
        // tool's cap is the parent's polygon *with the hole*, not the full
        // polygon (which would re-cover the void).
        let inner2: Vec<Vec<[f64; 2]>> = cell
            .inner_rings
            .iter()
            .map(|ring| ring.iter().copied().map(to_xy).collect())
            .collect();
        for k in 0..=nbands {
            let below = if k == 0 { false } else { resident[ci][k - 1] };
            let above = if k == nbands { false } else { resident[ci][k] };
            if below == above {
                continue;
            }
            let z = level_z(&bands, k);
            // The cap belongs to the resident voxel it bounds: the band below
            // (k−1) when material is below, else the band above (k).
            let comp = if below {
                comp_of(ci, k - 1, &mut vparent)
            } else {
                comp_of(ci, k, &mut vparent)
            };
            // below && !above: top of lower material, outward normal +d (up).
            // above && !below: bottom of upper material, outward normal −d.
            builder.interface(comp, &ring2, &inner2, z, /* up = */ below && !above);
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
/// Breakpoints are every operand interval endpoint merged within `tol`; each
/// resulting band records which operands it lies inside. Bands are kept whenever
/// the first operand (the base / positive side) is present, since the combining
/// rule only ever keeps material there — an empty band would emit no faces
/// anyway, so this prune is just an optimisation.
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

/// A deferred planar face: its outer ring, inner hole rings, and outward normal,
/// all in 3-D world coordinates. Faces are collected first, grouped into
/// connected components, and only *then* turned into B-rep topology — each
/// component interning its own vertices and curves. This keeps two solids that
/// merely touch at a corner (a checkerboard subtraction) geometrically
/// **independent**: their shared corner edge is interned twice (once per solid),
/// so no curve carries the four half-edges that would otherwise make the corner
/// non-manifold and trip the validator's sibling pairing.
struct FaceSpec {
    /// The connected component (solid) this face belongs to.
    comp: usize,
    outer: Vec<Point3>,
    holes: Vec<Vec<Point3>>,
    n_out: Vec3,
}

/// Number of axial bands, clamped to at least one for voxel indexing.
#[inline]
fn nbands_for(bands: &[Band]) -> usize {
    bands.len().max(1)
}

/// Union-find find with path halving.
fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// Union-find union.
fn uf_union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (uf_find(parent, a), uf_find(parent, b));
    if ra != rb {
        parent[ra] = rb;
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

    /// Record a planar face from an ordered 3-D ring with the given outward
    /// normal (deferred; built per-component in [`finish`]).
    fn planar_face(&mut self, comp: usize, ring: &[Point3], n_out: Vec3) {
        self.planar_face_with_holes(comp, ring, &[], n_out);
    }

    /// Record a planar face with an outer ring and any inner hole rings
    /// (deferred). Holes are already wound opposite to the outer ring for
    /// `n_out`.
    fn planar_face_with_holes(
        &mut self,
        comp: usize,
        ring: &[Point3],
        holes: &[Vec<Point3>],
        n_out: Vec3,
    ) {
        if ring.len() < 3 {
            return;
        }
        self.specs.push(FaceSpec {
            comp,
            outer: ring.to_vec(),
            holes: holes.iter().filter(|h| h.len() >= 3).cloned().collect(),
            n_out,
        });
    }

    /// Emit a vertical wall quad on the 2-D segment `p0→p1` spanning `[z0, z1]`.
    ///
    /// `material_on_left` says the resident cell is to the left of `p0→p1`; the
    /// wall's outward normal must point to the void (right) side.
    #[allow(clippy::too_many_arguments)]
    fn wall(
        &mut self,
        comp: usize,
        p0: [f64; 2],
        p1: [f64; 2],
        z0: f64,
        z1: f64,
        material_on_left: bool,
    ) {
        let b0 = self.lift(p0, z0);
        let b1 = self.lift(p1, z0);
        let t1 = self.lift(p1, z1);
        let t0 = self.lift(p0, z1);
        // Quad b0→b1→t1→t0 has, with d the up axis, an outward normal of
        // (edge × d) for the right side of p0→p1. If material is on the left,
        // the void is on the right and that is the correct outward direction;
        // otherwise reverse the loop so the normal flips.
        let edge = b1 - b0;
        let right_normal = edge.cross(self.frame.d);
        if material_on_left {
            self.planar_face(comp, &[b0, b1, t1, t0], right_normal);
        } else {
            self.planar_face(comp, &[b0, t0, t1, b1], -right_normal);
        }
    }

    /// Emit a horizontal interface face for a cell ring at axial level `z`.
    ///
    /// `up` selects an outward normal of `+d` (material below, face is the lid)
    /// versus `−d` (material above, face is the floor). The CCW ring is wound so
    /// its normal matches.
    fn interface(
        &mut self,
        comp: usize,
        ring2: &[[f64; 2]],
        inner2: &[Vec<[f64; 2]>],
        z: f64,
        up: bool,
    ) {
        let ring: Vec<Point3> = ring2.iter().map(|&xy| self.lift(xy, z)).collect();
        // The cell's inner rings arrive CW in the (e1, e2) frame (hole
        // boundaries), which is the correct *opposite-to-outer* hole sense for an
        // outer-CCW / normal-+d face — so the inner-ring edges are
        // reversed-coincident with the hole's vertical wall edges and sibling
        // pairing closes (the watertightness contract).
        let holes: Vec<Vec<Point3>> = inner2
            .iter()
            .map(|h| h.iter().map(|&xy| self.lift(xy, z)).collect())
            .collect();
        // ring2 is CCW in the (e1, e2) frame, whose normal is +d. So a CCW lift
        // has outward normal +d (up). For a down-facing floor, reverse every
        // ring so the outer is CW and the holes CCW, matching the −d normal.
        if up {
            self.planar_face_with_holes(comp, &ring, &holes, self.frame.d);
        } else {
            let rev: Vec<Point3> = ring.iter().rev().copied().collect();
            let rev_holes: Vec<Vec<Point3>> = holes
                .iter()
                .map(|h| h.iter().rev().copied().collect())
                .collect();
            self.planar_face_with_holes(comp, &rev, &rev_holes, -self.frame.d);
        }
    }

    /// Build the B-rep from the deferred face specs, one independent solid per
    /// connected component, and validate the result at `Full`.
    ///
    /// Each component interns its **own** vertices and curves, so two solids that
    /// touch at a corner (a checkerboard subtraction) do not share the corner
    /// curve: every curve carries exactly two half-edges within its component,
    /// and sibling pairing closes per solid.
    fn finish(self, tol: &Tol) -> Result<Brep, PrismError> {
        if self.specs.is_empty() {
            return Ok(Brep::new());
        }
        // Group spec indices by their component id (deterministic order).
        let mut comp_ids: Vec<usize> = self.specs.iter().map(|s| s.comp).collect();
        comp_ids.sort_unstable();
        comp_ids.dedup();

        let mut brep = Brep::new();
        let mut solids: Vec<Id<Solid>> = Vec::new();
        for comp in comp_ids {
            let mut cb = ComponentBuilder::new(&mut brep, self.tol);
            for spec in self.specs.iter().filter(|s| s.comp == comp) {
                cb.add_face(&spec.outer, &spec.holes, spec.n_out);
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

/// Builds one connected component's faces into a shared [`Brep`], interning its
/// own vertices and curves (a fresh namespace per component).
struct ComponentBuilder<'a> {
    brep: &'a mut Brep,
    tol: Tol,
    vert_by_key: HashMap<CoordKey, Id<Vertex>>,
    coord_by_key: HashMap<CoordKey, Point3>,
    line_by_key: HashMap<(CoordKey, CoordKey), (CurveId, CoordKey)>,
    faces: Vec<Id<Face>>,
}

impl<'a> ComponentBuilder<'a> {
    fn new(brep: &'a mut Brep, tol: Tol) -> Self {
        Self {
            brep,
            tol,
            vert_by_key: HashMap::new(),
            coord_by_key: HashMap::new(),
            line_by_key: HashMap::new(),
            faces: Vec::new(),
        }
    }

    fn vertex(&mut self, p: Point3) -> Id<Vertex> {
        let k = key(p);
        if let Some(&v) = self.vert_by_key.get(&k) {
            return v;
        }
        let pid = self.brep.geom.insert_point(VertexGeom::Explicit(p));
        let v = self.brep.topo.add_vertex(Vertex { point: pid });
        self.vert_by_key.insert(k, v);
        self.coord_by_key.insert(k, p);
        v
    }

    fn line_curve(&mut self, a: Point3, b: Point3) -> (CurveId, CoordKey) {
        let (ka, kb) = (key(a), key(b));
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

    fn half_edge(&mut self, a: Point3, b: Point3) -> Id<HalfEdge> {
        let start = self.vertex(a);
        let _ = self.vertex(b);
        let (curve, origin_key) = self.line_curve(a, b);
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

    fn add_face(&mut self, ring: &[Point3], holes: &[Vec<Point3>], n_out: Vec3) {
        let n = ring.len();
        if n < 3 {
            return;
        }
        let hes: Vec<Id<HalfEdge>> = (0..n)
            .map(|i| self.half_edge(ring[i], ring[(i + 1) % n]))
            .collect();
        let lp = self.brep.topo.add_loop(Loop { half_edges: hes });
        let mut inner_loops: Vec<Id<Loop>> = Vec::with_capacity(holes.len());
        for hole in holes {
            let m = hole.len();
            if m < 3 {
                continue;
            }
            let hhes: Vec<Id<HalfEdge>> = (0..m)
                .map(|i| self.half_edge(hole[i], hole[(i + 1) % m]))
                .collect();
            inner_loops.push(self.brep.topo.add_loop(Loop { half_edges: hhes }));
        }
        let plane = match Plane::new(ring[0], n_out) {
            Ok(p) => p,
            Err(_) => return,
        };
        let (surface, flipped) = self.brep.geom.insert_plane(plane, &self.tol);
        let sense = if flipped {
            Sense::Reversed
        } else {
            Sense::Same
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
