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
    let arr = Arrangement::build(&regions, tol)?;

    let bands = bands(operands, tol);
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
            // Material is on the side that is resident; the wall faces the void.
            // Edge a→b has its `left` cell on the left. Build the quad so its
            // outward normal points away from the material.
            builder.wall(
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
    for (ci, cell) in arr.cells.iter().enumerate() {
        let ring2: Vec<[f64; 2]> = cell
            .vertex_ids
            .iter()
            .map(|&v| {
                let p = arr.point(v);
                [p.x, p.y]
            })
            .collect();
        for k in 0..=nbands {
            let below = if k == 0 { false } else { resident[ci][k - 1] };
            let above = if k == nbands { false } else { resident[ci][k] };
            if below == above {
                continue;
            }
            let z = level_z(&bands, k);
            // below && !above: top of lower material, outward normal +d (up).
            // above && !below: bottom of upper material, outward normal −d.
            builder.interface(&ring2, z, /* up = */ below && !above);
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

/// Integer-quantised 3-D coordinate key (matches the extruder / cut scale 1e9).
type CoordKey = (i64, i64, i64);

fn key(p: Point3) -> CoordKey {
    let q = |x: f64| (x * 1.0e9_f64).round() as i64;
    (q(p.x), q(p.y), q(p.z))
}

/// Accumulates the 3-D faces with shared vertex / curve interning.
struct Builder {
    frame: Frame,
    tol: Tol,
    brep: Brep,
    vert_by_key: HashMap<CoordKey, Id<Vertex>>,
    coord_by_key: HashMap<CoordKey, Point3>,
    /// Shared straight curves keyed on the unordered endpoint key pair.
    line_by_key: HashMap<(CoordKey, CoordKey), (CurveId, CoordKey)>,
    faces: Vec<Id<Face>>,
}

impl Builder {
    fn new(frame: &Frame, tol: Tol) -> Self {
        Self {
            frame: *frame,
            tol,
            brep: Brep::new(),
            vert_by_key: HashMap::new(),
            coord_by_key: HashMap::new(),
            line_by_key: HashMap::new(),
            faces: Vec::new(),
        }
    }

    /// Lift a frame point at axial height `t` to 3-D.
    #[inline]
    fn lift(&self, xy: [f64; 2], t: f64) -> Point3 {
        self.frame.lift(xy, t)
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

    /// Get-or-create the shared line curve through `a`, `b`, returning the curve
    /// and the coord key chosen as its parameter origin.
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

    /// Add a straight half-edge from `a` to `b` on its shared line.
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

    /// Build a planar face from an ordered 3-D ring with the given outward
    /// normal, registering the plane canonically.
    fn planar_face(&mut self, ring: &[Point3], n_out: Vec3) {
        let n = ring.len();
        if n < 3 {
            return;
        }
        let hes: Vec<Id<HalfEdge>> = (0..n)
            .map(|i| self.half_edge(ring[i], ring[(i + 1) % n]))
            .collect();
        let lp = self.brep.topo.add_loop(Loop { half_edges: hes });
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
            inners: Vec::new(),
        });
        self.faces.push(f);
    }

    /// Emit a vertical wall quad on the 2-D segment `p0→p1` spanning `[z0, z1]`.
    ///
    /// `material_on_left` says the resident cell is to the left of `p0→p1`; the
    /// wall's outward normal must point to the void (right) side.
    fn wall(&mut self, p0: [f64; 2], p1: [f64; 2], z0: f64, z1: f64, material_on_left: bool) {
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
            self.planar_face(&[b0, b1, t1, t0], right_normal);
        } else {
            self.planar_face(&[b0, t0, t1, b1], -right_normal);
        }
    }

    /// Emit a horizontal interface face for a cell ring at axial level `z`.
    ///
    /// `up` selects an outward normal of `+d` (material below, face is the lid)
    /// versus `−d` (material above, face is the floor). The CCW ring is wound so
    /// its normal matches.
    fn interface(&mut self, ring2: &[[f64; 2]], z: f64, up: bool) {
        let ring: Vec<Point3> = ring2.iter().map(|&xy| self.lift(xy, z)).collect();
        // ring2 is CCW in the (e1, e2) frame, whose normal is +d. So a CCW lift
        // has outward normal +d (up). For a down-facing floor, reverse it.
        if up {
            self.planar_face(&ring, self.frame.d);
        } else {
            let rev: Vec<Point3> = ring.iter().rev().copied().collect();
            self.planar_face(&rev, -self.frame.d);
        }
    }

    /// Group the accumulated faces into connected components (by shared edge),
    /// wrap each as its own solid, and validate the result at `Full`.
    fn finish(mut self, tol: &Tol) -> Result<Brep, PrismError> {
        if self.faces.is_empty() {
            return Ok(Brep::new());
        }
        let groups = self.connected_components();
        let mut solids: Vec<Id<Solid>> = Vec::new();
        for group in groups {
            let shell = self.brep.topo.add_shell(Shell { faces: group });
            let solid = self.brep.topo.add_solid(Solid {
                shells: vec![shell],
            });
            solids.push(solid);
        }
        self.brep.solids = solids;
        let brep = std::mem::take(&mut self.brep);
        brep.validate(tol, ValidateLevel::Full)
            .map_err(PrismError::InvalidResult)?;
        Ok(brep)
    }

    /// Partition the faces into connected components by shared (sibling) edge.
    ///
    /// Two faces are adjacent when they share an undirected curve segment — i.e.
    /// a half-edge of one pairs with a sibling half-edge of the other. Walls and
    /// interfaces that border the same voxel share an edge, so a face set that
    /// bounds several disjoint solids splits cleanly here.
    fn connected_components(&self) -> Vec<Vec<Id<Face>>> {
        let nf = self.faces.len();
        // Map an undirected edge (curve, unordered boundary params) to the faces
        // that use it.
        let mut edge_faces: HashMap<(CurveId, i64, i64), Vec<usize>> = HashMap::new();
        for (fi, &face_id) in self.faces.iter().enumerate() {
            let Some(face) = self.brep.topo.faces.get(face_id) else {
                continue;
            };
            let Some(lp) = self.brep.topo.loops.get(face.outer) else {
                continue;
            };
            for &he_id in &lp.half_edges {
                let Some(he) = self.brep.topo.half_edges.get(he_id) else {
                    continue;
                };
                let q = |x: f64| (x * 1.0e9_f64).round() as i64;
                let (a, b) = (q(he.boundary[0]), q(he.boundary[1]));
                let key = (he.curve, a.min(b), a.max(b));
                edge_faces.entry(key).or_default().push(fi);
            }
        }

        // Union-find over faces.
        let mut parent: Vec<usize> = (0..nf).collect();
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for faces in edge_faces.values() {
            for w in faces.windows(2) {
                let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
                if a != b {
                    parent[a] = b;
                }
            }
        }

        let mut groups: HashMap<usize, Vec<Id<Face>>> = HashMap::new();
        for fi in 0..nf {
            let root = find(&mut parent, fi);
            groups.entry(root).or_default().push(self.faces[fi]);
        }
        groups.into_values().collect()
    }
}
