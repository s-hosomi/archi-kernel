//! Pure topological validation.
//!
//! These checks look only at the combinatorial structure — they never inspect
//! coordinates (the geometric checks live in [`crate::brep`]). The four checks
//! mirror `DESIGN.md` §7:
//!
//! 1. **Euler characteristic, with the ring term.** `V − E + F − (L − F) =
//!    2(S − G)` where `L` is the total number of loops and `E` the number of
//!    sibling pairs. The ring term `(L − F)` is essential: a face with an
//!    interior loop (a wall with a window) breaks the ring-free formula
//!    immediately. The genus `G` is solved from the equation and checked to be
//!    a non-negative integer; an optional `expected_genus` is matched when the
//!    caller knows it.
//! 2. **Sibling-pair completeness (watertight).** Every half-edge must pair
//!    with exactly one sibling: same curve, reversed boundary, with crossing
//!    endpoints.
//! 3. **Loop continuity.** Within a loop the end vertex of each half-edge must
//!    equal the start vertex of the next.
//! 4. **Resolvability.** Every handle must resolve in the store (no dangling
//!    references).

use std::collections::HashMap;

use crate::geom::{CurveId, PointId};
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::{Face, HalfEdge, Loop, Shell, Solid, TopoStore, Vertex};

/// How much of the validation to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum ValidateLevel {
    /// Topology only (Euler, siblings, continuity, resolvability).
    Light,
    /// Topology plus geometric consistency (handled by [`crate::brep`]).
    Full,
}

/// A detected validation failure.
///
/// Defects are self-contained: they carry the handles and, where available,
/// the human-readable quantities (parameter values, counts) needed to diagnose
/// the problem without re-deriving it (`DESIGN.md` §5.2, `synthesis.md` §2-15).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum Defect {
    /// A handle did not resolve in the store (dangling reference).
    DanglingReference {
        /// What kind of entity the handle was meant to name.
        kind: EntityKind,
        /// The unresolved handle's slot index.
        index: u32,
        /// The unresolved handle's generation.
        generation: u32,
    },
    /// The Euler characteristic did not yield a non-negative integer genus.
    EulerCharacteristic {
        /// Vertex count.
        v: i64,
        /// Edge (sibling-pair) count.
        e: i64,
        /// Face count.
        f: i64,
        /// Loop count.
        l: i64,
        /// Shell count.
        s: i64,
        /// The genus implied by the equation (may be fractional or negative,
        /// which is itself the defect).
        implied_genus_times_two: i64,
    },
    /// The implied genus did not match the caller's expectation.
    GenusMismatch {
        /// The genus the caller expected.
        expected: u32,
        /// The genus implied by the Euler characteristic.
        found: i64,
    },
    /// A half-edge has no sibling (the shell is not watertight).
    MissingSibling {
        /// The unpaired half-edge.
        half_edge: Id<HalfEdge>,
        /// The curve it runs along.
        curve: CurveId,
        /// Its boundary interval.
        boundary: [f64; 2],
    },
    /// A half-edge has more than one candidate sibling (local non-manifold).
    MultipleSiblings {
        /// The over-shared half-edge.
        half_edge: Id<HalfEdge>,
        /// The curve it runs along.
        curve: CurveId,
        /// How many candidate siblings were found.
        candidates: usize,
    },
    /// Within a loop, a half-edge's end vertex does not match the next
    /// half-edge's start vertex.
    LoopDiscontinuity {
        /// The loop containing the break.
        loop_id: Id<Loop>,
        /// The half-edge whose end does not meet the next start.
        half_edge: Id<HalfEdge>,
        /// The end vertex implied by this half-edge's sibling pairing.
        expected_start: Option<Id<Vertex>>,
        /// The actual start vertex of the next half-edge.
        actual_start: Id<Vertex>,
    },
    /// A loop has no half-edges.
    EmptyLoop {
        /// The empty loop.
        loop_id: Id<Loop>,
    },
    /// A vertex lies off the surface of a face it bounds (geometric, Full only).
    VertexOffSurface {
        /// The point that is off-surface.
        point: PointId,
        /// Its signed distance from the surface.
        distance: f64,
    },
    /// A half-edge boundary endpoint does not match its vertex (Full only).
    BoundaryVertexMismatch {
        /// The half-edge whose boundary endpoint is off.
        half_edge: Id<HalfEdge>,
        /// The distance between the evaluated boundary point and the vertex.
        distance: f64,
    },
    /// Geometrically, a loop's consecutive edges do not meet (Full only).
    LoopGeometryGap {
        /// The loop containing the gap.
        loop_id: Id<Loop>,
        /// The half-edge whose evaluated end does not meet the next start.
        half_edge: Id<HalfEdge>,
        /// The size of the gap.
        distance: f64,
    },
}

/// The kind of entity a dangling handle was meant to name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum EntityKind {
    /// A vertex.
    Vertex,
    /// A half-edge.
    HalfEdge,
    /// A loop.
    Loop,
    /// A face.
    Face,
    /// A shell.
    Shell,
    /// A solid.
    Solid,
    /// A point geometry handle.
    Point,
}

/// Run the pure topological validation over the given solids.
///
/// `tol` is used only for comparing boundary parameters when pairing siblings;
/// no coordinates are inspected. `expected_genus`, when `Some`, is matched
/// against the genus implied by the Euler characteristic.
///
/// Returns `Ok(())` if no defects are found, otherwise every defect detected.
pub fn validate_topology(
    store: &TopoStore,
    solids: &[Id<Solid>],
    tol: &Tol,
    expected_genus: Option<u32>,
) -> Result<(), Vec<Defect>> {
    let mut defects = Vec::new();

    // Gather the reachable entities. If anything dangles we record it and bail
    // on the counting checks (counts would be meaningless).
    let reachable = match gather(store, solids, &mut defects) {
        Some(r) => r,
        None => return Err(defects),
    };

    let (pairs, sibling_map) = check_siblings(store, &reachable, tol, &mut defects);
    check_loop_continuity(store, &reachable, &sibling_map, &mut defects);
    check_euler(&reachable, pairs, expected_genus, &mut defects);

    if defects.is_empty() {
        Ok(())
    } else {
        Err(defects)
    }
}

/// The set of entities reachable from the validated solids.
struct Reachable {
    vertices: Vec<Id<Vertex>>,
    half_edges: Vec<Id<HalfEdge>>,
    loops: Vec<Id<Loop>>,
    faces: Vec<Id<Face>>,
    shells: Vec<Id<Shell>>,
}

/// Walk `solids → shells → faces → loops → half-edges → vertices`, collecting
/// distinct reachable entities and recording any dangling handle as a defect.
///
/// Returns `None` if any handle failed to resolve, so callers skip the
/// counting checks.
fn gather(store: &TopoStore, solids: &[Id<Solid>], defects: &mut Vec<Defect>) -> Option<Reachable> {
    use std::collections::HashSet;
    let mut vertices: HashSet<Id<Vertex>> = HashSet::new();
    let mut half_edges: HashSet<Id<HalfEdge>> = HashSet::new();
    let mut loops: HashSet<Id<Loop>> = HashSet::new();
    let mut faces: HashSet<Id<Face>> = HashSet::new();
    let mut shells: HashSet<Id<Shell>> = HashSet::new();
    let mut ok = true;

    for &solid_id in solids {
        let Some(solid) = store.solids.get(solid_id) else {
            push_dangling(
                defects,
                EntityKind::Solid,
                solid_id.index(),
                solid_id.generation(),
            );
            ok = false;
            continue;
        };
        for &shell_id in &solid.shells {
            if !shells.insert(shell_id) {
                continue;
            }
            let Some(shell) = store.shells.get(shell_id) else {
                push_dangling(
                    defects,
                    EntityKind::Shell,
                    shell_id.index(),
                    shell_id.generation(),
                );
                ok = false;
                continue;
            };
            for &face_id in &shell.faces {
                if !faces.insert(face_id) {
                    continue;
                }
                let Some(face) = store.faces.get(face_id) else {
                    push_dangling(
                        defects,
                        EntityKind::Face,
                        face_id.index(),
                        face_id.generation(),
                    );
                    ok = false;
                    continue;
                };
                let mut loop_ids = Vec::with_capacity(1 + face.inners.len());
                loop_ids.push(face.outer);
                loop_ids.extend(face.inners.iter().copied());
                for loop_id in loop_ids {
                    if !loops.insert(loop_id) {
                        continue;
                    }
                    let Some(lp) = store.loops.get(loop_id) else {
                        push_dangling(
                            defects,
                            EntityKind::Loop,
                            loop_id.index(),
                            loop_id.generation(),
                        );
                        ok = false;
                        continue;
                    };
                    for &he_id in &lp.half_edges {
                        half_edges.insert(he_id);
                        let Some(he) = store.half_edges.get(he_id) else {
                            push_dangling(
                                defects,
                                EntityKind::HalfEdge,
                                he_id.index(),
                                he_id.generation(),
                            );
                            ok = false;
                            continue;
                        };
                        if store.vertices.get(he.start).is_some() {
                            vertices.insert(he.start);
                        } else {
                            push_dangling(
                                defects,
                                EntityKind::Vertex,
                                he.start.index(),
                                he.start.generation(),
                            );
                            ok = false;
                        }
                    }
                }
            }
        }
    }

    if !ok {
        return None;
    }
    Some(Reachable {
        vertices: vertices.into_iter().collect(),
        half_edges: half_edges.into_iter().collect(),
        loops: loops.into_iter().collect(),
        faces: faces.into_iter().collect(),
        shells: shells.into_iter().collect(),
    })
}

fn push_dangling(defects: &mut Vec<Defect>, kind: EntityKind, index: u32, generation: u32) {
    defects.push(Defect::DanglingReference {
        kind,
        index,
        generation,
    });
}

/// Pair up sibling half-edges and report unpaired / over-shared ones.
///
/// Returns `(pair count, sibling map)`. The pair count is the edge count `E`
/// for Euler; the map sends each successfully paired half-edge to its sibling
/// and is used by the loop-continuity check to derive end vertices. Two
/// half-edges are siblings when they share the same curve and have reversed
/// boundaries (`[a, b]` vs `[b, a]`, compared with `tol.eq_length`). Pairing is
/// greedy within a curve group; a 3-way match is reported as a defect.
fn check_siblings(
    store: &TopoStore,
    reachable: &Reachable,
    tol: &Tol,
    defects: &mut Vec<Defect>,
) -> (usize, HashMap<Id<HalfEdge>, Id<HalfEdge>>) {
    // Group half-edges by the curve they run along.
    let mut by_curve: HashMap<CurveId, Vec<Id<HalfEdge>>> = HashMap::new();
    for &he_id in &reachable.half_edges {
        if let Some(he) = store.half_edges.get(he_id) {
            by_curve.entry(he.curve).or_default().push(he_id);
        }
    }

    let mut pairs = 0usize;
    let mut sibling_map: HashMap<Id<HalfEdge>, Id<HalfEdge>> = HashMap::new();
    for (curve, group) in &by_curve {
        let mut matched = vec![false; group.len()];
        for i in 0..group.len() {
            if matched[i] {
                continue;
            }
            let he_i = store.half_edges.get(group[i]).expect("reachable");
            // Count and collect candidate siblings for `i`.
            let mut candidates: Vec<usize> = Vec::new();
            for (j, &_he_j_id) in group.iter().enumerate().skip(i + 1) {
                if matched[j] {
                    continue;
                }
                let he_j = store.half_edges.get(group[j]).expect("reachable");
                if is_sibling(he_i, he_j, tol) {
                    candidates.push(j);
                }
            }
            match candidates.len() {
                0 => {
                    defects.push(Defect::MissingSibling {
                        half_edge: group[i],
                        curve: *curve,
                        boundary: he_i.boundary,
                    });
                    matched[i] = true;
                }
                1 => {
                    matched[i] = true;
                    matched[candidates[0]] = true;
                    sibling_map.insert(group[i], group[candidates[0]]);
                    sibling_map.insert(group[candidates[0]], group[i]);
                    pairs += 1;
                }
                n => {
                    defects.push(Defect::MultipleSiblings {
                        half_edge: group[i],
                        curve: *curve,
                        candidates: n,
                    });
                    matched[i] = true;
                }
            }
        }
    }
    (pairs, sibling_map)
}

/// `true` if `b` is a valid sibling of `a`: reversed boundary and crossing
/// endpoints (a's start corresponds to b's end and vice versa).
///
/// Endpoint crossing is checked combinatorially through the loop structure via
/// the start-vertex equality of opposite half-edges; here we verify the
/// reversed-boundary condition on the shared curve, which is the geometric
/// signature of a sibling. The vertex-crossing condition is enforced by the
/// loop-continuity check together with this reversal.
fn is_sibling(a: &HalfEdge, b: &HalfEdge, tol: &Tol) -> bool {
    a.curve == b.curve
        && tol.eq_length(a.boundary[0], b.boundary[1])
        && tol.eq_length(a.boundary[1], b.boundary[0])
}

/// Check, purely topologically, that consecutive half-edges in each loop share
/// a vertex.
///
/// The end vertex of a half-edge is not stored, but it equals the start vertex
/// of its sibling (the sibling runs the same curve in reverse, so it starts
/// where this half-edge ends). For each loop we therefore require that
/// `start(sibling(he[i])) == start(he[(i+1) % n])`. When a half-edge has no
/// recorded sibling (already reported by [`check_siblings`]) we skip its check
/// rather than double-reporting. Empty loops are flagged here.
fn check_loop_continuity(
    store: &TopoStore,
    reachable: &Reachable,
    sibling_map: &HashMap<Id<HalfEdge>, Id<HalfEdge>>,
    defects: &mut Vec<Defect>,
) {
    for &loop_id in &reachable.loops {
        let Some(lp) = store.loops.get(loop_id) else {
            continue;
        };
        let n = lp.half_edges.len();
        if n == 0 {
            defects.push(Defect::EmptyLoop { loop_id });
            continue;
        }
        for i in 0..n {
            let cur = lp.half_edges[i];
            let next = lp.half_edges[(i + 1) % n];
            let Some(&sib) = sibling_map.get(&cur) else {
                continue; // no sibling: continuity is undefined, already reported
            };
            let expected_start = store.half_edges.get(sib).map(|he| he.start);
            let actual_start = match store.half_edges.get(next) {
                Some(he) => he.start,
                None => continue,
            };
            if expected_start != Some(actual_start) {
                defects.push(Defect::LoopDiscontinuity {
                    loop_id,
                    half_edge: cur,
                    expected_start,
                    actual_start,
                });
            }
        }
    }
}

/// Check the Euler characteristic with the ring term and the genus.
fn check_euler(
    reachable: &Reachable,
    edge_pairs: usize,
    expected_genus: Option<u32>,
    defects: &mut Vec<Defect>,
) {
    let v = reachable.vertices.len() as i64;
    let e = edge_pairs as i64;
    let f = reachable.faces.len() as i64;
    let l = reachable.loops.len() as i64;
    let s = reachable.shells.len() as i64;

    // V − E + F − (L − F) = 2(S − G)  ⟺  2G = 2S − (V − E + F − (L − F)).
    let chi = v - e + f - (l - f);
    let two_g = 2 * s - chi;

    if two_g < 0 || two_g % 2 != 0 {
        defects.push(Defect::EulerCharacteristic {
            v,
            e,
            f,
            l,
            s,
            implied_genus_times_two: two_g,
        });
        return;
    }
    let g = two_g / 2;
    if let Some(expected) = expected_genus {
        if g != expected as i64 {
            defects.push(Defect::GenusMismatch { expected, found: g });
        }
    }
}
