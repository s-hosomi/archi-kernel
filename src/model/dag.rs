//! The inter-member dependency DAG: edge extraction, dirty propagation, and a
//! topological order that isolates cycles.
//!
//! Dependency edges run **clipper → base**: a base that clips against a clipper
//! depends on that clipper's geometry being ready first. The edges are read out
//! of every [`CsgNode::Clip`](crate::csg::CsgNode::Clip) by walking the tree.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::csg::{CsgNode, Member, StableId};

/// Collect the clipper ids referenced anywhere in a member's CSG tree.
///
/// Each id is a dependency: the member depends on that clipper. Walks nested
/// nodes so a clip inside a union/difference is still seen.
pub(crate) fn clipper_ids(node: &CsgNode, out: &mut BTreeSet<StableId>) {
    match node {
        CsgNode::Clip { base, clippers, .. } => {
            for &c in clippers {
                out.insert(c);
            }
            clipper_ids(base, out);
        }
        CsgNode::OpeningSubtraction { base, .. } => clipper_ids(base, out),
        CsgNode::Union(nodes) => {
            for n in nodes {
                clipper_ids(n, out);
            }
        }
        CsgNode::Difference { positive, negative } => {
            clipper_ids(positive, out);
            clipper_ids(negative, out);
        }
        CsgNode::Extrude { .. } => {}
    }
}

/// The direct dependencies (clipper ids) of `id`, filtered to ids actually
/// present in the model (a dangling clipper is surfaced at evaluation time, not
/// here, so it does not perturb ordering).
fn deps_of(members: &BTreeMap<StableId, Member>, id: StableId) -> BTreeSet<StableId> {
    let mut set = BTreeSet::new();
    if let Some(m) = members.get(&id) {
        clipper_ids(m.csg(), &mut set);
    }
    set.retain(|d| members.contains_key(d));
    set
}

/// The transitive closure of *dependents* of `id`: every member that depends on
/// `id` directly or indirectly (excluding `id` itself).
///
/// Used by [`Model::mark_dirty`](crate::model::Model::mark_dirty): the dependents
/// are exactly the members whose cached B-rep is invalidated when `id` changes.
pub(crate) fn dependents_closure(
    members: &BTreeMap<StableId, Member>,
    id: StableId,
) -> BTreeSet<StableId> {
    // Build the reverse adjacency once: dependent_of[c] = bases that clip c.
    let mut rev: BTreeMap<StableId, Vec<StableId>> = BTreeMap::new();
    for &base in members.keys() {
        for dep in deps_of(members, base) {
            rev.entry(dep).or_default().push(base);
        }
    }
    let mut seen = BTreeSet::new();
    let mut queue = VecDeque::new();
    queue.push_back(id);
    while let Some(cur) = queue.pop_front() {
        if let Some(bases) = rev.get(&cur) {
            for &b in bases {
                if seen.insert(b) {
                    queue.push_back(b);
                }
            }
        }
    }
    seen.remove(&id);
    seen
}

/// The result of ordering the members for evaluation.
pub(crate) struct EvalOrder {
    /// Members with no cyclic dependency, in an order where every member's
    /// clippers precede it.
    pub acyclic: Vec<StableId>,
    /// The strongly-connected components of size > 1 (and self-loops): each is a
    /// set of members forming a dependency cycle, sorted ascending.
    pub cycles: Vec<Vec<StableId>>,
}

/// Order the members so that clippers precede the bases that deduct them
/// (Kahn's algorithm on the clipper → base edges), separating out any members
/// caught in a cycle.
///
/// A member is "in a cycle" if it cannot be placed by Kahn's algorithm because
/// its in-degree never reaches zero; the remaining members are grouped into
/// their strongly-connected components for precise reporting.
pub(crate) fn topological_order(members: &BTreeMap<StableId, Member>) -> EvalOrder {
    // Forward edges clipper → base, so a base waits for its clippers.
    // in_degree[base] = number of clippers it depends on.
    let ids: Vec<StableId> = members.keys().copied().collect();
    let mut deps: BTreeMap<StableId, BTreeSet<StableId>> = BTreeMap::new();
    let mut dependents: BTreeMap<StableId, Vec<StableId>> = BTreeMap::new();
    for &id in &ids {
        let d = deps_of(members, id);
        for &c in &d {
            dependents.entry(c).or_default().push(id);
        }
        deps.insert(id, d);
    }

    let mut in_degree: BTreeMap<StableId, usize> =
        ids.iter().map(|&id| (id, deps[&id].len())).collect();

    // Kahn: process zero-in-degree members in ascending id order for determinism.
    let mut ready: BTreeSet<StableId> = ids
        .iter()
        .copied()
        .filter(|id| in_degree[id] == 0)
        .collect();
    let mut acyclic = Vec::new();
    while let Some(&id) = ready.iter().next() {
        ready.remove(&id);
        acyclic.push(id);
        if let Some(bs) = dependents.get(&id) {
            for &b in bs {
                let e = in_degree.get_mut(&b).expect("base in map");
                *e -= 1;
                if *e == 0 {
                    ready.insert(b);
                }
            }
        }
    }

    // Whatever Kahn could not place is part of a cycle. Group the leftovers into
    // strongly-connected components so each reported cycle is the actual set.
    let placed: BTreeSet<StableId> = acyclic.iter().copied().collect();
    let leftover: BTreeSet<StableId> = ids
        .iter()
        .copied()
        .filter(|id| !placed.contains(id))
        .collect();
    let cycles = strongly_connected(&leftover, &deps);

    EvalOrder { acyclic, cycles }
}

/// Tarjan-style SCC over the leftover (cyclic) subgraph, returning each component
/// (sorted) — every component here has at least one cycle by construction.
fn strongly_connected(
    nodes: &BTreeSet<StableId>,
    deps: &BTreeMap<StableId, BTreeSet<StableId>>,
) -> Vec<Vec<StableId>> {
    // Iterative Tarjan to stay panic-free on deep graphs.
    let mut index_counter = 0usize;
    let mut indices: BTreeMap<StableId, usize> = BTreeMap::new();
    let mut lowlink: BTreeMap<StableId, usize> = BTreeMap::new();
    let mut on_stack: BTreeSet<StableId> = BTreeSet::new();
    let mut stack: Vec<StableId> = Vec::new();
    let mut result: Vec<Vec<StableId>> = Vec::new();

    // Edges restricted to the leftover subgraph (clipper → base direction is not
    // needed for SCC; use the dependency edges base → clipper, which give the
    // same SCCs).
    let succ = |n: StableId| -> Vec<StableId> {
        deps.get(&n)
            .map(|s| s.iter().copied().filter(|d| nodes.contains(d)).collect())
            .unwrap_or_default()
    };

    // Explicit work stack: (node, next-successor-index).
    enum Frame {
        Enter(StableId),
        Resume(StableId, usize, Vec<StableId>),
    }

    for &start in nodes {
        if indices.contains_key(&start) {
            continue;
        }
        let mut work: Vec<Frame> = vec![Frame::Enter(start)];
        while let Some(frame) = work.pop() {
            match frame {
                Frame::Enter(v) => {
                    indices.insert(v, index_counter);
                    lowlink.insert(v, index_counter);
                    index_counter += 1;
                    stack.push(v);
                    on_stack.insert(v);
                    let succs = succ(v);
                    work.push(Frame::Resume(v, 0, succs));
                }
                Frame::Resume(v, i, succs) => {
                    if i < succs.len() {
                        let w = succs[i];
                        work.push(Frame::Resume(v, i + 1, succs.clone()));
                        if !indices.contains_key(&w) {
                            work.push(Frame::Enter(w));
                        } else if on_stack.contains(&w) {
                            let lw = indices[&w];
                            let lv = lowlink[&v];
                            lowlink.insert(v, lv.min(lw));
                        }
                        continue;
                    }
                    // Done with v: propagate lowlink to whatever scheduled it, and
                    // if v roots an SCC, pop it.
                    if lowlink[&v] == indices[&v] {
                        let mut comp = Vec::new();
                        loop {
                            let w = stack.pop().expect("non-empty scc stack");
                            on_stack.remove(&w);
                            comp.push(w);
                            if w == v {
                                break;
                            }
                        }
                        comp.sort();
                        // Keep only true cycles: a singleton component is a cycle
                        // only if it has a self-edge.
                        let is_cycle =
                            comp.len() > 1 || (comp.len() == 1 && succ(comp[0]).contains(&comp[0]));
                        if is_cycle {
                            result.push(comp);
                        }
                    }
                    // Propagate v's lowlink upward to its parent frame (the
                    // Resume frame now on top, if any).
                    if let Some(Frame::Resume(parent, _, _)) = work.last() {
                        let p = *parent;
                        let lp = lowlink[&p];
                        let lv = lowlink[&v];
                        lowlink.insert(p, lp.min(lv));
                    }
                }
            }
        }
    }
    result.sort();
    result
}
