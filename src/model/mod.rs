//! Model layer — a member collection with an inter-member dependency DAG.
//!
//! A [`Model`] owns the building's members keyed by their caller-assigned
//! [`StableId`](crate::csg::StableId) and resolves the *inter-member* booleans
//! that a single [`Member`](crate::csg::Member) cannot express on its own: the
//! priority deductions of quantity take-off (`DESIGN.md` §5.1, §6-2). A girder
//! that is clipped by its columns depends on those columns; when a column moves,
//! the girder must be re-evaluated. That cross-member dirty propagation is what
//! the model adds over the per-member cache.
//!
//! # Dependency DAG
//!
//! Every [`CsgNode::Clip`](crate::csg::CsgNode::Clip) names its clippers by
//! stable id. Each such reference is a dependency edge **clipper → base** (the
//! base depends on the clipper). [`Model::mark_dirty`] walks these edges in the
//! *reverse* direction: marking a column dirty marks every girder that clips
//! against it dirty too, transitively.
//!
//! # Evaluation order and cycles
//!
//! [`Model::evaluate_all`] evaluates members in topological order so that a
//! clipper's geometry is ready before the base that subtracts it. A dependency
//! **cycle** (member A clips B while B clips A) has no topological order; rather
//! than fail the whole model, the members in each cycle are isolated with
//! [`EvalError::CyclicDependency`] and every member *not* in a cycle still
//! evaluates normally (`DESIGN.md` §5.1, local failure isolation).
//!
//! # Clip semantics
//!
//! See [`Model::evaluate`] for the precise residency rule the priority deduction
//! computes (the "keep where in the base and in no clipper" closure) and why it
//! is set-theoretically exact even when several priority levels chain.

use std::collections::BTreeMap;

use crate::boolean::prismatic::{self, ExtrudeLeaf};
use crate::brep::Brep;
use crate::csg::{CsgNode, EvalError, Member, Opening, StableId};
use crate::tolerance::Tol;
use crate::topo::validate::ValidateLevel;

mod dag;

/// A collection of building members plus their inter-member dependency DAG.
///
/// Members are keyed by a caller-assigned [`StableId`]; the kernel never invents
/// ids. Adding a member with an id already present is an error
/// ([`ModelError::DuplicateId`]).
#[derive(Debug, Clone, Default)]
pub struct Model {
    members: BTreeMap<StableId, Member>,
}

/// A failure manipulating the model's membership (distinct from the per-member
/// evaluation failures in [`EvalError`]).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ModelError {
    /// A member with this id is already present; ids must be unique.
    DuplicateId(StableId),
    /// No member with this id is present.
    UnknownId(StableId),
}

impl Model {
    /// Create an empty model.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a member under a caller-assigned id.
    ///
    /// # Errors
    ///
    /// [`ModelError::DuplicateId`] if a member with `id` is already present.
    pub fn insert(&mut self, id: StableId, member: Member) -> Result<(), ModelError> {
        if self.members.contains_key(&id) {
            return Err(ModelError::DuplicateId(id));
        }
        self.members.insert(id, member);
        Ok(())
    }

    /// Remove a member by id, returning it.
    ///
    /// # Errors
    ///
    /// [`ModelError::UnknownId`] if no member with `id` is present.
    pub fn remove(&mut self, id: StableId) -> Result<Member, ModelError> {
        self.members.remove(&id).ok_or(ModelError::UnknownId(id))
    }

    /// Borrow a member by id.
    pub fn get(&self, id: StableId) -> Option<&Member> {
        self.members.get(&id)
    }

    /// Mutably borrow a member by id.
    ///
    /// Mutating the returned member through its own [`csg_mut`](Member::csg_mut)
    /// marks *its* cache stale, but does **not** propagate to dependents; call
    /// [`mark_dirty`](Self::mark_dirty) on the model when a member that others
    /// clip against changes, so the dependent girders/beams are invalidated too.
    pub fn get_mut(&mut self, id: StableId) -> Option<&mut Member> {
        self.members.get_mut(&id)
    }

    /// The number of members in the model.
    pub fn len(&self) -> usize {
        self.members.len()
    }

    /// `true` if the model has no members.
    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    /// The ids of the members, in ascending order.
    pub fn ids(&self) -> impl Iterator<Item = StableId> + '_ {
        self.members.keys().copied()
    }

    /// Mark a member dirty and propagate the invalidation to everything that
    /// depends on it (transitively).
    ///
    /// The dependency edges run clipper → base, so a member's *dependents* are
    /// the bases that clip against it. Marking a column dirty therefore marks
    /// every girder clipping the column dirty, then every member clipping those
    /// girders, and so on — the column-moves-so-the-beam-follows rule
    /// (`DESIGN.md` §5.1). A missing id is a no-op.
    pub fn mark_dirty(&mut self, id: StableId) {
        let dependents = dag::dependents_closure(&self.members, id);
        for dep in dependents {
            if let Some(m) = self.members.get_mut(&dep) {
                m.mark_dirty();
            }
        }
    }

    /// Evaluate the B-rep of a single member, resolving its clip deductions.
    ///
    /// Members that this member clips against must already be evaluable (they
    /// are resolved to their *gross extrusion leaves* here, see below). For a
    /// take-off pipeline call [`evaluate_all`](Self::evaluate_all) first, which
    /// orders the evaluation; this method is the per-member resolver it drives
    /// and the entry point [`takeoff`](crate::model::takeoff) uses.
    ///
    /// # Clip / priority-deduction semantics
    ///
    /// [`ClipRule::Priority`](crate::csg::ClipRule::Priority) means: **deduct the
    /// occupied volume of every clipper from the base**. The kept region is the
    /// flat residency closure
    ///
    /// ```text
    ///   keep(x) = in_base(x) ∧ ¬in_opening₀(x) ∧ … ∧ ¬in_clipper₀(x) ∧ …
    /// ```
    ///
    /// evaluated in one shared prismatic arrangement (`DESIGN.md` §4.2). Two
    /// points make this set-theoretically exact rather than a heuristic:
    ///
    /// * **Clippers are each member's *gross* extrusion**, not its already-clipped
    ///   shape. A girder is clipped by the column's full prism; a small beam is
    ///   clipped by the column's *and* the girder's full prisms. One might worry
    ///   this double-counts where beam ∩ girder ∩ column overlap, but it does
    ///   not: the kept region is the single conjunction above — a voxel removed
    ///   because it lies in the column is simply removed once, regardless of how
    ///   many clippers also cover it. Subtraction of a union is idempotent on
    ///   overlap (`A ∖ (B ∪ C) = (A ∖ B) ∖ C`), so the flat closure is correct.
    /// * **Priority order is expressed by the caller**, by choosing *whose*
    ///   clipper list contains *whom*: column-priority means the girder lists the
    ///   columns as clippers, and the beam lists the columns and girders. The
    ///   rule here imposes no ordering of its own; it only deducts.
    ///
    /// # Errors
    ///
    /// * [`EvalError::UnknownClipper`] if a clipper id is absent from the model.
    /// * Any error the underlying prismatic boolean or extrusion raises
    ///   (unsupported direction, complexity limit, invalid result).
    pub fn evaluate(&self, id: StableId, tol: &Tol) -> Result<Brep, EvalError> {
        let member = self
            .members
            .get(&id)
            .ok_or(EvalError::UnknownClipper { clipper: id.0 })?;
        self.evaluate_node(member.csg(), tol)
    }

    /// Evaluate one CSG node, resolving any [`Clip`](CsgNode::Clip) against the
    /// model's members.
    fn evaluate_node(&self, node: &CsgNode, tol: &Tol) -> Result<Brep, EvalError> {
        match node {
            CsgNode::Clip {
                base,
                clippers,
                rule: _,
            } => {
                let base_leaf = occupancy_leaf(base).ok_or(EvalError::Unsupported3dBoolean {
                    reason: crate::csg::UnsupportedReason::NonLeafOperands,
                })?;
                // Openings declared on the base (an OpeningSubtraction base)
                // subtract alongside the clippers in the same shared pass.
                let opening_leaves = base_opening_leaves(base)?;
                let mut clipper_leaves = Vec::with_capacity(clippers.len());
                for &clip_id in clippers {
                    let clip_member = self
                        .members
                        .get(&clip_id)
                        .ok_or(EvalError::UnknownClipper { clipper: clip_id.0 })?;
                    let leaf = occupancy_leaf(clip_member.csg()).ok_or(
                        EvalError::Unsupported3dBoolean {
                            reason: crate::csg::UnsupportedReason::NonLeafOperands,
                        },
                    )?;
                    clipper_leaves.push(leaf);
                }
                let brep = prismatic::clip(&base_leaf, &opening_leaves, &clipper_leaves, tol)?;
                brep.validate(tol, ValidateLevel::Full)
                    .map_err(EvalError::InvalidResult)?;
                Ok(brep)
            }
            // Any non-clip node is a plain per-member evaluation: delegate to a
            // throwaway Member so the (already validated) single-member logic is
            // reused verbatim, keeping one evaluation path.
            other => {
                let mut m = Member::new(other.clone());
                m.brep(tol).cloned()
            }
        }
    }

    /// Evaluate every member in dependency order, returning a result per id.
    ///
    /// Members are processed so that a clipper is evaluated before any base that
    /// deducts it. Members caught in a dependency cycle are each returned as
    /// [`EvalError::CyclicDependency`] (carrying the cycle's member ids), while
    /// every acyclic member evaluates normally — a local failure isolation, not
    /// a whole-model abort.
    pub fn evaluate_all(&self, tol: &Tol) -> BTreeMap<StableId, Result<Brep, EvalError>> {
        let order = dag::topological_order(&self.members);
        let mut out: BTreeMap<StableId, Result<Brep, EvalError>> = BTreeMap::new();

        // Members in a cycle: isolate each with the cycle membership.
        for cycle in &order.cycles {
            let members: Vec<u64> = cycle.iter().map(|id| id.0).collect();
            for &id in cycle {
                out.insert(
                    id,
                    Err(EvalError::CyclicDependency {
                        members: members.clone(),
                    }),
                );
            }
        }

        // Acyclic members, in topological order.
        for &id in &order.acyclic {
            out.insert(id, self.evaluate(id, tol));
        }
        out
    }
}

/// The gross occupancy [`ExtrudeLeaf`] of a member's CSG tree — the prism it
/// occupies for deduction purposes, ignoring its own openings and clips.
///
/// For a clip deduction we subtract a clipper's *occupied volume*; openings the
/// clipper itself carries do not add material back into the base (a column with a
/// sleeve still displaces the girder over its full prism). So we descend through
/// [`OpeningSubtraction`](CsgNode::OpeningSubtraction) and [`Clip`](CsgNode::Clip)
/// to the underlying extrusion. Returns `None` if the base is not ultimately a
/// single extrusion (outside the 2.5-D fast path of this phase).
fn occupancy_leaf(node: &CsgNode) -> Option<ExtrudeLeaf> {
    match node {
        CsgNode::Extrude {
            profile,
            origin,
            axis,
            length,
        } => Some(ExtrudeLeaf {
            profile: *profile,
            origin: *origin,
            axis: *axis,
            length: *length,
        }),
        CsgNode::CurvedPanel(_) => None,
        CsgNode::OpeningSubtraction { base, .. } => occupancy_leaf(base),
        CsgNode::Clip { base, .. } => occupancy_leaf(base),
        _ => None,
    }
}

/// The opening leaves declared directly on a clip base, when the base is an
/// [`OpeningSubtraction`](CsgNode::OpeningSubtraction). These subtract alongside
/// the clippers. Returns an empty vector for a plain extrusion base.
fn base_opening_leaves(base: &CsgNode) -> Result<Vec<ExtrudeLeaf>, EvalError> {
    match base {
        CsgNode::OpeningSubtraction { openings, .. } => {
            let mut leaves = Vec::with_capacity(openings.len());
            for (_id, Opening { shape }) in openings {
                let leaf = occupancy_leaf(shape).ok_or(EvalError::Unsupported3dBoolean {
                    reason: crate::csg::UnsupportedReason::NonLeafOperands,
                })?;
                leaves.push(leaf);
            }
            Ok(leaves)
        }
        _ => Ok(Vec::new()),
    }
}

// ── quantity take-off ────────────────────────────────────────────────────────

mod takeoff;

pub use takeoff::{takeoff, FormworkArea, QuantityTakeoff};
