//! The per-member cache and lazy-evaluation machinery.

use crate::boolean::prismatic::{self, ExtrudeLeaf, Operand, PrismError};
use crate::brep::Brep;
use crate::build::extrude;
use crate::csg::node::{CsgNode, Opening};
use crate::error::KernelError;
use crate::primitives::Line3;
use crate::tolerance::Tol;
use crate::topo::validate::{Defect, ValidateLevel};

/// A machine-readable evaluation failure.
///
/// Failures are member-local and detectable, never silent data corruption
/// (`DESIGN.md` §2.3, §4.5). Each variant is descriptive enough to tell *which*
/// degeneracy or limit was hit, so the offending member can be isolated.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum EvalError {
    /// A full 3-D boolean was required but is not yet supported; the result is
    /// withheld rather than approximated (`DESIGN.md` §4.1). The `reason`
    /// records *why* the 2.5-D fast path declined, machine-readably.
    Unsupported3dBoolean {
        /// Why the prismatic reduction was not applicable.
        reason: UnsupportedReason,
    },
    /// The evaluation exceeded the configured complexity limit and was isolated
    /// (`DESIGN.md` §4.5). The measure that tripped the budget is carried.
    ComplexityLimit {
        /// The complexity measure that exceeded the budget.
        measure: usize,
        /// The budget it exceeded.
        budget: usize,
    },
    /// Evaluation produced a structurally invalid B-rep; the defects are
    /// carried for diagnosis.
    InvalidResult(Vec<Defect>),
    /// A constructor rejected the inputs (e.g. a non-positive extrusion length
    /// or a degenerate profile). The error's diagnostic message is carried as a
    /// string; the structured [`KernelError`] is not stored because it holds a
    /// `&'static str` field that cannot round-trip through `serde`.
    Construction(String),
    /// A [`Clip`](crate::csg::CsgNode::Clip) referenced a clipper
    /// [`StableId`](crate::csg::StableId) that is not present in the model, so
    /// the deduction cannot be resolved. The missing id is carried.
    UnknownClipper {
        /// The `u64` payload of the clipper id that did not resolve.
        clipper: u64,
    },
    /// This member is part of a dependency cycle (e.g. two members each clip the
    /// other), so it cannot be placed in topological order. Every member in the
    /// cycle is isolated with this error while members outside the cycle still
    /// evaluate normally (`DESIGN.md` §5.1, local failure isolation). The ids of
    /// the members forming the cycle are carried for diagnosis.
    CyclicDependency {
        /// The `u64` payloads of the [`StableId`](crate::csg::StableId)s that
        /// form the dependency cycle this member belongs to, sorted ascending.
        members: Vec<u64>,
    },
    /// This operation is not implemented in the current phase.
    NotYetImplemented,
}

/// Why a member could not be reduced to the 2.5-D prismatic fast path.
///
/// Fine-grained so a failing member can be isolated and reported precisely
/// (`DESIGN.md` §4.5, `synthesis.md` §2-15).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum UnsupportedReason {
    /// The two operands share no common prismatic direction (e.g. two H-sections
    /// crossed at a right angle), so the boolean does not fall to 2.5-D.
    NoCommonDirection,
    /// An operand has a circular cross-section whose curved boundary would land
    /// on the 2-D side as an arc (Phase 3c). The offending operand is recorded.
    CircularOperand {
        /// Which operand is circular.
        operand: Operand,
    },
    /// The 2-D engine encountered an arc it cannot yet handle (Phase 3c).
    ArcNotYetSupported,
    /// A boolean whose operands are not both single extruded leaves (e.g. a
    /// union of three or more solids), which needs a general `Brep × Brep`
    /// boolean outside this phase.
    NonLeafOperands,
    /// A curved panel is present in the CSG tree. Curved panels are tessellated
    /// through `curved`, not evaluated as B-rep solids in this phase.
    CurvedPanel,
}

impl From<KernelError> for EvalError {
    fn from(e: KernelError) -> Self {
        EvalError::Construction(e.to_string())
    }
}

impl From<PrismError> for EvalError {
    fn from(e: PrismError) -> Self {
        match e {
            PrismError::NoCommonDirection => EvalError::Unsupported3dBoolean {
                reason: UnsupportedReason::NoCommonDirection,
            },
            PrismError::DegenerateAxis => {
                EvalError::Construction("extrusion axis is degenerate (zero direction)".to_string())
            }
            PrismError::CircularInvolved { operand } => EvalError::Unsupported3dBoolean {
                reason: UnsupportedReason::CircularOperand { operand },
            },
            PrismError::ArcNotYetSupported => EvalError::Unsupported3dBoolean {
                reason: UnsupportedReason::ArcNotYetSupported,
            },
            PrismError::ComplexityLimit { measure, budget } => {
                EvalError::ComplexityLimit { measure, budget }
            }
            PrismError::InvalidResult(defects) => EvalError::InvalidResult(defects),
            // An internal 2-D failure is a bug surfaced as a string, not a panic.
            PrismError::Poly2(p) => EvalError::Construction(p.to_string()),
        }
    }
}

/// A building member: its CSG tree plus the cached lazily-evaluated B-rep.
///
/// Evaluation is *push-dirty / pull-clean* at member granularity
/// (`DESIGN.md` §5.2). When the tree changes, [`mark_dirty`](Self::mark_dirty)
/// invalidates the cache; [`brep`](Self::brep) rebuilds it on demand. The cache
/// key includes the [`Tol`] it was built with, and the previous valid B-rep is
/// retained so display can continue across a failed re-evaluation.
#[derive(Debug, Clone)]
pub struct Member {
    /// The CSG tree — the source of truth for this member. Private so a mutation
    /// cannot bypass the dirty flag: read it through [`csg`](Self::csg) and
    /// mutate it through [`csg_mut`](Self::csg_mut), which marks the cache stale
    /// automatically (the old `pub csg` field was a `mark_dirty`-forgetting
    /// footgun, `docs/design/progress.md`).
    csg: CsgNode,
    /// Whether the cache is stale.
    dirty: bool,
    /// The cached evaluation result (B-rep or error).
    cache: Option<Result<Brep, EvalError>>,
    /// The last successfully evaluated B-rep, kept as a display fallback.
    prev_valid: Option<Brep>,
    /// The tolerance the cache was built with.
    cached_tol: Option<Tol>,
}

impl Member {
    /// Create a member from its CSG tree. The cache starts empty and dirty.
    pub fn new(csg: CsgNode) -> Self {
        Self {
            csg,
            dirty: true,
            cache: None,
            prev_valid: None,
            cached_tol: None,
        }
    }

    /// Mark the cache stale.
    ///
    /// Usually unnecessary — [`csg_mut`](Self::csg_mut) marks it for you — but
    /// exposed for the model layer, which marks dependents dirty when a member
    /// they clip changes.
    pub fn mark_dirty(&mut self) {
        self.dirty = true;
    }

    /// Borrow the CSG tree (the source of truth for this member).
    pub fn csg(&self) -> &CsgNode {
        &self.csg
    }

    /// Mutably borrow the CSG tree, marking the cache stale.
    ///
    /// Any change reached through this borrow invalidates the cached B-rep, so a
    /// re-`brep` re-evaluates. This is the only mutable access to `csg`, which is
    /// why a `mark_dirty` can never be forgotten.
    pub fn csg_mut(&mut self) -> &mut CsgNode {
        self.dirty = true;
        &mut self.csg
    }

    /// `true` if the cache is stale or absent for the given tolerance.
    pub fn is_dirty(&self, tol: &Tol) -> bool {
        self.dirty || self.cached_tol.as_ref() != Some(tol)
    }

    /// The last successfully evaluated B-rep, if any.
    ///
    /// Used as the display fallback when the current evaluation fails.
    pub fn last_valid(&self) -> Option<&Brep> {
        self.prev_valid.as_ref()
    }

    /// Evaluate (or return the cached) B-rep for the given tolerance.
    ///
    /// Extrusions, prismatic differences, opening subtractions and unions are
    /// evaluated; anything outside the 2.5-D fast path yields a descriptive
    /// [`EvalError`]. The caching contract holds: a clean cache built with the
    /// same `tol` is returned directly, and a successful result updates
    /// [`last_valid`](Self::last_valid).
    pub fn brep(&mut self, tol: &Tol) -> Result<&Brep, EvalError> {
        if !self.is_dirty(tol) {
            if let Some(Ok(_)) = &self.cache {
                // Safe: matched Ok above.
                return Ok(self.cache.as_ref().unwrap().as_ref().unwrap());
            }
            if let Some(Err(e)) = &self.cache {
                return Err(e.clone());
            }
        }

        let result = self.evaluate(tol);
        if let Ok(brep) = &result {
            self.prev_valid = Some(brep.clone());
        }
        self.cache = Some(result);
        self.cached_tol = Some(*tol);
        self.dirty = false;

        match self.cache.as_ref().unwrap() {
            Ok(b) => Ok(b),
            Err(e) => Err(e.clone()),
        }
    }

    /// Evaluate the CSG tree into a B-rep.
    ///
    /// Implemented nodes:
    ///
    /// * [`CsgNode::Extrude`] — builds the extruded solid (Phase 2).
    /// * [`CsgNode::Difference`] of two extruded leaves — the prismatic 2.5-D
    ///   difference (`DESIGN.md` §4.2).
    /// * [`CsgNode::OpeningSubtraction`] of an extruded base by extruded
    ///   openings — the openings are fused with a 2-D union and removed in one
    ///   prismatic difference (`DESIGN.md` §4.5).
    /// * [`CsgNode::Union`] of two extruded leaves — the prismatic union.
    ///
    /// All results validate at [`ValidateLevel::Full`]. Anything outside the
    /// 2.5-D fast path returns [`EvalError::Unsupported3dBoolean`] (a
    /// data-preserving explicit error, never a wrong answer) or
    /// [`EvalError::NotYetImplemented`].
    fn evaluate(&self, tol: &Tol) -> Result<Brep, EvalError> {
        match &self.csg {
            CsgNode::Extrude {
                profile,
                origin,
                axis,
                length,
            } => {
                let line = Line3::new(*origin, *axis)?;
                let brep = extrude(profile, &line, *length, tol)?;
                brep.validate(tol, ValidateLevel::Full)
                    .map_err(EvalError::InvalidResult)?;
                Ok(brep)
            }
            CsgNode::Difference { positive, negative } => {
                let pos = extrude_leaf(positive).ok_or(EvalError::Unsupported3dBoolean {
                    reason: UnsupportedReason::NonLeafOperands,
                })?;
                let neg = extrude_leaf(negative).ok_or(EvalError::Unsupported3dBoolean {
                    reason: UnsupportedReason::NonLeafOperands,
                })?;
                Ok(prismatic::difference(&pos, &neg, tol)?)
            }
            CsgNode::OpeningSubtraction { base, openings } => {
                let base_leaf = extrude_leaf(base).ok_or(EvalError::Unsupported3dBoolean {
                    reason: UnsupportedReason::NonLeafOperands,
                })?;
                let mut opening_leaves = Vec::with_capacity(openings.len());
                for (_id, Opening { shape }) in openings {
                    let leaf = extrude_leaf(shape).ok_or(EvalError::Unsupported3dBoolean {
                        reason: UnsupportedReason::NonLeafOperands,
                    })?;
                    opening_leaves.push(leaf);
                }
                Ok(prismatic::opening_subtraction(
                    &base_leaf,
                    &opening_leaves,
                    tol,
                )?)
            }
            CsgNode::Union(nodes) => {
                let mut leaves = Vec::with_capacity(nodes.len());
                for node in nodes {
                    let leaf = extrude_leaf(node).ok_or(EvalError::Unsupported3dBoolean {
                        reason: UnsupportedReason::NonLeafOperands,
                    })?;
                    leaves.push(leaf);
                }
                Ok(prismatic::union(&leaves, tol)?)
            }
            CsgNode::CurvedPanel(_) => Err(EvalError::Unsupported3dBoolean {
                reason: UnsupportedReason::CurvedPanel,
            }),
            _ => Err(EvalError::NotYetImplemented),
        }
    }
}

/// Extract an [`ExtrudeLeaf`] from a CSG node when it is a single extrusion.
///
/// Returns `None` for any non-leaf node, so the boolean fast path can report
/// [`UnsupportedReason::NonLeafOperands`] cleanly.
fn extrude_leaf(node: &CsgNode) -> Option<ExtrudeLeaf> {
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
        _ => None,
    }
}
