//! The per-member cache and lazy-evaluation machinery.

use crate::brep::Brep;
use crate::csg::node::CsgNode;
use crate::tolerance::Tol;
use crate::topo::validate::Defect;

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
    /// withheld rather than approximated (`DESIGN.md` §4.1).
    Unsupported3dBoolean,
    /// The evaluation exceeded the configured complexity limit and was isolated
    /// (`DESIGN.md` §4.5).
    ComplexityLimit,
    /// Evaluation produced a structurally invalid B-rep; the defects are
    /// carried for diagnosis.
    InvalidResult(Vec<Defect>),
    /// This operation is not implemented in the current phase.
    NotYetImplemented,
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
    /// The CSG tree — the source of truth for this member.
    pub csg: CsgNode,
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

    /// Mark the cache stale (call after mutating [`csg`](Self::csg)).
    pub fn mark_dirty(&mut self) {
        self.dirty = true;
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
    /// Evaluation is not implemented in this phase, so this currently always
    /// yields [`EvalError::NotYetImplemented`]. The caching contract is in
    /// place: a clean cache built with the same `tol` is returned directly, and
    /// a successful result updates [`last_valid`](Self::last_valid).
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

    /// The actual evaluation. A stub until the build/boolean phases land.
    fn evaluate(&self, _tol: &Tol) -> Result<Brep, EvalError> {
        Err(EvalError::NotYetImplemented)
    }
}
