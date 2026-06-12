//! Geometric comparison tolerances and tolerant classification.
//!
//! Building-scale geometry is computed in metres in `f64`, and dimensions
//! range only from about 1e-3 m (plate thickness) to 1e2 m (building span),
//! so a single absolute length tolerance is sufficient — no scale-adaptive
//! tolerancing is needed. This is one of the major simplifications this
//! kernel gains from being domain-specific.
//!
//! Every geometric predicate in the kernel must go through [`Tol`] rather
//! than comparing floats directly. Coincident geometry (a slab face lying
//! exactly on a wall face) is the *common case* in buildings, so the 3-value
//! classification [`Sign3`] is part of the public API from day one.

/// Tolerances used by all geometric predicates.
///
/// `PartialEq` is derived so that a [`Tol`] can be used as part of a cache key
/// (the lazy CSG evaluator keys the cached B-rep on the tolerance it was built
/// with, `DESIGN.md` §5.2).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Tol {
    /// Absolute length tolerance in metres.
    pub length: f64,
    /// Angular tolerance in radians. Also used for comparing cosines and
    /// unit-vector components against 0 or ±1 (the mapping is linear near
    /// 0 and quadratic near ±1, which is acceptable for the parallel /
    /// perpendicular checks it guards).
    pub angular: f64,
}

impl Default for Tol {
    fn default() -> Self {
        Self {
            length: 1e-6,
            angular: 1e-9,
        }
    }
}

/// Tolerant sign of a signed scalar (distance, dot product, …).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Sign3 {
    /// Strictly negative beyond tolerance.
    Below,
    /// Within tolerance of zero.
    On,
    /// Strictly positive beyond tolerance.
    Above,
}

impl Tol {
    /// Classify a signed length against zero with the length tolerance.
    pub fn classify_length(&self, signed: f64) -> Sign3 {
        if signed.abs() <= self.length {
            Sign3::On
        } else if signed > 0.0 {
            Sign3::Above
        } else {
            Sign3::Below
        }
    }

    /// `true` if two lengths are equal within tolerance.
    pub fn eq_length(&self, a: f64, b: f64) -> bool {
        (a - b).abs() <= self.length
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_length_three_ways() {
        let tol = Tol::default();
        assert_eq!(tol.classify_length(1e-3_f64), Sign3::Above);
        assert_eq!(tol.classify_length(-1e-3_f64), Sign3::Below);
        assert_eq!(tol.classify_length(1e-9_f64), Sign3::On);
        assert_eq!(tol.classify_length(-1e-9_f64), Sign3::On);
        assert_eq!(tol.classify_length(0.0_f64), Sign3::On);
    }
}
