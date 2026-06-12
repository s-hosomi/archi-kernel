//! Face classification: the shared in/out classifier and the per-operation
//! keep/drop truth tables.
//!
//! Every boolean operation uses the **same** classifier — a face is selected
//! purely from `(inside_a, inside_b)`. Only the truth table differs. This is the
//! design the parent kernel mandates (DESIGN.md §4.4): "the in/out classifier is
//! shared across difference / intersection / union; only the keep table
//! changes."
//!
//! `inside` is derived from the winding number being non-zero. For the
//! well-formed regions this engine accepts (outer CCW, holes CW, no
//! self-intersection) the non-zero rule and the even-odd rule agree, and
//! non-zero is the natural fit for "outer minus hole" nesting.

/// The boolean operation to perform.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// `A − B`: keep points in A but not in B.
    Difference,
    /// `A ∪ B`: keep points in A or B.
    Union,
    /// `A ∩ B`: keep points in both A and B.
    Intersection,
}

impl Op {
    /// Keep/drop decision for a face given its insideness in each operand.
    ///
    /// This is the truth table. It is the *only* thing that differs between the
    /// three operations; the classifier feeding it is shared.
    #[inline]
    pub fn keep(self, inside_a: bool, inside_b: bool) -> bool {
        match self {
            Op::Difference => inside_a && !inside_b,
            Op::Union => inside_a || inside_b,
            Op::Intersection => inside_a && inside_b,
        }
    }
}

/// `true` if a winding number means "inside" under the non-zero rule.
#[inline]
pub fn inside(winding: i32) -> bool {
    winding != 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn difference_table() {
        assert!(Op::Difference.keep(true, false));
        assert!(!Op::Difference.keep(true, true));
        assert!(!Op::Difference.keep(false, true));
        assert!(!Op::Difference.keep(false, false));
    }

    #[test]
    fn union_table() {
        assert!(Op::Union.keep(true, false));
        assert!(Op::Union.keep(false, true));
        assert!(Op::Union.keep(true, true));
        assert!(!Op::Union.keep(false, false));
    }

    #[test]
    fn intersection_table() {
        assert!(Op::Intersection.keep(true, true));
        assert!(!Op::Intersection.keep(true, false));
        assert!(!Op::Intersection.keep(false, true));
        assert!(!Op::Intersection.keep(false, false));
    }

    #[test]
    fn inside_rule() {
        assert!(!inside(0));
        assert!(inside(1));
        assert!(inside(-1));
        assert!(inside(2));
    }
}
