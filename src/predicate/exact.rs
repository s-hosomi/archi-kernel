//! Exact-sign predicates (the "exact side" of the two-layer funnel).
//!
//! These wrap Shewchuk's adaptive-precision predicates (the `robust` crate) so
//! that the **exact** half of the predicate semantics described in
//! [`crate::predicate`] has a real implementation, not just a reserved
//! signature. The tolerant [`Sign3`] classification answers "is this *on* the
//! plane as a matter of design intent?"; the exact sign answers "which side is
//! it on, with no round-off mistake?" — used here to orient cap loops and to
//! decide 2-D point-in-polygon containment when building section caps.
//!
//! `robust` is used **only** inside this module (`DESIGN.md` §3.5,
//! `synthesis.md` §1): no other part of the kernel may `use robust`. Everything
//! crossing this boundary is plain `f64`, and the result is the kernel's own
//! [`Sign3`].

use crate::tolerance::Sign3;

/// Exact sign of the orientation of the 2-D triangle `(a, b, c)`.
///
/// Returns [`Sign3::Above`] when `a → b → c` turns counter-clockwise (positive
/// signed area), [`Sign3::Below`] when clockwise, and [`Sign3::On`] when the
/// three points are exactly collinear. The decision is made with Shewchuk's
/// adaptive `orient2d`, so it never mis-signs a near-degenerate triangle the way
/// a naive cross product can.
///
/// Unlike the tolerant layer there is no `On` *band*: `On` here means *exactly*
/// collinear (a true zero of the determinant). Snapping near-coincident input to
/// the same coordinates is the tolerant layer's job; once that is done this
/// gives the combinatorially-consistent sign.
pub fn orient2d_exact(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Sign3 {
    let d = robust::orient2d(
        robust::Coord { x: a[0], y: a[1] },
        robust::Coord { x: b[0], y: b[1] },
        robust::Coord { x: c[0], y: c[1] },
    );
    if d > 0.0 {
        Sign3::Above
    } else if d < 0.0 {
        Sign3::Below
    } else {
        Sign3::On
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ccw_triangle_is_above() {
        let a = [0.0_f64, 0.0_f64];
        let b = [1.0_f64, 0.0_f64];
        let c = [0.0_f64, 1.0_f64];
        assert_eq!(orient2d_exact(a, b, c), Sign3::Above);
    }

    #[test]
    fn cw_triangle_is_below() {
        let a = [0.0_f64, 0.0_f64];
        let b = [0.0_f64, 1.0_f64];
        let c = [1.0_f64, 0.0_f64];
        assert_eq!(orient2d_exact(a, b, c), Sign3::Below);
    }

    #[test]
    fn collinear_is_on() {
        let a = [0.0_f64, 0.0_f64];
        let b = [1.0_f64, 1.0_f64];
        let c = [2.0_f64, 2.0_f64];
        assert_eq!(orient2d_exact(a, b, c), Sign3::On);
    }
}
