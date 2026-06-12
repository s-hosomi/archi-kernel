//! Closed-form edge intersection with **edge-kind dispatch**.
//!
//! [`intersect`] dispatches on the pair of [`Edge2`] variants. Today only
//! `Seg × Seg` is implemented; the `Seg × Arc` and `Arc × Arc` arms return
//! [`Poly2Error::ArcNotYetSupported`]. This dispatch table is the seam that lets
//! arc support be *added* (one new function per arm) without rewriting any
//! caller: the arrangement asks "where do these two edges cross?" and gets back
//! a list of points, independent of the curve kind.

use crate::boolean::poly2d::error::Poly2Error;
use crate::boolean::poly2d::geom::{eps_sq, orient2d, Edge2, Orient, Point2};
use crate::tolerance::Tol;

/// Result of intersecting two edges: the crossing points (already de-duplicated
/// against shared endpoints by the caller's snapping).
///
/// For segments this is at most two points (a transversal crossing gives one;
/// a collinear overlap gives the two overlap endpoints).
#[derive(Debug, Clone, Default)]
pub struct EdgeCrossings {
    /// The intersection points.
    pub points: Vec<Point2>,
}

/// Intersect two edges, returning their crossing points.
///
/// Returns [`Poly2Error::ArcNotYetSupported`] if either edge is an arc.
pub fn intersect(a: &Edge2, b: &Edge2, tol: &Tol) -> Result<EdgeCrossings, Poly2Error> {
    match (a, b) {
        (Edge2::Seg { start: a0, end: a1 }, Edge2::Seg { start: b0, end: b1 }) => {
            Ok(seg_seg(*a0, *a1, *b0, *b1, tol))
        }
        // Arc-involving arms are *added* here; the structure is ready.
        _ => Err(Poly2Error::ArcNotYetSupported),
    }
}

/// Segment × segment intersection.
///
/// Uses the exact [`orient2d`] predicate to decide the topological case
/// (transversal / collinear / disjoint), then computes the crossing coordinate
/// in floating point. The coordinate is only ever *snapped* afterwards, so the
/// floating-point construction here does not threaten combinatorial robustness:
/// the exact predicate decides *whether* there is a crossing, the float decides
/// *roughly where*, and the snap decides the *final* shared vertex.
fn seg_seg(a0: Point2, a1: Point2, b0: Point2, b1: Point2, tol: &Tol) -> EdgeCrossings {
    let mut out = EdgeCrossings::default();

    let o0 = orient2d(a0, a1, b0);
    let o1 = orient2d(a0, a1, b1);
    let o2 = orient2d(b0, b1, a0);
    let o3 = orient2d(b0, b1, a1);

    // Collinear case: both endpoints of one segment are collinear with the other.
    if o0 == Orient::Collinear && o1 == Orient::Collinear {
        // Project onto the dominant axis of segment A and find the overlap.
        collinear_overlap(a0, a1, b0, b1, tol, &mut out);
        return out;
    }

    // Proper / improper transversal crossing: the standard straddle test using
    // exact orientations. Endpoints lying *on* the other segment count (closed
    // segments), which is exactly what we want for building geometry.
    let straddle_ab = orientations_straddle(o0, o1);
    let straddle_ba = orientations_straddle(o2, o3);
    if straddle_ab && straddle_ba {
        if let Some(p) = line_line_point(a0, a1, b0, b1) {
            out.points.push(p);
        }
        return out;
    }

    // No transversal crossing, but an endpoint of one may touch the interior of
    // the other (T-junction): e.g. o2 == Collinear means a0 is on line(b0,b1).
    // Add such touch points so the arrangement splits at them.
    add_touch_if_on_segment(b0, b1, a0, o2, tol, &mut out);
    add_touch_if_on_segment(b0, b1, a1, o3, tol, &mut out);
    add_touch_if_on_segment(a0, a1, b0, o0, tol, &mut out);
    add_touch_if_on_segment(a0, a1, b1, o1, tol, &mut out);

    out
}

/// Do the orientations of the two query points relative to the base line
/// indicate that they straddle it (opposite sides, or one/both on the line)?
#[inline]
fn orientations_straddle(o_first: Orient, o_second: Orient) -> bool {
    matches!(
        (o_first, o_second),
        (Orient::Left, Orient::Right)
            | (Orient::Right, Orient::Left)
            | (Orient::Collinear, _)
            | (_, Orient::Collinear)
    )
}

/// If `p` is collinear with segment `(s0,s1)` (signalled by `o == Collinear`)
/// and falls within its span, record it as a touch point.
fn add_touch_if_on_segment(
    s0: Point2,
    s1: Point2,
    p: Point2,
    o: Orient,
    tol: &Tol,
    out: &mut EdgeCrossings,
) {
    if o != Orient::Collinear {
        return;
    }
    if point_within_segment_span(s0, s1, p, tol) {
        out.points.push(p);
    }
}

/// `true` if `p`, assumed collinear with `(s0,s1)`, lies within the closed span.
fn point_within_segment_span(s0: Point2, s1: Point2, p: Point2, tol: &Tol) -> bool {
    let d = s0.to(s1);
    let len_sq = d.len_sq();
    if len_sq <= eps_sq(tol) {
        // Degenerate base segment: treat as the point s0.
        return p.coincident(s0, tol);
    }
    let t = s0.to(p).dot(d) / len_sq;
    // Allow a small relative slack so endpoints count.
    let slack = tol.length / len_sq.sqrt();
    t >= -slack && t <= 1.0 + slack
}

/// Compute the intersection point of two lines through `(a0,a1)` and `(b0,b1)`.
///
/// Returns `None` if the lines are (numerically) parallel; callers reach this
/// only after the exact straddle test, so a `None` here is a near-parallel
/// grazing case that the snapping will handle via endpoints instead.
fn line_line_point(a0: Point2, a1: Point2, b0: Point2, b1: Point2) -> Option<Point2> {
    let r = a0.to(a1);
    let s = b0.to(b1);
    let denom = r.cross(s);
    if denom == 0.0 {
        return None;
    }
    let t = a0.to(b0).cross(s) / denom;
    Some(a0.lerp(a1, t))
}

/// Collinear overlap of two segments: record the overlap endpoints that lie
/// strictly inside the *other* segment (the shared-endpoint coincidences are
/// handled by snapping, so we only need the interior split points here).
fn collinear_overlap(
    a0: Point2,
    a1: Point2,
    b0: Point2,
    b1: Point2,
    tol: &Tol,
    out: &mut EdgeCrossings,
) {
    // For each endpoint of B that lies within A's span, it splits A; and vice
    // versa. Recording all four (snapping removes duplicates) guarantees both
    // segments get cut at every shared breakpoint, which is what makes collinear
    // overlap dedup possible downstream.
    for &p in &[b0, b1] {
        if point_within_segment_span(a0, a1, p, tol) {
            out.points.push(p);
        }
    }
    for &p in &[a0, a1] {
        if point_within_segment_span(b0, b1, p, tol) {
            out.points.push(p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(ax: f64, ay: f64, bx: f64, by: f64) -> Edge2 {
        Edge2::seg(Point2::new(ax, ay), Point2::new(bx, by))
    }

    #[test]
    fn transversal_cross() {
        let tol = Tol::default();
        let a = seg(0.0_f64, 0.0_f64, 2.0_f64, 2.0_f64);
        let b = seg(0.0_f64, 2.0_f64, 2.0_f64, 0.0_f64);
        let r = intersect(&a, &b, &tol).unwrap();
        assert_eq!(r.points.len(), 1);
        let p = r.points[0];
        assert!((p.x - 1.0_f64).abs() <= 1e-12_f64 && (p.y - 1.0_f64).abs() <= 1e-12_f64);
    }

    #[test]
    fn t_junction_endpoint_on_edge() {
        let tol = Tol::default();
        let a = seg(0.0_f64, 0.0_f64, 4.0_f64, 0.0_f64);
        let b = seg(2.0_f64, 0.0_f64, 2.0_f64, 3.0_f64); // starts on a's interior
        let r = intersect(&a, &b, &tol).unwrap();
        assert!(r
            .points
            .iter()
            .any(|p| (p.x - 2.0_f64).abs() <= 1e-9_f64 && p.y.abs() <= 1e-9_f64));
    }

    #[test]
    fn collinear_overlap_breakpoints() {
        let tol = Tol::default();
        let a = seg(0.0_f64, 0.0_f64, 4.0_f64, 0.0_f64);
        let b = seg(2.0_f64, 0.0_f64, 6.0_f64, 0.0_f64);
        let r = intersect(&a, &b, &tol).unwrap();
        // Breakpoints at x=2 (b0 inside a) and x=4 (a1 inside b).
        assert!(r.points.iter().any(|p| (p.x - 2.0_f64).abs() <= 1e-9_f64));
        assert!(r.points.iter().any(|p| (p.x - 4.0_f64).abs() <= 1e-9_f64));
    }

    #[test]
    fn parallel_disjoint_no_crossing() {
        let tol = Tol::default();
        let a = seg(0.0_f64, 0.0_f64, 4.0_f64, 0.0_f64);
        let b = seg(0.0_f64, 1.0_f64, 4.0_f64, 1.0_f64);
        let r = intersect(&a, &b, &tol).unwrap();
        assert!(r.points.is_empty());
    }

    #[test]
    fn arc_returns_error() {
        use crate::boolean::poly2d::geom::Arc;
        let tol = Tol::default();
        let a = seg(0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64);
        let arc = Edge2::Arc(Arc::new(
            Point2::new(0.0_f64, 0.0_f64),
            1.0_f64,
            0.0_f64,
            std::f64::consts::PI,
        ));
        assert!(matches!(
            intersect(&a, &arc, &tol),
            Err(Poly2Error::ArcNotYetSupported)
        ));
    }
}
