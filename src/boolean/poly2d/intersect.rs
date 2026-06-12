//! Closed-form edge intersection with **edge-kind dispatch**.
//!
//! [`intersect`] dispatches on the pair of [`Edge2`] variants. Today only
//! `Seg × Seg` is implemented; the `Seg × Arc` and `Arc × Arc` arms return
//! [`Poly2Error::ArcNotYetSupported`]. This dispatch table is the seam that lets
//! arc support be *added* (one new function per arm) without rewriting any
//! caller: the arrangement asks "where do these two edges cross?" and gets back
//! a list of points, independent of the curve kind.

use crate::boolean::poly2d::error::Poly2Error;
use crate::boolean::poly2d::geom::{eps_sq, orient2d, Arc, Edge2, Orient, Point2};
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
        (Edge2::Seg { start, end }, Edge2::Arc(arc)) => seg_arc(*start, *end, arc, tol),
        (Edge2::Arc(arc), Edge2::Seg { start, end }) => seg_arc(*start, *end, arc, tol),
        (Edge2::Arc(a), Edge2::Arc(b)) => arc_arc(a, b, tol),
    }
}

/// Segment × arc intersection: line × circle closed form, then keep the points
/// that lie within **both** the segment span and the arc's angular sweep.
///
/// Tangency (discriminant ≈ 0 in length units) is normalised to a single
/// contact point. Endpoints that coincide with the arc are recorded too so the
/// arrangement splits at them.
fn seg_arc(s0: Point2, s1: Point2, arc: &Arc, tol: &Tol) -> Result<EdgeCrossings, Poly2Error> {
    let mut out = EdgeCrossings::default();
    let d = s0.to(s1);
    let len = d.len();
    if len <= tol.length {
        // Degenerate segment: treat as the point s0.
        if on_circle(s0, arc, tol) && arc.angle_of_point(s0, tol).is_some() {
            out.points.push(s0);
        }
        return Ok(out);
    }
    let r = arc.radius;
    // Solve |s0 + t·d − c|² = r² for t ∈ [0, 1] (line parameter, not unit).
    let fx = s0.x - arc.center.x;
    let fy = s0.y - arc.center.y;
    let a = d.dot(d);
    let b = 2.0 * (fx * d.x + fy * d.y);
    let cc = fx * fx + fy * fy - r * r;
    let disc = b * b - 4.0 * a * cc;
    // Tangency normalisation: when the perpendicular distance from the centre to
    // the line is within tol of r, treat the line as tangent (one contact).
    // |disc| relates to (dist² − r²); convert the discriminant threshold to a
    // length scale via the foot-of-perpendicular distance.
    let foot_t = -b / (2.0 * a);
    let foot = Point2::new(s0.x + d.x * foot_t, s0.y + d.y * foot_t);
    let foot_dist = foot.dist(arc.center);
    if (foot_dist - r).abs() <= tol.length {
        // Tangent contact: the line grazes the circle at a single point without
        // crossing. A pure tangency does **not** divide either region (it is a
        // measure-zero touch), so we do not emit it as a split point — splitting
        // there would pinch the arrangement and silently drop a face. The regions
        // simply touch; snap handles the (rare) case where the contact coincides
        // with an existing vertex.
        return Ok(out);
    }
    if disc < 0.0 {
        return Ok(out); // no real intersection
    }
    let sq = disc.sqrt();
    for &t in &[(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
        if t >= -tol.length / len && t <= 1.0 + tol.length / len {
            let raw = Point2::new(s0.x + d.x * t, s0.y + d.y * t);
            let p = snap_to_circle(raw, arc);
            if arc.angle_of_point(p, tol).is_some() {
                out.points.push(p);
            }
        }
    }
    Ok(out)
}

/// Arc × arc intersection: two-circle closed form, with same-circle angular
/// overlap and tangency handled explicitly.
fn arc_arc(a: &Arc, b: &Arc, tol: &Tol) -> Result<EdgeCrossings, Poly2Error> {
    let mut out = EdgeCrossings::default();
    let same_center = a.center.coincident(b.center, tol);
    let same_radius = (a.radius - b.radius).abs() <= tol.length;

    if same_center && same_radius {
        // The two arcs lie on the *same circle*. The "intersection" is the angular
        // overlap of their two sweeps. Record the overlap-interval endpoints that
        // fall strictly inside the other arc's sweep, plus shared endpoints (which
        // snapping will dedup), so each arc is split at every shared breakpoint —
        // the circular analogue of collinear segment overlap.
        same_circle_overlap(a, b, tol, &mut out);
        return Ok(out);
    }

    let dx = b.center.x - a.center.x;
    let dy = b.center.y - a.center.y;
    let dist = (dx * dx + dy * dy).sqrt();
    let (r0, r1) = (a.radius, b.radius);

    if dist <= tol.length {
        // Concentric but different radii: never intersect.
        return Ok(out);
    }

    // Tangency. Two sub-cases with very different topology:
    //
    // * **Internal tangency** (`d ≈ |r0 − r1|`): one circle sits inside the other
    //   and touches its boundary at a single point. The contained disc is still a
    //   proper sub-region, so the plane subdivides cleanly *without* a split at
    //   the contact (a pure touch divides nothing). This is the everyday "void
    //   touching the member face" case — we emit no split and let it through.
    //
    // * **External tangency** (`d ≈ r0 + r1`): the two discs touch from outside at
    //   one point. Here the contact point is a genuine **degree-4 pinch**: both
    //   circles share an identical tangent direction there, so neither leaving it
    //   as a non-vertex (the exterior of one disc then wrongly absorbs part of the
    //   other) nor splitting at it (the angular order at the pinch is ambiguous —
    //   the two circles' tangents coincide) yields a correct subdivision. Rather
    //   than return a silently wrong result, we report it as an explicit
    //   unsupported degeneracy (`DESIGN.md` §13-3): the caller can nudge the
    //   geometry off the exact tangent, which is the realistic remedy.
    let external_tan = (dist - (r0 + r1)).abs() <= tol.length;
    let internal_tan = (dist - (r0 - r1).abs()).abs() <= tol.length;
    if external_tan {
        return Err(Poly2Error::UnsupportedArcDegeneracy {
            what: "two circles are externally tangent (degree-4 pinch contact)",
        });
    }
    if internal_tan {
        // Internal tangency: a contained disc touches the enclosing one at a
        // single point. Like external tangency this is a degree-2 pinch contact
        // where the two circles' tangents coincide; the surrounding subdivision
        // and the winding ray-cast cannot be made unambiguous at the contact
        // (a ray grazing it double-counts). Rather than risk a silently collapsed
        // result, we report it explicitly (`DESIGN.md` §13-3) — nudging the inner
        // circle off the exact tangent recovers a normal contained-void case.
        return Err(Poly2Error::UnsupportedArcDegeneracy {
            what: "two circles are internally tangent (pinch contact)",
        });
    }

    if dist > r0 + r1 + tol.length || dist < (r0 - r1).abs() - tol.length {
        return Ok(out); // separate or one strictly inside the other
    }

    // Two intersection points (standard circle-circle):
    // a = (d² + r0² − r1²) / (2d); h = sqrt(r0² − a²).
    let ad = (dist * dist + r0 * r0 - r1 * r1) / (2.0 * dist);
    let h_sq = r0 * r0 - ad * ad;
    if h_sq < 0.0 {
        return Ok(out);
    }
    let h = h_sq.max(0.0).sqrt();
    // Midpoint along the line of centres.
    let mx = a.center.x + ad * dx / dist;
    let my = a.center.y + ad * dy / dist;
    // Perpendicular offset.
    let ox = -dy / dist * h;
    let oy = dx / dist * h;
    for &(px, py) in &[(mx + ox, my + oy), (mx - ox, my - oy)] {
        let p = snap_to_circle(Point2::new(px, py), a);
        if a.angle_of_point(p, tol).is_some() && b.angle_of_point(p, tol).is_some() {
            out.points.push(p);
        }
    }
    Ok(out)
}

/// Record the breakpoints where two arcs on the same circle overlap angularly.
fn same_circle_overlap(a: &Arc, b: &Arc, tol: &Tol, out: &mut EdgeCrossings) {
    // For each endpoint of one arc that lies strictly within the other's sweep,
    // it is a breakpoint. Snapping deduplicates shared endpoints.
    for ep in [b.start(), b.end()] {
        if a.angle_of_point(ep, tol).is_some() {
            out.points.push(snap_to_circle(ep, a));
        }
    }
    for ep in [a.start(), a.end()] {
        if b.angle_of_point(ep, tol).is_some() {
            out.points.push(snap_to_circle(ep, a));
        }
    }
}

/// `true` if `p` lies on the arc's circle within `tol`.
#[inline]
fn on_circle(p: Point2, arc: &Arc, tol: &Tol) -> bool {
    (p.dist(arc.center) - arc.radius).abs() <= tol.length
}

/// Project `p` radially onto the arc's circle (exact radius), so a snapped
/// crossing lands precisely on the curve.
#[inline]
fn snap_to_circle(p: Point2, arc: &Arc) -> Point2 {
    let d = arc.center.to(p);
    let l = d.len();
    if l <= 0.0 {
        return arc.point_at_angle(arc.start_angle);
    }
    Point2::new(
        arc.center.x + d.x / l * arc.radius,
        arc.center.y + d.y / l * arc.radius,
    )
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
    fn seg_arc_two_crossings() {
        // Horizontal line y=0 through a unit circle at origin: hits (±1, 0).
        let tol = Tol::default();
        let a = seg(-2.0_f64, 0.0_f64, 2.0_f64, 0.0_f64);
        // A CCW arc from −π/2 sweeping 3π/2 covers both (1,0) at angle 0 and
        // (−1,0) at angle π.
        let arc = Edge2::Arc(Arc::new(
            Point2::new(0.0_f64, 0.0_f64),
            1.0_f64,
            -std::f64::consts::FRAC_PI_2,
            1.5 * std::f64::consts::PI,
        ));
        let r = intersect(&a, &arc, &tol).unwrap();
        assert!(r.points.iter().any(|p| (p.x - 1.0).abs() <= 1e-9));
        assert!(r.points.iter().any(|p| (p.x + 1.0).abs() <= 1e-9));
    }

    #[test]
    fn seg_arc_tangent_emits_no_split() {
        // Line y=1 tangent to unit circle at (0,1): a pure tangency is not a
        // crossing, so it is not recorded as a split point (it would pinch the
        // arrangement). The regions merely touch.
        let tol = Tol::default();
        let a = seg(-2.0_f64, 1.0_f64, 2.0_f64, 1.0_f64);
        let arc = Edge2::Arc(Arc::new(
            Point2::new(0.0_f64, 0.0_f64),
            1.0_f64,
            0.0_f64,
            std::f64::consts::PI,
        ));
        let r = intersect(&a, &arc, &tol).unwrap();
        assert!(r.points.is_empty(), "tangency emits no split");
    }

    #[test]
    fn arc_arc_two_crossings() {
        use std::f64::consts::PI;
        let tol = Tol::default();
        // Two unit circles centred at (0,0) and (1,0): cross at x=0.5, y=±√0.75.
        let a = Edge2::Arc(Arc::new(
            Point2::new(0.0_f64, 0.0_f64),
            1.0_f64,
            0.0_f64,
            PI,
        ));
        let b = Edge2::Arc(Arc::new(
            Point2::new(1.0_f64, 0.0_f64),
            1.0_f64,
            0.0_f64,
            PI,
        ));
        let r = intersect(&a, &b, &tol).unwrap();
        // Only the upper crossing lies on both upper semicircles.
        assert!(r
            .points
            .iter()
            .any(|p| (p.x - 0.5).abs() <= 1e-9 && (p.y - 0.75_f64.sqrt()).abs() <= 1e-9));
    }

    #[test]
    fn arc_arc_external_tangent_reports_degeneracy() {
        use std::f64::consts::PI;
        let tol = Tol::default();
        // Circles at (0,0) r=1 and (2,0) r=1 touch at (1,0): an external tangent
        // pinch, reported as an explicit degeneracy rather than mis-answered.
        let a = Edge2::Arc(Arc::new(
            Point2::new(0.0_f64, 0.0_f64),
            1.0_f64,
            -PI / 2.0,
            PI,
        ));
        let b = Edge2::Arc(Arc::new(
            Point2::new(2.0_f64, 0.0_f64),
            1.0_f64,
            PI / 2.0,
            PI,
        ));
        assert!(matches!(
            intersect(&a, &b, &tol),
            Err(Poly2Error::UnsupportedArcDegeneracy { .. })
        ));
    }
}
