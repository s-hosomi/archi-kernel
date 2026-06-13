//! Shared UV-domain helpers for trimmed analytic panels.

use std::f64::consts::PI;

use crate::boolean::poly2d::geom::{Arc, Edge2, Point2};
use crate::boolean::poly2d::intersect::intersect;
use crate::tolerance::Tol;

use super::{CurvedError, TrimEdge2d, TrimLoop2d};

pub(crate) fn validate_range(name: &'static str, min: f64, max: f64) -> Result<(), CurvedError> {
    if !min.is_finite() || !max.is_finite() || max <= min {
        return Err(CurvedError::InvalidRange { name, min, max });
    }
    Ok(())
}

pub(crate) fn angular_step(radius: f64, chord_tol: f64) -> f64 {
    if chord_tol >= radius {
        PI / 8.0
    } else {
        (2.0_f64 * (1.0_f64 - chord_tol / radius).acos()).clamp(PI / 180.0, PI / 8.0)
    }
}

pub(crate) fn parameter_values(
    min: f64,
    max: f64,
    step: f64,
    extra: impl IntoIterator<Item = f64>,
    tol: &Tol,
) -> Vec<f64> {
    let span = max - min;
    let n = ((span / step).ceil() as usize).max(1);
    let mut values = Vec::with_capacity(n + 1);
    for i in 0..=n {
        values.push(min + span * (i as f64) / (n as f64));
    }
    for v in extra {
        if v > min + tol.length && v < max - tol.length && v.is_finite() {
            values.push(v);
        }
    }
    values.sort_by(|a, b| a.total_cmp(b));
    values.dedup_by(|a, b| (*a - *b).abs() <= tol.length);
    values
}

pub(crate) fn validate_holes_in_rect(
    holes: &[TrimLoop2d],
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    tol: &Tol,
) -> Result<(), CurvedError> {
    for hole in holes {
        hole.validate(tol)?;
        let (min, max) = hole.bounds();
        if min[0] < u_min - tol.length
            || max[0] > u_max + tol.length
            || min[1] < v_min - tol.length
            || max[1] > v_max + tol.length
        {
            return Err(CurvedError::HoleOutsidePanel);
        }
    }
    for i in 0..holes.len() {
        for j in (i + 1)..holes.len() {
            if loops_overlap(&holes[i], &holes[j], tol) {
                return Err(CurvedError::HoleOverlap);
            }
        }
    }
    Ok(())
}

pub(crate) fn loops_overlap(a: &TrimLoop2d, b: &TrimLoop2d, tol: &Tol) -> bool {
    if !bounds_overlap(a, b, tol) {
        return false;
    }
    for ea in &a.edges {
        let edge_a = edge2(*ea);
        for eb in &b.edges {
            let edge_b = edge2(*eb);
            let Ok(crossings) = intersect(&edge_a, &edge_b, tol) else {
                // Tangent arc degeneracies mean the holes touch; reject that as
                // overlap rather than accepting an ambiguous zero-width gap.
                return true;
            };
            if !crossings.points.is_empty() {
                return true;
            }
        }
    }
    let pa = a.vertices()[0];
    let pb = b.vertices()[0];
    a.contains_point(pb, tol) || b.contains_point(pa, tol)
}

fn bounds_overlap(a: &TrimLoop2d, b: &TrimLoop2d, tol: &Tol) -> bool {
    let (amin, amax) = a.bounds();
    let (bmin, bmax) = b.bounds();
    amin[0] <= bmax[0] + tol.length
        && amax[0] + tol.length >= bmin[0]
        && amin[1] <= bmax[1] + tol.length
        && amax[1] + tol.length >= bmin[1]
}

fn edge2(edge: TrimEdge2d) -> Edge2 {
    match edge {
        TrimEdge2d::Line { start, end } => {
            Edge2::seg(Point2::new(start[0], start[1]), Point2::new(end[0], end[1]))
        }
        TrimEdge2d::Arc {
            center,
            radius,
            start_angle,
            end_angle,
        } => Edge2::Arc(Arc::new(
            Point2::new(center[0], center[1]),
            radius,
            start_angle,
            end_angle - start_angle,
        )),
    }
}
