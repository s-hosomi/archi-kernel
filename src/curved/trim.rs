//! UV-space trim loops for curved panels.

use std::f64::consts::TAU;

use crate::boolean::poly2d::geom::Point2;
use crate::tolerance::Tol;

use super::CurvedError;

/// One directed edge of a UV-space trim loop.
///
/// Coordinates are in the surface parameter domain. For a cylinder panel,
/// `x = theta` in radians and `y = z` in metres.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum TrimEdge2d {
    /// A straight trim segment.
    Line {
        /// Segment start `[u, v]`.
        start: [f64; 2],
        /// Segment end `[u, v]`.
        end: [f64; 2],
    },
    /// A circular arc in UV space.
    Arc {
        /// Arc centre `[u, v]`.
        center: [f64; 2],
        /// Arc radius in UV units.
        radius: f64,
        /// Start angle in radians.
        start_angle: f64,
        /// End angle in radians.
        end_angle: f64,
    },
}

impl TrimEdge2d {
    /// Construct a straight trim segment.
    #[inline]
    pub fn line(start: [f64; 2], end: [f64; 2]) -> Self {
        Self::Line { start, end }
    }

    /// Start point.
    #[inline]
    pub fn start(self) -> [f64; 2] {
        match self {
            TrimEdge2d::Line { start, .. } => start,
            TrimEdge2d::Arc {
                center,
                radius,
                start_angle,
                ..
            } => [
                center[0] + radius * start_angle.cos(),
                center[1] + radius * start_angle.sin(),
            ],
        }
    }

    /// End point.
    #[inline]
    pub fn end(self) -> [f64; 2] {
        match self {
            TrimEdge2d::Line { end, .. } => end,
            TrimEdge2d::Arc {
                center,
                radius,
                end_angle,
                ..
            } => [
                center[0] + radius * end_angle.cos(),
                center[1] + radius * end_angle.sin(),
            ],
        }
    }

    pub(crate) fn validate(self) -> Result<(), CurvedError> {
        match self {
            TrimEdge2d::Line { start, end } => {
                if finite_point(start) && finite_point(end) {
                    Ok(())
                } else {
                    Err(CurvedError::InvalidTrimEdge)
                }
            }
            TrimEdge2d::Arc {
                center,
                radius,
                start_angle,
                end_angle,
            } => {
                if finite_point(center)
                    && radius.is_finite()
                    && radius > 0.0
                    && start_angle.is_finite()
                    && end_angle.is_finite()
                    && (end_angle - start_angle).abs() > 0.0
                    && (end_angle - start_angle).abs() <= TAU + 1e-12_f64
                {
                    Ok(())
                } else {
                    Err(CurvedError::InvalidTrimEdge)
                }
            }
        }
    }

    fn sample_points(self, chord_tolerance: f64) -> Vec<[f64; 2]> {
        match self {
            TrimEdge2d::Line { start, .. } => vec![start],
            TrimEdge2d::Arc {
                center,
                radius,
                start_angle,
                end_angle,
            } => {
                let sweep = end_angle - start_angle;
                let n = arc_segment_count(radius, sweep.abs(), chord_tolerance).max(1);
                (0..n)
                    .map(|i| {
                        let t = (i as f64) / (n as f64);
                        let a = start_angle + sweep * t;
                        [center[0] + radius * a.cos(), center[1] + radius * a.sin()]
                    })
                    .collect()
            }
        }
    }
}

/// A closed UV-space boundary loop.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TrimLoop2d {
    /// Directed edges in loop order.
    pub edges: Vec<TrimEdge2d>,
}

impl TrimLoop2d {
    /// Construct a trim loop after validating closure and non-degenerate area.
    pub fn new(edges: Vec<TrimEdge2d>, tol: &Tol) -> Result<Self, CurvedError> {
        let loop_ = Self { edges };
        loop_.validate(tol)?;
        Ok(loop_)
    }

    /// Construct an axis-aligned rectangular loop in CCW order.
    pub fn rectangle(
        u_min: f64,
        u_max: f64,
        v_min: f64,
        v_max: f64,
        tol: &Tol,
    ) -> Result<Self, CurvedError> {
        if !u_min.is_finite() || !u_max.is_finite() || u_max <= u_min {
            return Err(CurvedError::InvalidRange {
                name: "u",
                min: u_min,
                max: u_max,
            });
        }
        if !v_min.is_finite() || !v_max.is_finite() || v_max <= v_min {
            return Err(CurvedError::InvalidRange {
                name: "v",
                min: v_min,
                max: v_max,
            });
        }
        Self::from_points(
            &[
                [u_min, v_min],
                [u_max, v_min],
                [u_max, v_max],
                [u_min, v_max],
            ],
            tol,
        )
    }

    /// Construct a circular trim loop in UV space.
    pub fn circle(center: [f64; 2], radius: f64, tol: &Tol) -> Result<Self, CurvedError> {
        Self::new(
            vec![TrimEdge2d::Arc {
                center,
                radius,
                start_angle: 0.0,
                end_angle: TAU,
            }],
            tol,
        )
    }

    /// Construct a straight-edged loop from vertices.
    ///
    /// The closing vertex must be omitted; the loop closes implicitly.
    pub fn from_points(points: &[[f64; 2]], tol: &Tol) -> Result<Self, CurvedError> {
        let n = points.len();
        if n < 3 {
            return Err(CurvedError::DegenerateLoop { area: 0.0 });
        }
        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            edges.push(TrimEdge2d::line(points[i], points[(i + 1) % n]));
        }
        Self::new(edges, tol)
    }

    /// The loop vertices, one per edge.
    pub fn vertices(&self) -> Vec<[f64; 2]> {
        self.edges.iter().map(|e| e.start()).collect()
    }

    /// Signed UV-space area. Positive means CCW.
    pub fn signed_area(&self) -> f64 {
        let mut acc = 0.0_f64;
        for edge in &self.edges {
            let a = edge.start();
            let b = edge.end();
            acc += a[0] * b[1] - b[0] * a[1];
        }
        let mut area = 0.5 * acc;
        for edge in &self.edges {
            if let TrimEdge2d::Arc {
                radius,
                start_angle,
                end_angle,
                ..
            } = *edge
            {
                let sweep = end_angle - start_angle;
                area += 0.5 * radius * radius * (sweep - sweep.sin());
            }
        }
        area
    }

    /// Absolute UV-space area.
    #[inline]
    pub fn area(&self) -> f64 {
        self.signed_area().abs()
    }

    /// `true` if this loop contains an arc edge.
    pub fn has_arc(&self) -> bool {
        self.edges
            .iter()
            .any(|e| matches!(e, TrimEdge2d::Arc { .. }))
    }

    /// Return a reversed copy.
    pub fn reversed(&self) -> Self {
        let mut edges = Vec::with_capacity(self.edges.len());
        for edge in self.edges.iter().rev() {
            edges.push(match *edge {
                TrimEdge2d::Line { start, end } => TrimEdge2d::Line {
                    start: end,
                    end: start,
                },
                TrimEdge2d::Arc {
                    center,
                    radius,
                    start_angle,
                    end_angle,
                } => TrimEdge2d::Arc {
                    center,
                    radius,
                    start_angle: end_angle,
                    end_angle: start_angle,
                },
            });
        }
        Self { edges }
    }

    pub(crate) fn validate(&self, tol: &Tol) -> Result<(), CurvedError> {
        if self.edges.is_empty() {
            return Err(CurvedError::EmptyLoop);
        }
        for edge in &self.edges {
            edge.validate()?;
        }
        for i in 0..self.edges.len() {
            let a = self.edges[i].end();
            let b = self.edges[(i + 1) % self.edges.len()].start();
            if !point_coincident(a, b, tol) {
                return Err(CurvedError::OpenLoop);
            }
        }
        let area = self.signed_area();
        if !area.is_finite() || area.abs() <= tol.length * tol.length {
            return Err(CurvedError::DegenerateLoop { area });
        }
        Ok(())
    }

    pub(crate) fn contains_point(&self, p: [f64; 2], tol: &Tol) -> bool {
        let q = Point2::new(p[0], p[1]);
        let mut inside = false;
        for edge in &self.edges {
            match *edge {
                TrimEdge2d::Line { start, end } => {
                    let a = Point2::new(start[0], start[1]);
                    let b = Point2::new(end[0], end[1]);
                    if point_on_segment(q, a, b, tol) {
                        return true;
                    }
                    let crosses = (a.y > q.y) != (b.y > q.y);
                    if crosses {
                        let x = a.x + (q.y - a.y) * (b.x - a.x) / (b.y - a.y);
                        if x >= q.x - tol.length {
                            inside = !inside;
                        }
                    }
                }
                TrimEdge2d::Arc {
                    center,
                    radius,
                    start_angle,
                    end_angle,
                } => {
                    if point_on_arc(q, center, radius, start_angle, end_angle, tol) {
                        return true;
                    }
                    for x in ray_arc_crossings(q, center, radius, start_angle, end_angle, tol) {
                        if x >= q.x - tol.length {
                            inside = !inside;
                        }
                    }
                }
            }
        }
        inside
    }

    pub(crate) fn bounds(&self) -> ([f64; 2], [f64; 2]) {
        let mut min = [f64::INFINITY, f64::INFINITY];
        let mut max = [f64::NEG_INFINITY, f64::NEG_INFINITY];
        for p in self.bound_points() {
            min[0] = min[0].min(p[0]);
            min[1] = min[1].min(p[1]);
            max[0] = max[0].max(p[0]);
            max[1] = max[1].max(p[1]);
        }
        (min, max)
    }

    pub(crate) fn sample_points(&self, chord_tolerance: f64) -> Vec<[f64; 2]> {
        self.edges
            .iter()
            .flat_map(|e| e.sample_points(chord_tolerance))
            .collect()
    }

    fn bound_points(&self) -> Vec<[f64; 2]> {
        let mut points = self.vertices();
        for edge in &self.edges {
            if let TrimEdge2d::Arc {
                center,
                radius,
                start_angle,
                end_angle,
            } = *edge
            {
                for a in [0.0_f64, TAU / 4.0, TAU / 2.0, 3.0 * TAU / 4.0] {
                    if angle_on_sweep(a, start_angle, end_angle, 0.0) {
                        points.push([center[0] + radius * a.cos(), center[1] + radius * a.sin()]);
                    }
                }
            }
        }
        points
    }
}

pub(crate) fn point_coincident(a: [f64; 2], b: [f64; 2], tol: &Tol) -> bool {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy <= tol.length * tol.length
}

fn point_on_segment(p: Point2, a: Point2, b: Point2, tol: &Tol) -> bool {
    let ab = a.to(b);
    let ap = a.to(p);
    let cross = ab.cross(ap).abs();
    if cross > tol.length * ab.len().max(1.0) {
        return false;
    }
    let dot = ap.dot(ab);
    dot >= -tol.length && dot <= ab.len_sq() + tol.length
}

fn finite_point(p: [f64; 2]) -> bool {
    p[0].is_finite() && p[1].is_finite()
}

fn point_on_arc(
    p: Point2,
    center: [f64; 2],
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    tol: &Tol,
) -> bool {
    let dx = p.x - center[0];
    let dy = p.y - center[1];
    let d = (dx * dx + dy * dy).sqrt();
    if (d - radius).abs() > tol.length {
        return false;
    }
    let angle = dy.atan2(dx);
    angle_on_sweep(angle, start_angle, end_angle, tol.length / radius)
}

fn ray_arc_crossings(
    p: Point2,
    center: [f64; 2],
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    tol: &Tol,
) -> Vec<f64> {
    let dy = p.y - center[1];
    if dy.abs() >= radius - tol.length {
        return Vec::new();
    }
    let dx = (radius * radius - dy * dy).sqrt();
    let mut xs = Vec::new();
    for x in [center[0] - dx, center[0] + dx] {
        let angle = (p.y - center[1]).atan2(x - center[0]);
        if angle_on_sweep(angle, start_angle, end_angle, tol.length / radius) {
            xs.push(x);
        }
    }
    xs
}

fn angle_on_sweep(angle: f64, start: f64, end: f64, slack: f64) -> bool {
    let sweep = end - start;
    if sweep >= 0.0 {
        let off = (angle - start).rem_euclid(TAU);
        off <= sweep.abs() + slack || (sweep.abs() >= TAU - slack && off <= TAU)
    } else {
        let off = (start - angle).rem_euclid(TAU);
        off <= sweep.abs() + slack || (sweep.abs() >= TAU - slack && off <= TAU)
    }
}

fn arc_segment_count(radius: f64, sweep: f64, chord_tolerance: f64) -> usize {
    if chord_tolerance <= 0.0 || chord_tolerance >= radius {
        return ((sweep / (TAU / 32.0)).ceil() as usize).max(8);
    }
    let step = 2.0_f64 * (1.0_f64 - chord_tolerance / radius).acos();
    ((sweep / step).ceil() as usize).max(8)
}
