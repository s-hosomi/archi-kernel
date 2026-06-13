//! Cylindrical panels with UV-space holes.

use std::f64::consts::{PI, TAU};

use crate::math::{Point3, Vec3};
use crate::primitives::{plane_basis, Cylinder};
use crate::tolerance::Tol;

use super::mesh::SurfaceMeshBuilder;
use super::{CurvedError, SurfaceMesh, TrimLoop2d};

/// A finite patch of a circular cylinder with holes defined in `(theta, z)`.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CylinderPanel {
    /// The supporting infinite cylinder.
    pub cylinder: Cylinder,
    /// Minimum unwrapped angular parameter in radians.
    pub theta_min: f64,
    /// Maximum unwrapped angular parameter in radians.
    pub theta_max: f64,
    /// Minimum axis parameter in metres.
    pub z_min: f64,
    /// Maximum axis parameter in metres.
    pub z_max: f64,
    /// Hole loops in the `(theta, z)` domain.
    pub holes: Vec<TrimLoop2d>,
}

impl CylinderPanel {
    /// Construct a cylindrical panel with UV-space holes.
    ///
    /// This phase accepts straight-edged trim loops only. Loops must be contained
    /// in the unwrapped rectangle `[theta_min, theta_max] × [z_min, z_max]` and
    /// must not cross or overlap.
    pub fn new(
        cylinder: Cylinder,
        theta_min: f64,
        theta_max: f64,
        z_min: f64,
        z_max: f64,
        holes: Vec<TrimLoop2d>,
        tol: &Tol,
    ) -> Result<Self, CurvedError> {
        validate_range("theta", theta_min, theta_max)?;
        validate_range("z", z_min, z_max)?;
        if theta_max - theta_min > TAU + tol.length {
            return Err(CurvedError::SeamCrossing);
        }

        let panel = Self {
            cylinder,
            theta_min,
            theta_max,
            z_min,
            z_max,
            holes,
        };
        panel.validate_holes(tol)?;
        Ok(panel)
    }

    /// UV-space area of the material region, in radians-metres.
    pub fn uv_area(&self) -> f64 {
        let outer = (self.theta_max - self.theta_min) * (self.z_max - self.z_min);
        let holes: f64 = self.holes.iter().map(TrimLoop2d::area).sum();
        outer - holes
    }

    /// Exact cylindrical surface area represented by this panel.
    pub fn surface_area(&self) -> f64 {
        self.cylinder.radius() * self.uv_area()
    }

    /// Map `(theta, z)` to a world-space point on the cylinder.
    pub fn point_at(&self, theta: f64, z: f64) -> Point3 {
        let axis = self.cylinder.axis();
        let origin = axis.point_at(z);
        let (u, v) = plane_basis(axis.dir());
        origin
            + u * (self.cylinder.radius() * theta.cos())
            + v * (self.cylinder.radius() * theta.sin())
    }

    fn validate_holes(&self, tol: &Tol) -> Result<(), CurvedError> {
        for hole in &self.holes {
            hole.validate(tol)?;
            let (min, max) = hole.bounds();
            if min[0] < self.theta_min - tol.length
                || max[0] > self.theta_max + tol.length
                || min[1] < self.z_min - tol.length
                || max[1] > self.z_max + tol.length
            {
                return Err(CurvedError::HoleOutsidePanel);
            }
            if min[0] <= self.theta_min + tol.length && max[0] >= self.theta_max - tol.length {
                return Err(CurvedError::SeamCrossing);
            }
        }
        for i in 0..self.holes.len() {
            for j in (i + 1)..self.holes.len() {
                if loops_overlap(&self.holes[i], &self.holes[j], tol) {
                    return Err(CurvedError::HoleOverlap);
                }
            }
        }
        Ok(())
    }

    fn material_contains(&self, theta: f64, z: f64, tol: &Tol) -> bool {
        let p = [theta, z];
        !self.holes.iter().any(|h| h.contains_point(p, tol))
    }
}

/// Tessellation controls for cylindrical panels.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct CylinderPanelOptions {
    /// Maximum chord error along the angular direction, in metres.
    pub chord_tolerance: f64,
}

impl Default for CylinderPanelOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: 1e-3_f64,
        }
    }
}

impl CylinderPanelOptions {
    /// Options with the given chord tolerance, in metres.
    pub fn with_chord_tolerance(chord_tolerance: f64) -> Self {
        Self { chord_tolerance }
    }
}

/// Tessellate a [`CylinderPanel`] into a display [`SurfaceMesh`].
///
/// The output is a surface mesh only; it is not a closed solid and is not used by
/// the B-rep CSG evaluator.
pub fn tessellate_cylinder_panel(
    panel: &CylinderPanel,
    opts: &CylinderPanelOptions,
    tol: &Tol,
) -> Result<SurfaceMesh, CurvedError> {
    if opts.chord_tolerance <= 0.0 || !opts.chord_tolerance.is_finite() {
        return Err(CurvedError::NonPositiveChordTolerance {
            value: opts.chord_tolerance,
        });
    }

    let theta_values = parameter_values(
        panel.theta_min,
        panel.theta_max,
        angular_step(panel.cylinder.radius(), opts.chord_tolerance),
        panel
            .holes
            .iter()
            .flat_map(|h| h.vertices().into_iter().map(|p| p[0]).collect::<Vec<f64>>()),
        tol,
    );
    let z_values = parameter_values(
        panel.z_min,
        panel.z_max,
        opts.chord_tolerance.max(tol.length),
        panel
            .holes
            .iter()
            .flat_map(|h| h.vertices().into_iter().map(|p| p[1]).collect::<Vec<f64>>()),
        tol,
    );

    let mut builder = SurfaceMeshBuilder::default();
    let mut grid = vec![vec![0_u32; z_values.len()]; theta_values.len()];
    for (i, theta) in theta_values.iter().enumerate() {
        for (j, z) in z_values.iter().enumerate() {
            grid[i][j] = builder.vertex(panel.point_at(*theta, *z));
        }
    }

    for i in 0..theta_values.len() - 1 {
        for j in 0..z_values.len() - 1 {
            let theta_mid = 0.5 * (theta_values[i] + theta_values[i + 1]);
            let z_mid = 0.5 * (z_values[j] + z_values[j + 1]);
            if !panel.material_contains(theta_mid, z_mid, tol) {
                continue;
            }
            let a = grid[i][j];
            let b = grid[i + 1][j];
            let c = grid[i + 1][j + 1];
            let d = grid[i][j + 1];
            builder.triangle(a, b, c);
            builder.triangle(a, c, d);
        }
    }

    Ok(builder.finish())
}

fn validate_range(name: &'static str, min: f64, max: f64) -> Result<(), CurvedError> {
    if !min.is_finite() || !max.is_finite() || max <= min {
        return Err(CurvedError::InvalidRange { name, min, max });
    }
    Ok(())
}

fn angular_step(radius: f64, chord_tol: f64) -> f64 {
    if chord_tol >= radius {
        PI / 8.0
    } else {
        (2.0_f64 * (1.0_f64 - chord_tol / radius).acos()).clamp(PI / 180.0, PI / 8.0)
    }
}

fn parameter_values(
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

fn loops_overlap(a: &TrimLoop2d, b: &TrimLoop2d, tol: &Tol) -> bool {
    if !bounds_overlap(a, b, tol) {
        return false;
    }
    for ea in &a.edges {
        let super::TrimEdge2d::Line { start: a0, end: a1 } = *ea else {
            return true;
        };
        for eb in &b.edges {
            let super::TrimEdge2d::Line { start: b0, end: b1 } = *eb else {
                return true;
            };
            if segments_intersect(a0, a1, b0, b1, tol) {
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

fn segments_intersect(a0: [f64; 2], a1: [f64; 2], b0: [f64; 2], b1: [f64; 2], tol: &Tol) -> bool {
    let o1 = orient(a0, a1, b0);
    let o2 = orient(a0, a1, b1);
    let o3 = orient(b0, b1, a0);
    let o4 = orient(b0, b1, a1);
    if o1.abs() <= tol.length && point_in_box(b0, a0, a1, tol) {
        return true;
    }
    if o2.abs() <= tol.length && point_in_box(b1, a0, a1, tol) {
        return true;
    }
    if o3.abs() <= tol.length && point_in_box(a0, b0, b1, tol) {
        return true;
    }
    if o4.abs() <= tol.length && point_in_box(a1, b0, b1, tol) {
        return true;
    }
    (o1 > 0.0) != (o2 > 0.0) && (o3 > 0.0) != (o4 > 0.0)
}

fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn point_in_box(p: [f64; 2], a: [f64; 2], b: [f64; 2], tol: &Tol) -> bool {
    p[0] >= a[0].min(b[0]) - tol.length
        && p[0] <= a[0].max(b[0]) + tol.length
        && p[1] >= a[1].min(b[1]) - tol.length
        && p[1] <= a[1].max(b[1]) + tol.length
}

#[allow(dead_code)]
fn _normal_at(axis: Vec3, theta: f64) -> Vec3 {
    let Some(unit) = axis.try_unit() else {
        return Vec3::ZERO;
    };
    let (u, v) = plane_basis(unit);
    u * theta.cos() + v * theta.sin()
}
