//! Conical panels with UV-space holes.

use std::f64::consts::{FRAC_PI_2, TAU};

use crate::math::{Point3, Unit3};
use crate::primitives::plane_basis;
use crate::tolerance::Tol;

use super::domain::{angular_step, parameter_values, validate_holes_in_rect, validate_range};
use super::mesh::SurfaceMeshBuilder;
use super::{CurvedError, SurfaceMesh, TrimLoop2d};

/// A finite frustum patch of a right circular cone with holes in `(theta, h)`.
///
/// `h` is distance from the apex along [`axis`](Self::axis), in metres.
/// The cone radius at height `h` is `h * tan(half_angle)`.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ConePanel {
    /// Cone apex in metres.
    pub apex: Point3,
    /// Cone axis, pointing from the apex toward increasing height.
    pub axis: Unit3,
    /// Cone half-angle in radians.
    pub half_angle: f64,
    /// Minimum angular parameter in radians.
    pub theta_min: f64,
    /// Maximum angular parameter in radians.
    pub theta_max: f64,
    /// Minimum height from apex in metres.
    pub height_min: f64,
    /// Maximum height from apex in metres.
    pub height_max: f64,
    /// Hole loops in the `(theta, h)` domain.
    pub holes: Vec<TrimLoop2d>,
}

/// Parameters for constructing a [`ConePanel`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ConePanelSpec {
    /// Cone apex in metres.
    pub apex: Point3,
    /// Cone axis, pointing from the apex toward increasing height.
    pub axis: Unit3,
    /// Cone half-angle in radians.
    pub half_angle: f64,
    /// Minimum angular parameter in radians.
    pub theta_min: f64,
    /// Maximum angular parameter in radians.
    pub theta_max: f64,
    /// Minimum height from apex in metres.
    pub height_min: f64,
    /// Maximum height from apex in metres.
    pub height_max: f64,
}

impl ConePanel {
    /// Construct a conical frustum panel with UV-space holes.
    ///
    /// This phase rejects domains touching the apex because longitude
    /// collapses there.
    pub fn new(
        spec: ConePanelSpec,
        holes: Vec<TrimLoop2d>,
        tol: &Tol,
    ) -> Result<Self, CurvedError> {
        let ConePanelSpec {
            apex,
            axis,
            half_angle,
            theta_min,
            theta_max,
            height_min,
            height_max,
        } = spec;
        if !half_angle.is_finite() || half_angle <= 0.0 || half_angle >= FRAC_PI_2 - tol.length {
            return Err(CurvedError::InvalidConeAngle { value: half_angle });
        }
        validate_range("theta", theta_min, theta_max)?;
        validate_range("height", height_min, height_max)?;
        if height_min <= tol.length {
            return Err(CurvedError::ApexCrossing);
        }
        if theta_max - theta_min > TAU + tol.length {
            return Err(CurvedError::SeamCrossing);
        }
        validate_holes_in_rect(&holes, theta_min, theta_max, height_min, height_max, tol)?;
        for hole in &holes {
            let (min, max) = hole.bounds();
            if min[0] <= theta_min + tol.length && max[0] >= theta_max - tol.length {
                return Err(CurvedError::SeamCrossing);
            }
        }
        Ok(Self {
            apex,
            axis,
            half_angle,
            theta_min,
            theta_max,
            height_min,
            height_max,
            holes,
        })
    }

    /// Radius at height `h` along the cone axis.
    pub fn radius_at_height(&self, h: f64) -> f64 {
        h * self.half_angle.tan()
    }

    /// Exact surface area of the untrimmed frustum rectangle.
    pub fn untrimmed_surface_area(&self) -> f64 {
        0.5_f64 * self.half_angle.tan() / self.half_angle.cos()
            * (self.height_max * self.height_max - self.height_min * self.height_min)
            * (self.theta_max - self.theta_min)
    }

    /// Map `(theta, h)` to a world-space point on the cone.
    pub fn point_at(&self, theta: f64, h: f64) -> Point3 {
        let (u, v) = plane_basis(self.axis);
        let radius = self.radius_at_height(h);
        self.apex + self.axis.as_vec() * h + (u * theta.cos() + v * theta.sin()) * radius
    }

    fn material_contains(&self, theta: f64, height: f64, tol: &Tol) -> bool {
        let p = [theta, height];
        !self.holes.iter().any(|h| h.contains_point(p, tol))
    }
}

/// Tessellation controls for conical panels.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct ConePanelOptions {
    /// Maximum chord error along the angular direction, in metres.
    pub chord_tolerance: f64,
}

impl Default for ConePanelOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: 1e-3_f64,
        }
    }
}

impl ConePanelOptions {
    /// Options with the given chord tolerance, in metres.
    pub fn with_chord_tolerance(chord_tolerance: f64) -> Self {
        Self { chord_tolerance }
    }
}

/// Tessellate a [`ConePanel`] into a display [`SurfaceMesh`].
pub fn tessellate_cone_panel(
    panel: &ConePanel,
    opts: &ConePanelOptions,
    tol: &Tol,
) -> Result<SurfaceMesh, CurvedError> {
    if opts.chord_tolerance <= 0.0 || !opts.chord_tolerance.is_finite() {
        return Err(CurvedError::NonPositiveChordTolerance {
            value: opts.chord_tolerance,
        });
    }

    let extras = panel
        .holes
        .iter()
        .flat_map(|h| h.sample_points(opts.chord_tolerance))
        .collect::<Vec<_>>();
    let theta_values = parameter_values(
        panel.theta_min,
        panel.theta_max,
        angular_step(
            panel.radius_at_height(panel.height_max),
            opts.chord_tolerance,
        ),
        extras.iter().map(|p| p[0]),
        tol,
    );
    let height_step = (opts.chord_tolerance / panel.half_angle.cos()).max(tol.length);
    let height_values = parameter_values(
        panel.height_min,
        panel.height_max,
        height_step,
        extras.iter().map(|p| p[1]),
        tol,
    );

    let mut builder = SurfaceMeshBuilder::default();
    let mut grid = vec![vec![0_u32; height_values.len()]; theta_values.len()];
    for (i, theta) in theta_values.iter().enumerate() {
        for (j, height) in height_values.iter().enumerate() {
            grid[i][j] = builder.vertex(panel.point_at(*theta, *height));
        }
    }

    for i in 0..theta_values.len() - 1 {
        for j in 0..height_values.len() - 1 {
            let theta_mid = 0.5 * (theta_values[i] + theta_values[i + 1]);
            let height_mid = 0.5 * (height_values[j] + height_values[j + 1]);
            if !panel.material_contains(theta_mid, height_mid, tol) {
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
