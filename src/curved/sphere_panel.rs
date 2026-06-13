//! Spherical panels with UV-space holes.

use std::f64::consts::{FRAC_PI_2, TAU};

use crate::math::{Point3, Unit3};
use crate::primitives::plane_basis;
use crate::tolerance::Tol;

use super::domain::{angular_step, parameter_values, validate_holes_in_rect, validate_range};
use super::mesh::SurfaceMeshBuilder;
use super::{CurvedError, SurfaceMesh, TrimLoop2d};

/// A finite latitude-longitude patch of a sphere with holes in `(theta, phi)`.
///
/// `theta` is longitude around [`pole`](Self::pole), in radians. `phi` is
/// latitude, also in radians, with `0` on the equator and `±π/2` at the poles.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpherePanel {
    /// Sphere centre in metres.
    pub center: Point3,
    /// Sphere radius in metres.
    pub radius: f64,
    /// Positive latitude direction.
    pub pole: Unit3,
    /// Minimum longitude in radians.
    pub theta_min: f64,
    /// Maximum longitude in radians.
    pub theta_max: f64,
    /// Minimum latitude in radians.
    pub phi_min: f64,
    /// Maximum latitude in radians.
    pub phi_max: f64,
    /// Hole loops in the `(theta, phi)` domain.
    pub holes: Vec<TrimLoop2d>,
}

/// Parameters for constructing a [`SpherePanel`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpherePanelSpec {
    /// Sphere centre in metres.
    pub center: Point3,
    /// Sphere radius in metres.
    pub radius: f64,
    /// Positive latitude direction.
    pub pole: Unit3,
    /// Minimum longitude in radians.
    pub theta_min: f64,
    /// Maximum longitude in radians.
    pub theta_max: f64,
    /// Minimum latitude in radians.
    pub phi_min: f64,
    /// Maximum latitude in radians.
    pub phi_max: f64,
}

impl SpherePanel {
    /// Construct a spherical panel with UV-space holes.
    ///
    /// This first phase rejects patches that include a pole because longitude
    /// collapses there. Use separate cap-specific primitives for exact polar
    /// domes in a later phase.
    pub fn new(
        spec: SpherePanelSpec,
        holes: Vec<TrimLoop2d>,
        tol: &Tol,
    ) -> Result<Self, CurvedError> {
        let SpherePanelSpec {
            center,
            radius,
            pole,
            theta_min,
            theta_max,
            phi_min,
            phi_max,
        } = spec;
        if radius <= 0.0 || !radius.is_finite() {
            return Err(CurvedError::NonPositiveRadius { value: radius });
        }
        validate_range("theta", theta_min, theta_max)?;
        validate_range("phi", phi_min, phi_max)?;
        if theta_max - theta_min > TAU + tol.length {
            return Err(CurvedError::SeamCrossing);
        }
        if phi_min <= -FRAC_PI_2 + tol.length || phi_max >= FRAC_PI_2 - tol.length {
            return Err(CurvedError::PoleCrossing);
        }
        validate_holes_in_rect(&holes, theta_min, theta_max, phi_min, phi_max, tol)?;
        for hole in &holes {
            let (min, max) = hole.bounds();
            if min[0] <= theta_min + tol.length && max[0] >= theta_max - tol.length {
                return Err(CurvedError::SeamCrossing);
            }
        }
        Ok(Self {
            center,
            radius,
            pole,
            theta_min,
            theta_max,
            phi_min,
            phi_max,
            holes,
        })
    }

    /// Exact surface area of the untrimmed latitude-longitude rectangle.
    pub fn untrimmed_surface_area(&self) -> f64 {
        self.radius
            * self.radius
            * (self.theta_max - self.theta_min)
            * (self.phi_max.sin() - self.phi_min.sin())
    }

    /// Map `(theta, phi)` to a world-space point on the sphere.
    pub fn point_at(&self, theta: f64, phi: f64) -> Point3 {
        self.point_at_radius(self.radius, theta, phi)
    }

    fn point_at_radius(&self, radius: f64, theta: f64, phi: f64) -> Point3 {
        let (u, v) = plane_basis(self.pole);
        let equator = u * theta.cos() + v * theta.sin();
        self.center + equator * (radius * phi.cos()) + self.pole.as_vec() * (radius * phi.sin())
    }

    fn material_contains(&self, theta: f64, phi: f64, tol: &Tol) -> bool {
        let p = [theta, phi];
        !self.holes.iter().any(|h| h.contains_point(p, tol))
    }
}

/// A spherical panel with radial thickness.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ThickSpherePanel {
    /// The mid-surface panel that defines the UV domain and holes.
    pub mid: SpherePanel,
    /// Radial thickness in metres, centred on [`mid`](Self::mid).
    pub thickness: f64,
}

impl ThickSpherePanel {
    /// Construct a thick spherical panel by offsetting the mid-surface radius
    /// by `± thickness / 2`.
    pub fn new(mid: SpherePanel, thickness: f64) -> Result<Self, CurvedError> {
        if thickness <= 0.0 || !thickness.is_finite() {
            return Err(CurvedError::NonPositiveThickness { value: thickness });
        }
        let inner = mid.radius - 0.5 * thickness;
        if inner <= 0.0 {
            return Err(CurvedError::NonPositiveInnerRadius { radius: inner });
        }
        Ok(Self { mid, thickness })
    }

    fn inner_radius(&self) -> f64 {
        self.mid.radius - 0.5 * self.thickness
    }

    fn outer_radius(&self) -> f64 {
        self.mid.radius + 0.5 * self.thickness
    }
}

/// Tessellation controls for spherical panels.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct SpherePanelOptions {
    /// Maximum chord error along angular directions, in metres.
    pub chord_tolerance: f64,
}

impl Default for SpherePanelOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: 1e-3_f64,
        }
    }
}

impl SpherePanelOptions {
    /// Options with the given chord tolerance, in metres.
    pub fn with_chord_tolerance(chord_tolerance: f64) -> Self {
        Self { chord_tolerance }
    }
}

/// Tessellate a [`SpherePanel`] into a display [`SurfaceMesh`].
pub fn tessellate_sphere_panel(
    panel: &SpherePanel,
    opts: &SpherePanelOptions,
    tol: &Tol,
) -> Result<SurfaceMesh, CurvedError> {
    validate_chord_tolerance(opts.chord_tolerance)?;
    let (theta_values, phi_values) = sphere_parameter_values(panel, opts.chord_tolerance, tol);

    let mut builder = SurfaceMeshBuilder::default();
    let mut grid = vec![vec![0_u32; phi_values.len()]; theta_values.len()];
    for (i, theta) in theta_values.iter().enumerate() {
        for (j, phi) in phi_values.iter().enumerate() {
            grid[i][j] = builder.vertex(panel.point_at(*theta, *phi));
        }
    }

    for i in 0..theta_values.len() - 1 {
        for j in 0..phi_values.len() - 1 {
            let theta_mid = 0.5 * (theta_values[i] + theta_values[i + 1]);
            let phi_mid = 0.5 * (phi_values[j] + phi_values[j + 1]);
            if !panel.material_contains(theta_mid, phi_mid, tol) {
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

/// Tessellate a thick spherical panel into a closed display mesh.
pub fn tessellate_thick_sphere_panel(
    panel: &ThickSpherePanel,
    opts: &SpherePanelOptions,
    tol: &Tol,
) -> Result<SurfaceMesh, CurvedError> {
    validate_chord_tolerance(opts.chord_tolerance)?;
    let mid = &panel.mid;
    if mid.holes.iter().any(TrimLoop2d::has_arc) {
        return Err(CurvedError::UnsupportedArcTrim);
    }
    let (theta_values, phi_values) = sphere_parameter_values(mid, opts.chord_tolerance, tol);
    let nt = theta_values.len() - 1;
    let np = phi_values.len() - 1;
    let mut material = vec![vec![false; np]; nt];
    for i in 0..nt {
        for j in 0..np {
            let theta_mid = 0.5 * (theta_values[i] + theta_values[i + 1]);
            let phi_mid = 0.5 * (phi_values[j] + phi_values[j + 1]);
            material[i][j] = mid.material_contains(theta_mid, phi_mid, tol);
        }
    }

    let mut builder = SurfaceMeshBuilder::default();
    let ri = panel.inner_radius();
    let ro = panel.outer_radius();
    for i in 0..nt {
        for j in 0..np {
            if !material[i][j] {
                continue;
            }
            let t0 = theta_values[i];
            let t1 = theta_values[i + 1];
            let p0 = phi_values[j];
            let p1 = phi_values[j + 1];
            emit_sphere_cell(&mut builder, mid, ro, [t0, t1], [p0, p1], false);
            emit_sphere_cell(&mut builder, mid, ri, [t0, t1], [p0, p1], true);

            if i == 0 || !material[i - 1][j] {
                emit_theta_side(&mut builder, mid, [ri, ro], t0, [p0, p1], false);
            }
            if i + 1 == nt || !material[i + 1][j] {
                emit_theta_side(&mut builder, mid, [ri, ro], t1, [p0, p1], true);
            }
            if j == 0 || !material[i][j - 1] {
                emit_phi_side(&mut builder, mid, [ri, ro], [t0, t1], p0, false);
            }
            if j + 1 == np || !material[i][j + 1] {
                emit_phi_side(&mut builder, mid, [ri, ro], [t0, t1], p1, true);
            }
        }
    }

    Ok(builder.finish())
}

fn validate_chord_tolerance(value: f64) -> Result<(), CurvedError> {
    if value <= 0.0 || !value.is_finite() {
        return Err(CurvedError::NonPositiveChordTolerance { value });
    }
    Ok(())
}

fn sphere_parameter_values(
    panel: &SpherePanel,
    chord_tolerance: f64,
    tol: &Tol,
) -> (Vec<f64>, Vec<f64>) {
    let extras = panel
        .holes
        .iter()
        .flat_map(|h| h.sample_points(chord_tolerance))
        .collect::<Vec<_>>();
    let theta_values = parameter_values(
        panel.theta_min,
        panel.theta_max,
        angular_step(panel.radius, chord_tolerance),
        extras.iter().map(|p| p[0]),
        tol,
    );
    let phi_values = parameter_values(
        panel.phi_min,
        panel.phi_max,
        angular_step(panel.radius, chord_tolerance),
        extras.iter().map(|p| p[1]),
        tol,
    );
    (theta_values, phi_values)
}

fn emit_sphere_cell(
    b: &mut SurfaceMeshBuilder,
    panel: &SpherePanel,
    radius: f64,
    theta: [f64; 2],
    phi: [f64; 2],
    inward: bool,
) {
    let [t0, t1] = theta;
    let [p0, p1] = phi;
    let a = b.vertex(panel.point_at_radius(radius, t0, p0));
    let bb = b.vertex(panel.point_at_radius(radius, t1, p0));
    let c = b.vertex(panel.point_at_radius(radius, t1, p1));
    let d = b.vertex(panel.point_at_radius(radius, t0, p1));
    if inward {
        b.triangle(a, c, bb);
        b.triangle(a, d, c);
    } else {
        b.triangle(a, bb, c);
        b.triangle(a, c, d);
    }
}

fn emit_theta_side(
    b: &mut SurfaceMeshBuilder,
    panel: &SpherePanel,
    radii: [f64; 2],
    theta: f64,
    phi: [f64; 2],
    theta_max_side: bool,
) {
    let [ri, ro] = radii;
    let [p0, p1] = phi;
    let a = b.vertex(panel.point_at_radius(ri, theta, p0));
    let bb = b.vertex(panel.point_at_radius(ro, theta, p0));
    let c = b.vertex(panel.point_at_radius(ro, theta, p1));
    let d = b.vertex(panel.point_at_radius(ri, theta, p1));
    if theta_max_side {
        b.triangle(a, bb, c);
        b.triangle(a, c, d);
    } else {
        b.triangle(a, c, bb);
        b.triangle(a, d, c);
    }
}

fn emit_phi_side(
    b: &mut SurfaceMeshBuilder,
    panel: &SpherePanel,
    radii: [f64; 2],
    theta: [f64; 2],
    phi: f64,
    phi_max_side: bool,
) {
    let [ri, ro] = radii;
    let [t0, t1] = theta;
    let a = b.vertex(panel.point_at_radius(ri, t0, phi));
    let bb = b.vertex(panel.point_at_radius(ro, t0, phi));
    let c = b.vertex(panel.point_at_radius(ro, t1, phi));
    let d = b.vertex(panel.point_at_radius(ri, t1, phi));
    if phi_max_side {
        b.triangle(a, bb, c);
        b.triangle(a, c, d);
    } else {
        b.triangle(a, c, bb);
        b.triangle(a, d, c);
    }
}
