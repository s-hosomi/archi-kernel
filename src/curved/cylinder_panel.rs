//! Cylindrical panels with UV-space holes.

use std::f64::consts::TAU;

use crate::math::{Point3, Vec3};
use crate::primitives::{plane_basis, Cylinder};
use crate::tolerance::Tol;

use super::domain::{angular_step, parameter_values, validate_holes_in_rect, validate_range};
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
    /// Loops must be contained in the unwrapped rectangle
    /// `[theta_min, theta_max] × [z_min, z_max]` and must not cross or overlap.
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
        self.point_at_radius(self.cylinder.radius(), theta, z)
    }

    fn point_at_radius(&self, radius: f64, theta: f64, z: f64) -> Point3 {
        let axis = self.cylinder.axis();
        let origin = axis.point_at(z);
        let (u, v) = plane_basis(axis.dir());
        origin + u * (radius * theta.cos()) + v * (radius * theta.sin())
    }

    fn validate_holes(&self, tol: &Tol) -> Result<(), CurvedError> {
        validate_holes_in_rect(
            &self.holes,
            self.theta_min,
            self.theta_max,
            self.z_min,
            self.z_max,
            tol,
        )?;
        for hole in &self.holes {
            let (min, max) = hole.bounds();
            if min[0] <= self.theta_min + tol.length && max[0] >= self.theta_max - tol.length {
                return Err(CurvedError::SeamCrossing);
            }
        }
        Ok(())
    }

    fn material_contains(&self, theta: f64, z: f64, tol: &Tol) -> bool {
        let p = [theta, z];
        !self.holes.iter().any(|h| h.contains_point(p, tol))
    }
}

/// A cylindrical panel with finite radial thickness.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ThickCylinderPanel {
    /// The mid-surface panel that defines the UV domain and holes.
    pub mid: CylinderPanel,
    /// Radial thickness in metres, centred on [`mid`](Self::mid).
    pub thickness: f64,
}

impl ThickCylinderPanel {
    /// Construct a thick cylindrical panel by offsetting the mid-surface radius
    /// by `± thickness / 2`.
    pub fn new(mid: CylinderPanel, thickness: f64) -> Result<Self, CurvedError> {
        if thickness <= 0.0 || !thickness.is_finite() {
            return Err(CurvedError::NonPositiveThickness { value: thickness });
        }
        let inner = mid.cylinder.radius() - 0.5 * thickness;
        if inner <= 0.0 {
            return Err(CurvedError::NonPositiveInnerRadius { radius: inner });
        }
        Ok(Self { mid, thickness })
    }

    /// Exact volume represented by the thick panel.
    pub fn volume(&self) -> f64 {
        self.mid.cylinder.radius() * self.thickness * self.mid.uv_area()
    }

    fn inner_radius(&self) -> f64 {
        self.mid.cylinder.radius() - 0.5 * self.thickness
    }

    fn outer_radius(&self) -> f64 {
        self.mid.cylinder.radius() + 0.5 * self.thickness
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
        panel.holes.iter().flat_map(|h| {
            h.sample_points(opts.chord_tolerance)
                .into_iter()
                .map(|p| p[0])
                .collect::<Vec<f64>>()
        }),
        tol,
    );
    let z_values = parameter_values(
        panel.z_min,
        panel.z_max,
        opts.chord_tolerance.max(tol.length),
        panel.holes.iter().flat_map(|h| {
            h.sample_points(opts.chord_tolerance)
                .into_iter()
                .map(|p| p[1])
                .collect::<Vec<f64>>()
        }),
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

/// Tessellate a thick cylindrical panel into a closed display mesh.
pub fn tessellate_thick_cylinder_panel(
    panel: &ThickCylinderPanel,
    opts: &CylinderPanelOptions,
    tol: &Tol,
) -> Result<SurfaceMesh, CurvedError> {
    if opts.chord_tolerance <= 0.0 || !opts.chord_tolerance.is_finite() {
        return Err(CurvedError::NonPositiveChordTolerance {
            value: opts.chord_tolerance,
        });
    }

    let mid = &panel.mid;
    if mid.holes.iter().any(TrimLoop2d::has_arc) {
        return Err(CurvedError::UnsupportedArcTrim);
    }
    let theta_values = parameter_values(
        mid.theta_min,
        mid.theta_max,
        angular_step(mid.cylinder.radius(), opts.chord_tolerance),
        mid.holes.iter().flat_map(|h| {
            h.sample_points(opts.chord_tolerance)
                .into_iter()
                .map(|p| p[0])
                .collect::<Vec<f64>>()
        }),
        tol,
    );
    let z_values = parameter_values(
        mid.z_min,
        mid.z_max,
        opts.chord_tolerance.max(tol.length),
        mid.holes.iter().flat_map(|h| {
            h.sample_points(opts.chord_tolerance)
                .into_iter()
                .map(|p| p[1])
                .collect::<Vec<f64>>()
        }),
        tol,
    );
    let nt = theta_values.len() - 1;
    let nz = z_values.len() - 1;
    let mut material = vec![vec![false; nz]; nt];
    for i in 0..nt {
        for j in 0..nz {
            let theta_mid = 0.5 * (theta_values[i] + theta_values[i + 1]);
            let z_mid = 0.5 * (z_values[j] + z_values[j + 1]);
            material[i][j] = mid.material_contains(theta_mid, z_mid, tol);
        }
    }

    let mut builder = SurfaceMeshBuilder::default();
    let ri = panel.inner_radius();
    let ro = panel.outer_radius();
    for i in 0..nt {
        for j in 0..nz {
            if !material[i][j] {
                continue;
            }
            let t0 = theta_values[i];
            let t1 = theta_values[i + 1];
            let z0 = z_values[j];
            let z1 = z_values[j + 1];
            emit_outer_cell(&mut builder, mid, ro, t0, t1, z0, z1);
            emit_inner_cell(&mut builder, mid, ri, t0, t1, z0, z1);

            if i == 0 || !material[i - 1][j] {
                emit_radial_side(&mut builder, mid, [ri, ro], t0, [z0, z1], false);
            }
            if i + 1 == nt || !material[i + 1][j] {
                emit_radial_side(&mut builder, mid, [ri, ro], t1, [z0, z1], true);
            }
            if j == 0 || !material[i][j - 1] {
                emit_z_side(&mut builder, mid, [ri, ro], [t0, t1], z0, false);
            }
            if j + 1 == nz || !material[i][j + 1] {
                emit_z_side(&mut builder, mid, [ri, ro], [t0, t1], z1, true);
            }
        }
    }

    Ok(builder.finish())
}

fn emit_outer_cell(
    b: &mut SurfaceMeshBuilder,
    panel: &CylinderPanel,
    radius: f64,
    t0: f64,
    t1: f64,
    z0: f64,
    z1: f64,
) {
    let a = b.vertex(panel.point_at_radius(radius, t0, z0));
    let bb = b.vertex(panel.point_at_radius(radius, t1, z0));
    let c = b.vertex(panel.point_at_radius(radius, t1, z1));
    let d = b.vertex(panel.point_at_radius(radius, t0, z1));
    b.triangle(a, bb, c);
    b.triangle(a, c, d);
}

fn emit_inner_cell(
    b: &mut SurfaceMeshBuilder,
    panel: &CylinderPanel,
    radius: f64,
    t0: f64,
    t1: f64,
    z0: f64,
    z1: f64,
) {
    let a = b.vertex(panel.point_at_radius(radius, t0, z0));
    let bb = b.vertex(panel.point_at_radius(radius, t1, z0));
    let c = b.vertex(panel.point_at_radius(radius, t1, z1));
    let d = b.vertex(panel.point_at_radius(radius, t0, z1));
    b.triangle(a, c, bb);
    b.triangle(a, d, c);
}

fn emit_radial_side(
    b: &mut SurfaceMeshBuilder,
    panel: &CylinderPanel,
    radii: [f64; 2],
    theta: f64,
    z: [f64; 2],
    theta_max_side: bool,
) {
    let [ri, ro] = radii;
    let [z0, z1] = z;
    let a = b.vertex(panel.point_at_radius(ri, theta, z0));
    let bb = b.vertex(panel.point_at_radius(ro, theta, z0));
    let c = b.vertex(panel.point_at_radius(ro, theta, z1));
    let d = b.vertex(panel.point_at_radius(ri, theta, z1));
    if theta_max_side {
        b.triangle(a, bb, c);
        b.triangle(a, c, d);
    } else {
        b.triangle(a, c, bb);
        b.triangle(a, d, c);
    }
}

fn emit_z_side(
    b: &mut SurfaceMeshBuilder,
    panel: &CylinderPanel,
    radii: [f64; 2],
    theta: [f64; 2],
    z: f64,
    z_max_side: bool,
) {
    let [ri, ro] = radii;
    let [t0, t1] = theta;
    let a = b.vertex(panel.point_at_radius(ri, t0, z));
    let bb = b.vertex(panel.point_at_radius(ro, t0, z));
    let c = b.vertex(panel.point_at_radius(ro, t1, z));
    let d = b.vertex(panel.point_at_radius(ri, t1, z));
    if z_max_side {
        b.triangle(a, bb, c);
        b.triangle(a, c, d);
    } else {
        b.triangle(a, c, bb);
        b.triangle(a, d, c);
    }
}

#[allow(dead_code)]
fn _normal_at(axis: Vec3, theta: f64) -> Vec3 {
    let Some(unit) = axis.try_unit() else {
        return Vec3::ZERO;
    };
    let (u, v) = plane_basis(unit);
    u * theta.cos() + v * theta.sin()
}
