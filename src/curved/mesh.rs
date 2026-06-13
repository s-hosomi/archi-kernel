//! Display mesh for trimmed curved panels.

use crate::math::Point3;

/// Indexed triangle mesh emitted by curved-panel tessellation.
///
/// Unlike [`crate::tess::Mesh`], this mesh is not tied to B-rep face arena IDs.
/// It is a display/export surface mesh for trimmed panels.
#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SurfaceMesh {
    /// Flat `xyz` vertex coordinates in metres.
    pub positions: Vec<f64>,
    /// Triangle corner indices, three per triangle.
    pub indices: Vec<u32>,
}

impl SurfaceMesh {
    /// Number of vertices.
    #[inline]
    pub fn vertex_count(&self) -> usize {
        self.positions.len() / 3
    }

    /// Number of triangles.
    #[inline]
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Surface area computed from triangle geometry.
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0_f64;
        for k in 0..self.triangle_count() {
            let a = self.position(self.indices[3 * k] as usize);
            let b = self.position(self.indices[3 * k + 1] as usize);
            let c = self.position(self.indices[3 * k + 2] as usize);
            area += 0.5 * (b - a).cross(c - a).norm();
        }
        area
    }

    fn position(&self, i: usize) -> Point3 {
        Point3::new(
            self.positions[3 * i],
            self.positions[3 * i + 1],
            self.positions[3 * i + 2],
        )
    }
}

#[derive(Default)]
pub(crate) struct SurfaceMeshBuilder {
    mesh: SurfaceMesh,
}

impl SurfaceMeshBuilder {
    pub(crate) fn vertex(&mut self, p: Point3) -> u32 {
        let i = self.mesh.vertex_count() as u32;
        self.mesh.positions.push(p.x);
        self.mesh.positions.push(p.y);
        self.mesh.positions.push(p.z);
        i
    }

    pub(crate) fn triangle(&mut self, a: u32, b: u32, c: u32) {
        if a == b || b == c || c == a {
            return;
        }
        self.mesh.indices.push(a);
        self.mesh.indices.push(b);
        self.mesh.indices.push(c);
    }

    pub(crate) fn finish(self) -> SurfaceMesh {
        self.mesh
    }
}
