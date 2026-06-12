//! The mesh builder: vertex interning on quantised coordinates and triangle
//! accumulation.
//!
//! Interning is the mechanism that makes the mesh watertight. Every emitted
//! position is keyed on its quantised coordinate (the kernel-wide
//! [`QUANT_SCALE`](crate::boolean::support), the same scale the extruder and the
//! cut use), so two faces that reach the same 3-D point obtain the *same* vertex
//! index. Combined with per-curve edge sampling (so adjacent faces produce
//! *identical* points along a shared edge), this guarantees each mesh edge is
//! shared by exactly two oppositely-wound triangles (`DESIGN.md` §7).

use std::collections::HashMap;

use crate::boolean::support::{key, CoordKey};
use crate::math::Point3;

use super::Mesh;

/// Accumulates interned vertices and triangles into a [`Mesh`].
pub(crate) struct MeshBuilder {
    positions: Vec<f64>,
    indices: Vec<u32>,
    face_of: Vec<u32>,
    /// Quantised coordinate → vertex index, so identical points de-duplicate.
    index_of: HashMap<CoordKey, u32>,
}

impl MeshBuilder {
    /// A fresh, empty builder.
    pub(crate) fn new() -> Self {
        Self {
            positions: Vec::new(),
            indices: Vec::new(),
            face_of: Vec::new(),
            index_of: HashMap::new(),
        }
    }

    /// Intern `p`, returning the (possibly shared) vertex index.
    ///
    /// The first call for a given quantised coordinate stores `p` verbatim and
    /// mints a new index; later calls with the same key return it. The stored
    /// coordinate is whichever point arrived first — the two are within the
    /// quantisation step (≪ `Tol::length`), so the choice is immaterial.
    pub(crate) fn vertex(&mut self, p: Point3) -> u32 {
        let k = key(p);
        if let Some(&i) = self.index_of.get(&k) {
            return i;
        }
        let i = (self.positions.len() / 3) as u32;
        self.positions.push(p.x);
        self.positions.push(p.y);
        self.positions.push(p.z);
        self.index_of.insert(k, i);
        i
    }

    /// Emit a triangle `(a, b, c)` already in the desired winding, tagged with
    /// the source face's arena index.
    ///
    /// A degenerate triangle (two corners interned to the same vertex) is
    /// dropped: it contributes nothing to the surface and would create a
    /// spurious one-sided edge in the watertight count.
    pub(crate) fn triangle(&mut self, a: u32, b: u32, c: u32, face_tag: u32) {
        if a == b || b == c || c == a {
            return;
        }
        self.indices.push(a);
        self.indices.push(b);
        self.indices.push(c);
        self.face_of.push(face_tag);
    }

    /// Consume the builder into the finished [`Mesh`].
    pub(crate) fn finish(self) -> Mesh {
        Mesh {
            positions: self.positions,
            indices: self.indices,
            face_of: self.face_of,
        }
    }
}
