//! The arena collection that owns every topology entity.

use crate::topo::arena::{Arena, Id};
use crate::topo::{Face, HalfEdge, Loop, Shell, Solid, Vertex};

/// Owns all topology entities, each in its own generational arena.
///
/// Handles ([`Id`]) into one arena cannot be used against another, so the kind
/// of every reference is checked at compile time. The store holds no geometry;
/// geometry lives in [`GeomStore`](crate::geom::GeomStore).
#[derive(Debug, Clone, Default)]
pub struct TopoStore {
    /// Vertices.
    pub vertices: Arena<Vertex>,
    /// Half-edges.
    pub half_edges: Arena<HalfEdge>,
    /// Loops.
    pub loops: Arena<Loop>,
    /// Faces.
    pub faces: Arena<Face>,
    /// Shells.
    pub shells: Arena<Shell>,
    /// Solids.
    pub solids: Arena<Solid>,
}

impl TopoStore {
    /// Create an empty store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a vertex.
    pub fn add_vertex(&mut self, v: Vertex) -> Id<Vertex> {
        self.vertices.insert(v)
    }

    /// Insert a half-edge.
    pub fn add_half_edge(&mut self, he: HalfEdge) -> Id<HalfEdge> {
        self.half_edges.insert(he)
    }

    /// Insert a loop.
    pub fn add_loop(&mut self, l: Loop) -> Id<Loop> {
        self.loops.insert(l)
    }

    /// Insert a face.
    pub fn add_face(&mut self, f: Face) -> Id<Face> {
        self.faces.insert(f)
    }

    /// Insert a shell.
    pub fn add_shell(&mut self, s: Shell) -> Id<Shell> {
        self.shells.insert(s)
    }

    /// Insert a solid.
    pub fn add_solid(&mut self, s: Solid) -> Id<Solid> {
        self.solids.insert(s)
    }
}
