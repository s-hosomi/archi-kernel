//! Vertex snapping / merging — the heart of degeneracy handling.
//!
//! Building geometry has "exact coincidence as the common case": shared edges,
//! a vertex landing on another edge, an opening touching a wall edge. The
//! arrangement only stays robust if all those near-coincident points become
//! **one** vertex before any topology is built. This module owns that merge.
//!
//! # Snap / merge specification
//!
//! * **Merge radius.** Two points are the same vertex iff they are within
//!   `tol.length` (Euclidean). Snapping is by *absorption*: each new point either
//!   joins an existing vertex (and is reported as that vertex's id) or creates a
//!   a new one. The stored coordinate of a vertex is the coordinate of the
//!   **first** point that created it — it is never moved by later absorptions.
//!   Fixing the representative this way makes the merge **order-stable for a
//!   given insertion order** and avoids the classic snap-rounding hazard of a
//!   moving target dragging in ever more points (`eps`-chaining is bounded:
//!   a point only merges if it is within `eps` of an *existing representative*,
//!   not of the transitive cluster).
//! * **No scale adaptation.** A single absolute `eps`, per the kernel's
//!   tolerance policy.
//! * **Lookup cost.** A uniform spatial hash with cell size `eps` keeps the
//!   query to its 3×3 neighbour cells, so insertion is expected O(1) and the
//!   whole snap pass is O(n) rather than O(n²).
//!
//! Non-transitivity of `eps`-equality is real (A≈B, B≈C, A≉C) and cannot be
//! erased here; it is *bounded* by fixing representatives and is then absorbed
//! downstream by the exact `orient2d` predicate, which never disagrees with
//! itself about a fixed vertex set.

use std::collections::HashMap;

use crate::boolean::poly2d::geom::{eps_sq, Point2};
use crate::tolerance::Tol;

/// Identifier of a merged vertex within a [`VertexStore`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VertexId(pub usize);

/// A spatial-hash-backed store that merges near-coincident points into a single
/// representative vertex.
#[derive(Debug, Clone)]
pub struct VertexStore {
    tol: Tol,
    inv_cell: f64,
    points: Vec<Point2>,
    /// Maps a grid cell to the vertex ids whose representative falls in it.
    grid: HashMap<(i64, i64), Vec<VertexId>>,
}

impl VertexStore {
    /// Create an empty store using the given tolerance as the merge radius.
    pub fn new(tol: Tol) -> Self {
        // Cell size = eps so any point within eps of a representative lies in
        // one of the 3×3 cells around the query cell.
        let cell = tol.length.max(f64::MIN_POSITIVE);
        Self {
            tol,
            inv_cell: 1.0 / cell,
            points: Vec::new(),
            grid: HashMap::new(),
        }
    }

    #[inline]
    fn cell_of(&self, p: Point2) -> (i64, i64) {
        (
            (p.x * self.inv_cell).floor() as i64,
            (p.y * self.inv_cell).floor() as i64,
        )
    }

    /// Insert a point, returning the id of the vertex it snapped to (existing or
    /// newly created).
    pub fn insert(&mut self, p: Point2) -> VertexId {
        let (cx, cy) = self.cell_of(p);
        let esq = eps_sq(&self.tol);
        let mut best: Option<(VertexId, f64)> = None;
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(ids) = self.grid.get(&(cx + dx, cy + dy)) {
                    for &id in ids {
                        let d = self.points[id.0].dist_sq(p);
                        if d <= esq {
                            match best {
                                Some((_, bd)) if bd <= d => {}
                                _ => best = Some((id, d)),
                            }
                        }
                    }
                }
            }
        }
        if let Some((id, _)) = best {
            return id;
        }
        let id = VertexId(self.points.len());
        self.points.push(p);
        self.grid.entry((cx, cy)).or_default().push(id);
        id
    }

    /// Coordinate of a stored vertex.
    #[inline]
    pub fn point(&self, id: VertexId) -> Point2 {
        self.points[id.0]
    }

    /// Number of distinct vertices after merging. Part of the store's natural
    /// API surface (and exercised by the merge tests); the lib build itself does
    /// not call it, hence the allow.
    #[inline]
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// `true` if the store has no vertices.
    #[inline]
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn near_points_merge() {
        let mut s = VertexStore::new(Tol::default());
        let a = s.insert(Point2::new(0.0_f64, 0.0_f64));
        let b = s.insert(Point2::new(1e-9_f64, 1e-9_f64));
        assert_eq!(a, b);
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn far_points_separate() {
        let mut s = VertexStore::new(Tol::default());
        let a = s.insert(Point2::new(0.0_f64, 0.0_f64));
        let b = s.insert(Point2::new(1.0_f64, 0.0_f64));
        assert_ne!(a, b);
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn representative_is_first_point() {
        let mut s = VertexStore::new(Tol::default());
        let a = s.insert(Point2::new(0.0_f64, 0.0_f64));
        let _ = s.insert(Point2::new(5e-7_f64, 0.0_f64));
        // Representative is the first point, never moved.
        assert_eq!(s.point(a), Point2::new(0.0_f64, 0.0_f64));
    }

    #[test]
    fn grid_neighbours_caught_across_cell_boundary() {
        // Two points straddling a cell boundary but within eps must still merge.
        let tol = Tol {
            length: 1e-6_f64,
            angular: 1e-9_f64,
        };
        let mut s = VertexStore::new(tol);
        let a = s.insert(Point2::new(0.9999996_f64, 0.0_f64));
        let b = s.insert(Point2::new(1.0000002_f64, 0.0_f64)); // 6e-7 apart < eps
        assert_eq!(a, b);
    }
}
