//! `Brep` — the unit that owns a topology store together with its geometry.
//!
//! Topology ([`TopoStore`]) and geometry ([`GeomStore`]) are deliberately kept
//! in separate stores so the "topology carries no geometry" invariant holds
//! mechanically. A [`Brep`] bundles the two plus the list of top-level solids,
//! and provides the one-stop [`validate`](Brep::validate) entry point.

use crate::geom::{GeomStore, SurfaceGeom, VertexGeom};
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::validate::{validate_topology, Defect, ValidateLevel};
use crate::topo::{Sense, Solid, TopoStore};

/// A boundary representation: topology, geometry and the top-level solids.
#[derive(Debug, Clone, Default)]
pub struct Brep {
    /// The combinatorial topology.
    pub topo: TopoStore,
    /// The analytic geometry referenced by the topology.
    pub geom: GeomStore,
    /// The top-level solids that make up this B-rep.
    pub solids: Vec<Id<Solid>>,
}

impl Brep {
    /// Create an empty B-rep.
    pub fn new() -> Self {
        Self::default()
    }

    /// Validate this B-rep at the requested level.
    ///
    /// [`ValidateLevel::Light`] runs the pure topological checks only (Euler
    /// with the ring term, sibling completeness, loop continuity,
    /// resolvability). [`ValidateLevel::Full`] additionally runs the geometric
    /// consistency checks:
    ///
    /// * each half-edge boundary endpoint evaluates to its vertex coordinate;
    /// * consecutive half-edges in a loop meet geometrically;
    /// * every loop vertex lies on the face's surface.
    ///
    /// Returns every defect found, or `Ok(())` if there are none.
    pub fn validate(&self, tol: &Tol, level: ValidateLevel) -> Result<(), Vec<Defect>> {
        // Topology is always checked first; its genus is not asserted here
        // (callers that know the expected genus use the topology entry point
        // directly).
        let mut defects = match validate_topology(&self.topo, &self.solids, tol, None) {
            Ok(()) => Vec::new(),
            Err(d) => d,
        };

        if level == ValidateLevel::Full {
            self.check_geometry(tol, &mut defects);
        }

        if defects.is_empty() {
            Ok(())
        } else {
            Err(defects)
        }
    }

    /// Run the geometric consistency checks (Full only).
    fn check_geometry(&self, tol: &Tol, defects: &mut Vec<Defect>) {
        for solid_id in &self.solids {
            let Some(solid) = self.topo.solids.get(*solid_id) else {
                continue;
            };
            for shell_id in &solid.shells {
                let Some(shell) = self.topo.shells.get(*shell_id) else {
                    continue;
                };
                for face_id in &shell.faces {
                    let Some(face) = self.topo.faces.get(*face_id) else {
                        continue;
                    };
                    let surface = self.geom.surface(face.surface);
                    let mut loops = Vec::with_capacity(1 + face.inners.len());
                    loops.push(face.outer);
                    loops.extend(face.inners.iter().copied());
                    for loop_id in loops {
                        self.check_loop_geometry(loop_id, surface, tol, defects);
                    }
                    // `sense` does not affect coordinate checks; it is carried
                    // for the boolean phase. Silence the unused read.
                    let _ = matches!(face.sense, Sense::Same | Sense::Reversed);
                }
            }
        }
    }

    /// Geometric checks for a single loop.
    fn check_loop_geometry(
        &self,
        loop_id: Id<crate::topo::Loop>,
        surface: Option<&SurfaceGeom>,
        tol: &Tol,
        defects: &mut Vec<Defect>,
    ) {
        let Some(lp) = self.topo.loops.get(loop_id) else {
            return;
        };
        let n = lp.half_edges.len();
        if n == 0 {
            return;
        }
        for i in 0..n {
            let he_id = lp.half_edges[i];
            let next_id = lp.half_edges[(i + 1) % n];
            let (Some(he), Some(next)) = (
                self.topo.half_edges.get(he_id),
                self.topo.half_edges.get(next_id),
            ) else {
                continue;
            };
            let Some(curve) = self.geom.curve(he.curve) else {
                continue;
            };
            let start_pt = self.vertex_point(he.start);
            let end_pt = curve.point_at(he.boundary[1]);
            let begin_pt = curve.point_at(he.boundary[0]);

            // (b) boundary start endpoint matches the start vertex.
            if let Some(sp) = start_pt {
                let d = (begin_pt - sp).norm();
                if d > tol.length {
                    defects.push(Defect::BoundaryVertexMismatch {
                        half_edge: he_id,
                        distance: d,
                    });
                }
            }

            // (a) boundary end meets the next half-edge's start vertex.
            if let Some(nsp) = self.vertex_point(next.start) {
                let d = (end_pt - nsp).norm();
                if d > tol.length {
                    defects.push(Defect::LoopGeometryGap {
                        loop_id,
                        half_edge: he_id,
                        distance: d,
                    });
                }
            }

            // (c) the start vertex lies on the face surface.
            if let (Some(surf), Some(sp)) = (surface, start_pt) {
                let d = surf.signed_distance(sp).abs();
                if d > tol.length {
                    if let Some(point) = self.topo.vertices.get(he.start) {
                        defects.push(Defect::VertexOffSurface {
                            point: point.point,
                            distance: d,
                        });
                    }
                }
            }
        }
    }

    /// Resolve a vertex handle to its explicit coordinate, if available.
    fn vertex_point(&self, v: Id<crate::topo::Vertex>) -> Option<crate::math::Point3> {
        let vert = self.topo.vertices.get(v)?;
        match self.geom.point(vert.point)? {
            VertexGeom::Explicit(p) => Some(*p),
        }
    }
}
