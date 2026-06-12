//! Watertight tessellation — Phase 6 (`DESIGN.md` §6-5, §7, §10 Phase 6).
//!
//! Turns a closed [`Brep`] solid into an indexed triangle [`Mesh`] for display
//! (Three.js) and as the front end of FEM mesh generation. The output is
//! **watertight**: every interior edge of the mesh is shared by exactly two
//! triangles (`DESIGN.md` §7, "各エッジちょうど 2 三角形").
//!
//! # Why the result is structurally watertight
//!
//! The watertight property is not checked-and-hoped-for; it is built in. The
//! two ingredients (`DESIGN.md` §6-5):
//!
//! 1. **Per-curve edge discretisation.** A curved edge (a circle / ellipse arc)
//!    is discretised once *per curve*, from a sample count derived in closed
//!    form from the chord tolerance. The two sibling half-edges that share a
//!    curve have reversed boundary intervals over the *same* parameter range, so
//!    sampling the curve at the shared parameters yields the *same* point
//!    sequence for both faces — only the traversal direction differs.
//! 2. **Coordinate interning.** Every emitted position is interned on its
//!    quantised coordinate ([`crate::boolean::support`]'s shared `QUANT_SCALE`),
//!    so two faces that reach the same 3-D point reach the same vertex index.
//!    Adjacent faces therefore reference identical indices along a shared edge,
//!    which makes every mesh edge appear in exactly two triangles with opposite
//!    orientation.
//!
//! Straight edges trivially share their two endpoints; the interning handles
//! them with no special case.
//!
//! # Orientation
//!
//! Each face is triangulated with the winding that agrees with the face's
//! outward normal (the surface normal folded through the face
//! [`Sense`](crate::topo::Sense)). The signed volume of the resulting mesh is
//! then positive for a correctly oriented closed solid, the same orientation
//! test Phase 2 uses (`signed_volume(mesh) > 0`).
//!
//! # Scope
//!
//! Planar faces (with holes) and cylinder faces (straight or obliquely cut,
//! with elliptical rims) are tessellated — the surfaces the extruder, the cut
//! and the prismatic engine produce. A surface kind outside that set is
//! reported as [`TessError::UnsupportedSurface`] rather than silently skipped
//! (`DESIGN.md` §6-4: no silent zero).

mod cylinder;
mod intern;
mod plane;

use crate::brep::Brep;
use crate::geom::SurfaceGeom;
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::Solid;

use intern::MeshBuilder;

/// An indexed triangle mesh with a flat coordinate layout for the wasm /
/// Three.js boundary (`DESIGN.md` §8: 境界はフラット f64 配列).
///
/// `positions` is a flat `[x0, y0, z0, x1, y1, z1, …]` array; vertex `i`
/// occupies `positions[3i .. 3i + 3]`. `indices` lists triangle corners three at
/// a time (every triangle is `indices[3k], indices[3k+1], indices[3k+2]`), each
/// an index into the vertex array. `face_of[k]` is the arena index of the
/// [`Face`](crate::topo::Face) triangle `k` came from, for FEM face tagging and
/// per-face display (`DESIGN.md` §6-5).
#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Mesh {
    /// Flat `xyz` vertex coordinates, in metres (SI).
    pub positions: Vec<f64>,
    /// Triangle corner indices, three per triangle.
    pub indices: Vec<u32>,
    /// For each triangle, the arena index of the source [`Face`](crate::topo::Face).
    pub face_of: Vec<u32>,
}

impl Mesh {
    /// The number of vertices (`positions.len() / 3`).
    pub fn vertex_count(&self) -> usize {
        self.positions.len() / 3
    }

    /// The number of triangles (`indices.len() / 3`).
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// The signed volume of the mesh by the divergence theorem, in cubic metres.
    ///
    /// `V = (1/6) Σ_tri (a · (b × c))` over the triangles `(a, b, c)`. For a
    /// closed, outward-oriented mesh this is positive and equals the enclosed
    /// volume; the sign is the orientation check (`DESIGN.md` §6-5,
    /// `signed_volume(mesh) > 0`).
    pub fn signed_volume(&self) -> f64 {
        let mut acc = 0.0_f64;
        for k in 0..self.triangle_count() {
            let [ia, ib, ic] = [
                self.indices[3 * k] as usize,
                self.indices[3 * k + 1] as usize,
                self.indices[3 * k + 2] as usize,
            ];
            let a = self.position(ia);
            let b = self.position(ib);
            let c = self.position(ic);
            // a · (b × c)
            let cross = [
                b[1] * c[2] - b[2] * c[1],
                b[2] * c[0] - b[0] * c[2],
                b[0] * c[1] - b[1] * c[0],
            ];
            acc += a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2];
        }
        acc / 6.0_f64
    }

    /// The `[x, y, z]` coordinate of vertex `i`.
    fn position(&self, i: usize) -> [f64; 3] {
        [
            self.positions[3 * i],
            self.positions[3 * i + 1],
            self.positions[3 * i + 2],
        ]
    }
}

/// Options controlling tessellation density.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct TessOptions {
    /// Maximum allowed chord error for curved edges and surfaces, in metres.
    ///
    /// An arc of radius `r` is split so that the deviation of each chord from
    /// the arc (the sagitta) does not exceed this value; the closed-form sample
    /// step is `Δφ ≤ 2·acos(1 − chord_tolerance / r)` (see
    /// [`arc_segment_count`]). Smaller values give a denser, more faithful mesh
    /// at the cost of more triangles.
    pub chord_tolerance: f64,
}

impl TessOptions {
    /// Options with the given chord tolerance (metres), all other knobs at their
    /// defaults.
    ///
    /// `TessOptions` is `#[non_exhaustive]` (future density knobs are
    /// semver-additive), so callers outside the crate build it through this
    /// constructor or [`Default`] rather than a struct literal.
    pub fn with_chord_tolerance(chord_tolerance: f64) -> Self {
        Self { chord_tolerance }
    }
}

impl Default for TessOptions {
    fn default() -> Self {
        // 1e-3 m (1 mm). For display the eye does not resolve a sub-millimetre
        // chord deviation at architectural viewing distances, and a typical
        // round column (r ≈ 0.2–0.6 m) then gets a few dozen facets — smooth
        // enough to read as round without flooding the GPU. This is a *display*
        // default; an FEM front end dials it down via `chord_tolerance`
        // (`DESIGN.md` §6-5).
        Self {
            chord_tolerance: 1e-3_f64,
        }
    }
}

/// A face configuration the tessellator cannot triangulate.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum TessError {
    /// A face's surface handle did not resolve, or its surface kind is not
    /// supported by the tessellator.
    UnsupportedSurface,
    /// A planar face could not be triangulated: its boundary, projected to the
    /// face plane, is degenerate (fewer than three distinct vertices on the
    /// outer loop, or a self-touching ring the ear-clipper could not resolve).
    DegenerateFace,
    /// A handle (loop, half-edge, vertex, curve) referenced by a face did not
    /// resolve in the store.
    DanglingReference,
    /// A `chord_tolerance` that is not strictly positive (it would imply an
    /// unbounded number of segments).
    NonPositiveChordTolerance {
        /// The offending value.
        value: f64,
    },
}

/// Tessellate `solid` of `brep` into a watertight triangle [`Mesh`].
///
/// Every face reachable from `solid` is triangulated with a winding consistent
/// with its outward normal, sharing vertices on the quantised-coordinate intern
/// so the result is structurally watertight (module docs). The mesh's signed
/// volume is positive for a correctly oriented closed solid.
///
/// # Errors
///
/// * [`TessError::NonPositiveChordTolerance`] if `opts.chord_tolerance ≤ 0`.
/// * [`TessError::UnsupportedSurface`] for a surface kind the tessellator does
///   not handle.
/// * [`TessError::DegenerateFace`] for a planar face whose boundary is
///   degenerate after projection.
/// * [`TessError::DanglingReference`] if a handle does not resolve.
pub fn tessellate(
    brep: &Brep,
    solid: Id<Solid>,
    opts: &TessOptions,
    tol: &Tol,
) -> Result<Mesh, TessError> {
    // Reject a non-positive (or NaN) chord tolerance: it would imply an
    // unbounded segment count. Written as `>= ` of the negation so NaN, which
    // fails every ordered comparison, is also rejected.
    if opts.chord_tolerance <= 0.0 || opts.chord_tolerance.is_nan() {
        return Err(TessError::NonPositiveChordTolerance {
            value: opts.chord_tolerance,
        });
    }

    let mut builder = MeshBuilder::new();
    let solid_ref = brep
        .topo
        .solids
        .get(solid)
        .ok_or(TessError::DanglingReference)?;

    for &shell_id in &solid_ref.shells {
        let shell = brep
            .topo
            .shells
            .get(shell_id)
            .ok_or(TessError::DanglingReference)?;
        for &face_id in &shell.faces {
            let face = brep
                .topo
                .faces
                .get(face_id)
                .ok_or(TessError::DanglingReference)?;
            let surface = brep
                .geom
                .surface(face.surface)
                .ok_or(TessError::UnsupportedSurface)?;
            let face_tag = face_id.index();
            match surface {
                SurfaceGeom::Plane(plane) => {
                    plane::tessellate_plane_face(
                        brep,
                        face,
                        plane,
                        &mut builder,
                        face_tag,
                        opts,
                        tol,
                    )?;
                }
                SurfaceGeom::Cylinder(cyl) => {
                    cylinder::tessellate_cylinder_face(
                        brep,
                        face,
                        cyl,
                        &mut builder,
                        face_tag,
                        opts,
                        tol,
                    )?;
                }
            }
        }
    }

    Ok(builder.finish())
}

/// Number of straight segments to split an arc of radius `r` sweeping `|sweep|`
/// radians into, so each chord deviates from the arc by at most `chord_tol`.
///
/// The sagitta of a chord subtending `Δφ` on radius `r` is `r(1 − cos(Δφ/2))`;
/// requiring it `≤ chord_tol` gives `Δφ ≤ 2·acos(1 − chord_tol / r)` (the
/// closed form in `DESIGN.md` §6-5). The segment count is
/// `⌈|sweep| / Δφ_max⌉`, clamped to at least 1. When `chord_tol ≥ r` (the chord
/// may deviate by the whole radius) the bound saturates and a single segment per
/// quadrant-ish suffices; we still return at least 1.
pub(crate) fn arc_segment_count(radius: f64, sweep: f64, chord_tol: f64) -> usize {
    let sweep = sweep.abs();
    if sweep <= 0.0 {
        return 1;
    }
    if radius <= 0.0 {
        return 1;
    }
    // 1 − chord_tol/r; if chord_tol ≥ r the arccos argument would leave [−1, 1],
    // meaning one chord already satisfies the tolerance.
    let x = 1.0_f64 - chord_tol / radius;
    if x <= -1.0 {
        return 1;
    }
    let dphi_max = 2.0_f64 * x.clamp(-1.0, 1.0).acos();
    if dphi_max <= 0.0 {
        // chord_tol ≈ 0: fall back to a fine but finite split rather than ∞.
        return ((sweep / 1e-3_f64).ceil() as usize).max(1);
    }
    ((sweep / dphi_max).ceil() as usize).max(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arc_segment_count_is_finite_and_monotone() {
        // Finer tolerance ⇒ at least as many segments.
        let coarse = arc_segment_count(0.5_f64, std::f64::consts::PI, 1e-2_f64);
        let fine = arc_segment_count(0.5_f64, std::f64::consts::PI, 1e-4_f64);
        assert!(fine >= coarse);
        assert!(coarse >= 1);
    }

    #[test]
    fn arc_segment_count_handles_degenerate_inputs() {
        assert_eq!(arc_segment_count(0.0_f64, 1.0_f64, 1e-3_f64), 1);
        assert_eq!(arc_segment_count(1.0_f64, 0.0_f64, 1e-3_f64), 1);
        // chord tolerance larger than radius: one segment is enough.
        assert_eq!(arc_segment_count(0.1_f64, std::f64::consts::PI, 1.0_f64), 1);
    }
}
