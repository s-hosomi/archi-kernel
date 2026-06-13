//! wasm-bindgen bindings for archi-kernel.
//!
//! The boundary follows `DESIGN.md` §8: geometry crosses into JavaScript as
//! flat typed arrays (`Float32Array` positions, `Uint32Array` indices), and
//! everything structured (CSG trees in, quantities/diagnostics out) crosses as
//! JSON strings — the kernel's own `serde` representation is the wire format.
//!
//! The kernel itself stays oblivious to wasm: this crate is a thin adapter and
//! deliberately contains no geometry logic beyond arc sampling for display.

use archi_kernel::clash::{clash_check, ClashKind, ClashOptions};
use archi_kernel::csg::{CsgNode, Member, StableId};
use archi_kernel::curved::{
    tessellate_cone_panel, tessellate_cylinder_panel, tessellate_sphere_panel,
    tessellate_thick_cylinder_panel, tessellate_thick_sphere_panel, ConePanel, ConePanelOptions,
    ConePanelSpec, CylinderPanel, CylinderPanelOptions, SpherePanel, SpherePanelOptions,
    SpherePanelSpec, ThickCylinderPanel, ThickSpherePanel, TrimLoop2d,
};
use archi_kernel::model::{takeoff, Model};
use archi_kernel::primitives::{Cylinder, Line3, Plane};
use archi_kernel::section::{section, SectionEdge, SectionResult};
use archi_kernel::tess::{tessellate, TessOptions};
use archi_kernel::tolerance::Tol;
use archi_kernel::{Point3, Vec3};
use std::collections::HashMap;
use wasm_bindgen::prelude::*;

/// One member's merged render mesh (all solids of the member, one buffer).
#[wasm_bindgen]
pub struct MeshData {
    positions: Vec<f32>,
    indices: Vec<u32>,
}

#[wasm_bindgen]
impl MeshData {
    /// Vertex positions, `xyz` interleaved, metres.
    #[wasm_bindgen(getter)]
    pub fn positions(&self) -> js_sys::Float32Array {
        js_sys::Float32Array::from(&self.positions[..])
    }

    /// Triangle indices into [`positions`](Self::positions).
    #[wasm_bindgen(getter)]
    pub fn indices(&self) -> js_sys::Uint32Array {
        js_sys::Uint32Array::from(&self.indices[..])
    }
}

/// Serializable summary of one member's evaluation, for the UI.
#[derive(serde::Serialize)]
struct MemberStatus {
    id: u64,
    ok: bool,
    error: Option<String>,
    volume: Option<f64>,
}

/// Serializable section output: world-space polylines per profile.
#[derive(serde::Serialize)]
struct SectionLoopOut {
    /// Closed polyline, `xyz` triples (the closing edge back to the first
    /// point is implicit). Arcs are sampled for display.
    points: Vec<[f64; 3]>,
}

#[derive(serde::Serialize)]
struct SectionProfileOut {
    outer: SectionLoopOut,
    holes: Vec<SectionLoopOut>,
}

#[derive(serde::Serialize)]
struct MemberSectionOut {
    id: u64,
    profiles: Vec<SectionProfileOut>,
}

#[derive(serde::Serialize)]
struct SectionErrorOut {
    id: u64,
    error: String,
}

#[derive(serde::Serialize)]
struct SectionAllOut {
    members: Vec<MemberSectionOut>,
    errors: Vec<SectionErrorOut>,
}

#[derive(serde::Serialize)]
struct ClashOut {
    a: u64,
    b: u64,
    kind: String,
    volume: Option<f64>,
}

/// A kernel model held behind the wasm boundary.
///
/// Members are inserted as JSON-encoded [`CsgNode`] trees (the kernel's serde
/// representation), evaluated lazily, and read back as meshes, sections,
/// quantity take-offs and clash reports.
#[wasm_bindgen]
pub struct KernelModel {
    model: Model,
    curved: HashMap<StableId, CurvedPanelInput>,
    tol: Tol,
}

impl Default for KernelModel {
    fn default() -> Self {
        Self::new()
    }
}

#[wasm_bindgen]
impl KernelModel {
    /// An empty model with the default kernel tolerance (`1e-6 m`).
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            model: Model::new(),
            curved: HashMap::new(),
            tol: Tol::default(),
        }
    }

    /// Insert a member from its JSON-encoded CSG tree. `id` must be fresh.
    pub fn insert(&mut self, id: u64, csg_json: &str) -> Result<(), JsError> {
        let node: CsgNode =
            serde_json::from_str(csg_json).map_err(|e| JsError::new(&format!("bad CSG: {e}")))?;
        self.model
            .insert(StableId(id), Member::new(node))
            .map_err(|e| JsError::new(&format!("insert: {e:?}")))
    }

    /// Insert or replace a curved analytic panel from JSON.
    ///
    /// Curved panels are a renderable surface layer at this phase: they can be
    /// tessellated through [`curved_mesh`](Self::curved_mesh), but are not yet
    /// part of the CSG/B-rep evaluator used by sections, take-off or clashes.
    pub fn insert_curved(&mut self, id: u64, curved_json: &str) -> Result<(), JsError> {
        let node: CurvedPanelInput = serde_json::from_str(curved_json)
            .map_err(|e| JsError::new(&format!("bad curved panel: {e}")))?;
        self.curved.insert(StableId(id), node);
        Ok(())
    }

    /// Member ids currently in the model, ascending.
    pub fn ids(&self) -> Vec<u64> {
        self.model.ids().map(|s| s.0).collect()
    }

    /// Curved panel ids currently registered in the render layer, ascending.
    pub fn curved_ids(&self) -> Vec<u64> {
        let mut ids: Vec<u64> = self.curved.keys().map(|s| s.0).collect();
        ids.sort_unstable();
        ids
    }

    /// Evaluate every member; returns a JSON array of per-member statuses
    /// (id, ok, error, volume). Failures are member-local by design.
    pub fn evaluate_all(&mut self) -> Result<String, JsError> {
        let mut out = Vec::new();
        for (id, result) in self.model.evaluate_all(&self.tol) {
            let status = match result {
                Ok(brep) => MemberStatus {
                    id: id.0,
                    ok: true,
                    error: None,
                    volume: archi_kernel::mass::signed_volume_checked(&brep).ok(),
                },
                Err(e) => MemberStatus {
                    id: id.0,
                    ok: false,
                    error: Some(format!("{e:?}")),
                    volume: None,
                },
            };
            out.push(status);
        }
        serde_json::to_string(&out).map_err(|e| JsError::new(&e.to_string()))
    }

    /// Tessellate one member (all of its solids) into a render mesh.
    ///
    /// `chord_tolerance` is the arc-to-chord deviation in metres (display
    /// quality; `0.002` is a good default for building scale).
    pub fn mesh(&mut self, id: u64, chord_tolerance: f64) -> Result<MeshData, JsError> {
        let brep = self
            .model
            .evaluate(StableId(id), &self.tol)
            .map_err(|e| JsError::new(&format!("evaluate {id}: {e:?}")))?;
        let opts = TessOptions::with_chord_tolerance(chord_tolerance);
        let mut positions: Vec<f32> = Vec::new();
        let mut indices: Vec<u32> = Vec::new();
        for &solid in &brep.solids {
            let mesh = tessellate(&brep, solid, &opts, &self.tol)
                .map_err(|e| JsError::new(&format!("tessellate {id}: {e:?}")))?;
            let base = (positions.len() / 3) as u32;
            positions.extend(mesh.positions.iter().map(|&v| v as f32));
            indices.extend(mesh.indices.iter().map(|&i| i + base));
        }
        Ok(MeshData { positions, indices })
    }

    /// Tessellate one registered curved analytic panel into a render mesh.
    pub fn curved_mesh(&self, id: u64, chord_tolerance: f64) -> Result<MeshData, JsError> {
        let node = self
            .curved
            .get(&StableId(id))
            .ok_or_else(|| JsError::new(&format!("unknown curved panel id: {id}")))?;
        if chord_tolerance <= 0.0 || !chord_tolerance.is_finite() {
            return Err(JsError::new("chord_tolerance must be positive and finite"));
        }
        let mesh = node.mesh(chord_tolerance, &self.tol)?;
        Ok(surface_mesh_to_mesh_data(mesh))
    }

    /// Quantity take-off of one member, as JSON
    /// (`concrete_volume`, `formwork_side`, `formwork_bottom`).
    pub fn takeoff(&mut self, id: u64) -> Result<String, JsError> {
        let q = takeoff(&mut self.model, StableId(id), &self.tol)
            .map_err(|e| JsError::new(&format!("takeoff {id}: {e:?}")))?;
        serde_json::to_string(&q).map_err(|e| JsError::new(&e.to_string()))
    }

    /// Section every member with the plane through `(px,py,pz)` with normal
    /// `(nx,ny,nz)`. Returns a JSON array of per-member profiles whose loops
    /// are world-space closed polylines (arcs sampled for display).
    /// Members the plane misses contribute no entry; members that fail to
    /// evaluate are skipped (they already report through `evaluate_all`).
    pub fn section_all(
        &mut self,
        px: f64,
        py: f64,
        pz: f64,
        nx: f64,
        ny: f64,
        nz: f64,
    ) -> Result<String, JsError> {
        let plane = Plane::new(Point3::new(px, py, pz), Vec3::new(nx, ny, nz))
            .map_err(|e| JsError::new(&format!("bad plane: {e}")))?;
        let ids: Vec<u64> = self.ids();
        let mut members: Vec<MemberSectionOut> = Vec::new();
        let mut errors: Vec<SectionErrorOut> = Vec::new();
        for id in ids {
            let Ok(brep) = self.model.evaluate(StableId(id), &self.tol) else {
                continue; // evaluation failures already report via evaluate_all
            };
            let mut profiles = Vec::new();
            for &solid in &brep.solids {
                match section(&brep, solid, &plane, &self.tol) {
                    Ok(result) => collect_profiles(&result, &mut profiles),
                    Err(e) => errors.push(SectionErrorOut {
                        id,
                        error: format!("{e:?}"),
                    }),
                }
            }
            if !profiles.is_empty() {
                members.push(MemberSectionOut { id, profiles });
            }
        }
        let out = SectionAllOut { members, errors };
        serde_json::to_string(&out).map_err(|e| JsError::new(&e.to_string()))
    }

    /// Run the model-wide clash check; returns a JSON array of
    /// `{a, b, kind, volume}` (kind: `"hard" | "touching" | "potential"`).
    pub fn clash(&mut self) -> Result<String, JsError> {
        let report = clash_check(&mut self.model, &self.tol, &ClashOptions::default());
        let out: Vec<ClashOut> = report
            .clashes
            .iter()
            .map(|c| {
                let (kind, volume) = match c.kind {
                    ClashKind::HardClash { volume } => ("hard", Some(volume)),
                    ClashKind::Touching => ("touching", None),
                    ClashKind::PotentialClash => ("potential", None),
                    _ => ("other", None),
                };
                ClashOut {
                    a: c.a.0,
                    b: c.b.0,
                    kind: kind.to_string(),
                    volume,
                }
            })
            .collect();
        serde_json::to_string(&out).map_err(|e| JsError::new(&e.to_string()))
    }
}

fn surface_mesh_to_mesh_data(mesh: archi_kernel::curved::SurfaceMesh) -> MeshData {
    MeshData {
        positions: mesh.positions.iter().map(|&v| v as f32).collect(),
        indices: mesh.indices,
    }
}

#[derive(Debug, Clone, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
enum CurvedPanelInput {
    Cylinder {
        axis_origin: Point3Input,
        axis_dir: Point3Input,
        radius: f64,
        theta_min: f64,
        theta_max: f64,
        z_min: f64,
        z_max: f64,
        #[serde(default)]
        holes: Vec<TrimLoopInput>,
        thickness: Option<f64>,
    },
    Sphere {
        center: Point3Input,
        radius: f64,
        pole: Point3Input,
        theta_min: f64,
        theta_max: f64,
        phi_min: f64,
        phi_max: f64,
        #[serde(default)]
        holes: Vec<TrimLoopInput>,
        thickness: Option<f64>,
    },
    Cone {
        apex: Point3Input,
        axis: Point3Input,
        half_angle: f64,
        theta_min: f64,
        theta_max: f64,
        height_min: f64,
        height_max: f64,
        #[serde(default)]
        holes: Vec<TrimLoopInput>,
        thickness: Option<f64>,
    },
}

impl CurvedPanelInput {
    fn mesh(
        &self,
        chord_tolerance: f64,
        tol: &Tol,
    ) -> Result<archi_kernel::curved::SurfaceMesh, JsError> {
        match self {
            CurvedPanelInput::Cylinder {
                axis_origin,
                axis_dir,
                radius,
                theta_min,
                theta_max,
                z_min,
                z_max,
                holes,
                thickness,
            } => {
                let axis = Line3::new(axis_origin.point(), axis_dir.vec())
                    .map_err(|e| JsError::new(&format!("cylinder axis: {e}")))?;
                let cylinder = Cylinder::new(axis, *radius)
                    .map_err(|e| JsError::new(&format!("cylinder: {e}")))?;
                let panel = CylinderPanel::new(
                    cylinder,
                    *theta_min,
                    *theta_max,
                    *z_min,
                    *z_max,
                    trim_loops(holes, tol)?,
                    tol,
                )
                .map_err(curved_error)?;
                let opts = CylinderPanelOptions::with_chord_tolerance(chord_tolerance);
                if let Some(thickness) = thickness {
                    let thick = ThickCylinderPanel::new(panel, *thickness).map_err(curved_error)?;
                    tessellate_thick_cylinder_panel(&thick, &opts, tol).map_err(curved_error)
                } else {
                    tessellate_cylinder_panel(&panel, &opts, tol).map_err(curved_error)
                }
            }
            CurvedPanelInput::Sphere {
                center,
                radius,
                pole,
                theta_min,
                theta_max,
                phi_min,
                phi_max,
                holes,
                thickness,
            } => {
                let pole = pole
                    .vec()
                    .try_unit()
                    .ok_or_else(|| JsError::new("sphere pole must be a finite non-zero vector"))?;
                let panel = SpherePanel::new(
                    SpherePanelSpec {
                        center: center.point(),
                        radius: *radius,
                        pole,
                        theta_min: *theta_min,
                        theta_max: *theta_max,
                        phi_min: *phi_min,
                        phi_max: *phi_max,
                    },
                    trim_loops(holes, tol)?,
                    tol,
                )
                .map_err(curved_error)?;
                let opts = SpherePanelOptions::with_chord_tolerance(chord_tolerance);
                if let Some(thickness) = thickness {
                    let thick = ThickSpherePanel::new(panel, *thickness).map_err(curved_error)?;
                    tessellate_thick_sphere_panel(&thick, &opts, tol).map_err(curved_error)
                } else {
                    tessellate_sphere_panel(&panel, &opts, tol).map_err(curved_error)
                }
            }
            CurvedPanelInput::Cone {
                apex,
                axis,
                half_angle,
                theta_min,
                theta_max,
                height_min,
                height_max,
                holes,
                thickness,
            } => {
                if thickness.is_some() {
                    return Err(JsError::new("thick cone panels are not supported yet"));
                }
                let axis = axis
                    .vec()
                    .try_unit()
                    .ok_or_else(|| JsError::new("cone axis must be a finite non-zero vector"))?;
                let panel = ConePanel::new(
                    ConePanelSpec {
                        apex: apex.point(),
                        axis,
                        half_angle: *half_angle,
                        theta_min: *theta_min,
                        theta_max: *theta_max,
                        height_min: *height_min,
                        height_max: *height_max,
                    },
                    trim_loops(holes, tol)?,
                    tol,
                )
                .map_err(curved_error)?;
                tessellate_cone_panel(
                    &panel,
                    &ConePanelOptions::with_chord_tolerance(chord_tolerance),
                    tol,
                )
                .map_err(curved_error)
            }
        }
    }
}

#[derive(Debug, Clone, Copy, serde::Deserialize)]
struct Point3Input {
    x: f64,
    y: f64,
    z: f64,
}

impl Point3Input {
    fn point(self) -> Point3 {
        Point3::new(self.x, self.y, self.z)
    }

    fn vec(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }
}

#[derive(Debug, Clone, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
enum TrimLoopInput {
    Rectangle {
        u_min: f64,
        u_max: f64,
        v_min: f64,
        v_max: f64,
        #[serde(default)]
        reverse: bool,
    },
    Circle {
        center: [f64; 2],
        radius: f64,
        #[serde(default)]
        reverse: bool,
    },
    Polygon {
        points: Vec<[f64; 2]>,
        #[serde(default)]
        reverse: bool,
    },
}

fn trim_loops(inputs: &[TrimLoopInput], tol: &Tol) -> Result<Vec<TrimLoop2d>, JsError> {
    inputs.iter().map(|input| trim_loop(input, tol)).collect()
}

fn trim_loop(input: &TrimLoopInput, tol: &Tol) -> Result<TrimLoop2d, JsError> {
    let (loop_, reverse) = match input {
        TrimLoopInput::Rectangle {
            u_min,
            u_max,
            v_min,
            v_max,
            reverse,
        } => (
            TrimLoop2d::rectangle(*u_min, *u_max, *v_min, *v_max, tol).map_err(curved_error)?,
            *reverse,
        ),
        TrimLoopInput::Circle {
            center,
            radius,
            reverse,
        } => (
            TrimLoop2d::circle(*center, *radius, tol).map_err(curved_error)?,
            *reverse,
        ),
        TrimLoopInput::Polygon { points, reverse } => (
            TrimLoop2d::from_points(points, tol).map_err(curved_error)?,
            *reverse,
        ),
    };
    Ok(if reverse { loop_.reversed() } else { loop_ })
}

fn curved_error(e: archi_kernel::curved::CurvedError) -> JsError {
    JsError::new(&format!("curved panel: {e}"))
}

/// Flatten a kernel [`SectionResult`] into world-space display polylines.
fn collect_profiles(result: &SectionResult, out: &mut Vec<SectionProfileOut>) {
    let Some(frame) = &result.frame else {
        return;
    };
    let sample = |edges: &[SectionEdge]| -> SectionLoopOut {
        let mut points: Vec<[f64; 3]> = Vec::new();
        for edge in edges {
            match edge {
                SectionEdge::Line { start, .. } => {
                    let p = frame.to_3d(*start);
                    points.push([p.x, p.y, p.z]);
                }
                SectionEdge::Arc {
                    center,
                    radius,
                    start_angle,
                    end_angle,
                    ..
                } => {
                    // Sample the arc finely enough for smooth display.
                    let sweep = end_angle - start_angle;
                    let n = ((sweep.abs() / (std::f64::consts::TAU / 64.0)).ceil() as usize).max(2);
                    for k in 0..n {
                        let t = *start_angle + sweep * (k as f64) / (n as f64);
                        let x = center[0] + radius * t.cos();
                        let y = center[1] + radius * t.sin();
                        let p = frame.to_3d([x, y]);
                        points.push([p.x, p.y, p.z]);
                    }
                }
                _ => {}
            }
        }
        SectionLoopOut { points }
    };
    for profile in &result.profiles {
        out.push(SectionProfileOut {
            outer: sample(&profile.outer.edges),
            holes: profile.holes.iter().map(|h| sample(&h.edges)).collect(),
        });
    }
}
