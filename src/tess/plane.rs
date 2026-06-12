//! Planar-face tessellation: ear-clipping with classic hole bridging.
//!
//! A planar face is an outer loop plus zero or more interior (hole) loops. We
//!
//! 1. polyline each loop in 3-D — straight edges contribute their two
//!    endpoints, arc edges are sampled at the shared per-curve density so the
//!    points coincide with the adjacent face's (`DESIGN.md` §6-5);
//! 2. project the rings to the face plane's 2-D frame;
//! 3. orient the outer ring CCW and each hole ring CW *in the frame whose
//!    `+z = outward normal*`, so the ear-clipper's CCW output triangles wind to
//!    the face's outward normal;
//! 4. bridge every hole into the outer ring (classic hole-bridging: connect a
//!    hole's rightmost vertex to a visible outer vertex with a doubled
//!    seam edge), turning the holed polygon into one simple polygon;
//! 5. ear-clip the simple polygon and emit each triangle through the shared
//!    coordinate intern.
//!
//! Every 2-D vertex keeps the 3-D point it came from so emission interns the
//! true coordinate (shared with adjacent faces) rather than a re-lift.

use crate::brep::Brep;
use crate::geom::CurveGeom;
use crate::math::{Point3, Vec3};
use crate::primitives::Plane;
use crate::tolerance::Tol;
use crate::topo::arena::Id;
use crate::topo::{Face, Loop, Sense};

use super::intern::MeshBuilder;
use super::{arc_segment_count, TessError, TessOptions};

/// A boundary vertex carrying both its 2-D frame coordinate and its 3-D point.
#[derive(Debug, Clone, Copy)]
struct RingVertex {
    xy: [f64; 2],
    p: Point3,
}

/// Tessellate one planar face into the builder.
pub(crate) fn tessellate_plane_face(
    brep: &Brep,
    face: &Face,
    plane: &Plane,
    builder: &mut MeshBuilder,
    face_tag: u32,
    opts: &TessOptions,
    _tol: &Tol,
) -> Result<(), TessError> {
    // Outward normal: the plane normal, flipped when the face sense reverses it.
    let n_out = match face.sense {
        Sense::Same => plane.normal().as_vec(),
        Sense::Reversed => -plane.normal().as_vec(),
    };

    // A 2-D frame (u, v) with u × v = n_out, so CCW-in-frame ⇒ CCW-about-n_out.
    let (u, v) = frame_for(n_out);
    let origin = plane.point();
    let project = |p: Point3| -> [f64; 2] {
        let d = p - origin;
        [d.dot(u), d.dot(v)]
    };

    // Polyline the outer loop and each hole into 3-D rings, projected to 2-D.
    let mut outer = ring_of(brep, face.outer, &project, opts)?;
    if outer.len() < 3 {
        return Err(TessError::DegenerateFace);
    }
    // Outer must be CCW in the n_out frame; reverse if not.
    if signed_area(&outer) < 0.0 {
        outer.reverse();
    }

    let mut holes: Vec<Vec<RingVertex>> = Vec::with_capacity(face.inners.len());
    for &inner in &face.inners {
        let mut hole = ring_of(brep, inner, &project, opts)?;
        if hole.len() < 3 {
            // A degenerate hole encloses no area; skip it rather than fail (it
            // contributes no surface).
            continue;
        }
        // Holes must be CW in the n_out frame (opposite the outer).
        if signed_area(&hole) > 0.0 {
            hole.reverse();
        }
        holes.push(hole);
    }

    // Bridge every hole into the outer polygon, producing one simple polygon.
    let polygon = bridge_holes(outer, holes);

    // Ear-clip and emit. Triangles come out CCW in the n_out frame, i.e. wound
    // to the face's outward normal.
    ear_clip(&polygon, builder, face_tag)
}

/// Build an orthonormal `(u, v)` with `u × v = n` (a unit-length `n`).
///
/// Deterministic seed pick (least-aligned axis) mirroring
/// [`plane_basis`](crate::primitives::plane_basis), but here keyed to an
/// arbitrary outward normal rather than a stored plane, so the winding frame is
/// well-defined for either face sense.
fn frame_for(n: Vec3) -> (Vec3, Vec3) {
    let nu = n.try_unit().map(|u| u.as_vec()).unwrap_or(Vec3::Z);
    let seed = if nu.x.abs() <= nu.y.abs() && nu.x.abs() <= nu.z.abs() {
        Vec3::X
    } else if nu.y.abs() <= nu.z.abs() {
        Vec3::Y
    } else {
        Vec3::Z
    };
    let u = nu
        .cross(seed)
        .try_unit()
        .map(|u| u.as_vec())
        .unwrap_or(Vec3::X);
    let v = nu.cross(u);
    (u, v)
}

/// Polyline a loop into a ring of [`RingVertex`]es (no duplicated closing
/// point).
///
/// Each half-edge contributes its start vertex and, for an arc, the interior
/// sample points up to (but not including) its end — the end is the next
/// half-edge's start, so the ring closes without duplication. The arc sample
/// count comes from the shared per-curve density so a sibling face samples the
/// same points.
fn ring_of(
    brep: &Brep,
    loop_id: Id<Loop>,
    project: &impl Fn(Point3) -> [f64; 2],
    opts: &TessOptions,
) -> Result<Vec<RingVertex>, TessError> {
    let lp = brep
        .topo
        .loops
        .get(loop_id)
        .ok_or(TessError::DanglingReference)?;
    let mut ring: Vec<RingVertex> = Vec::new();
    for &he_id in &lp.half_edges {
        let he = brep
            .topo
            .half_edges
            .get(he_id)
            .ok_or(TessError::DanglingReference)?;
        let curve = brep
            .geom
            .curve(he.curve)
            .ok_or(TessError::DanglingReference)?;
        let [a, b] = he.boundary;
        match curve {
            CurveGeom::Line(line) => {
                let p = line.point_at(a);
                ring.push(RingVertex { xy: project(p), p });
            }
            CurveGeom::Circle(_) | CurveGeom::Ellipse(_) => {
                let radius = match curve {
                    CurveGeom::Circle(c) => c.radius(),
                    CurveGeom::Ellipse(e) => e.semi_major(),
                    CurveGeom::Line(_) => unreachable!(),
                };
                let segs = arc_segment_count(radius, b - a, opts.chord_tolerance);
                // Sample t at a, a+step, … (segs points: a .. just before b); the
                // endpoint b is the next half-edge's start.
                for s in 0..segs {
                    let t = a + (b - a) * (s as f64) / (segs as f64);
                    let p = curve.point_at(t);
                    ring.push(RingVertex { xy: project(p), p });
                }
            }
        }
    }
    Ok(ring)
}

/// Signed area of a 2-D ring (shoelace); positive for CCW.
fn signed_area(ring: &[RingVertex]) -> f64 {
    let n = ring.len();
    let mut acc = 0.0_f64;
    for i in 0..n {
        let a = ring[i].xy;
        let b = ring[(i + 1) % n].xy;
        acc += a[0] * b[1] - b[0] * a[1];
    }
    acc / 2.0_f64
}

/// Merge every hole into the outer ring by a bridge edge, yielding one simple
/// polygon (classic hole-bridging).
///
/// Holes are processed in decreasing order of their bridge vertex's x so an
/// inner hole never bridges across an already-merged outer one. For each hole we
/// pick its rightmost vertex `M`, find a mutually visible vertex `P` on the
/// current outer ring, and splice the hole in as
/// `… P, M, hole(rotated to start at M)…, M, P …`. The doubled `P`–`M` seam has
/// zero area, so it does not perturb the triangulation, and the ear-clipper
/// simply walks past it.
fn bridge_holes(outer: Vec<RingVertex>, mut holes: Vec<Vec<RingVertex>>) -> Vec<RingVertex> {
    // Bridge the hole whose rightmost vertex is farthest right first.
    holes.sort_by(|a, b| {
        rightmost(b)
            .partial_cmp(&rightmost(a))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut poly = outer;
    for hole in holes {
        poly = bridge_one(poly, &hole);
    }
    poly
}

/// The maximum x of a ring's vertices (its rightmost extent).
fn rightmost(ring: &[RingVertex]) -> f64 {
    ring.iter()
        .map(|rv| rv.xy[0])
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Splice a single hole into `outer` via a visible bridge.
fn bridge_one(outer: Vec<RingVertex>, hole: &[RingVertex]) -> Vec<RingVertex> {
    // M = the hole's rightmost vertex (ties broken by larger y for determinism).
    let m_idx = (0..hole.len())
        .max_by(|&i, &j| {
            let a = hole[i].xy;
            let b = hole[j].xy;
            a[0].partial_cmp(&b[0])
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a[1].partial_cmp(&b[1]).unwrap_or(std::cmp::Ordering::Equal))
        })
        .unwrap_or(0);
    let m = hole[m_idx];

    // Find the outer vertex that best serves as the bridge endpoint: the one
    // closest to M to the right (a ray-cast toward +x picks the nearest visible
    // edge; we then snap to that edge's vertex with the larger x that does not
    // make the bridge cross any outer edge). A robust-enough heuristic for the
    // building cases: choose the visible outer vertex minimising the bridge
    // length, breaking ties by smaller angle.
    let p_idx = best_bridge_vertex(&outer, hole, m);

    // Build: outer[0..=p], M and the hole rotated to start at M (full loop),
    // back to M, then P again, then outer[p+1..].
    let mut out = Vec::with_capacity(outer.len() + hole.len() + 2);
    out.extend_from_slice(&outer[..=p_idx]);
    // Hole walked starting at M, all the way round, ending back at M.
    for k in 0..hole.len() {
        out.push(hole[(m_idx + k) % hole.len()]);
    }
    out.push(m); // close the hole loop back to M
    out.push(outer[p_idx]); // return seam to P
    out.extend_from_slice(&outer[p_idx + 1..]);
    out
}

/// Choose the outer-ring vertex to bridge a hole vertex `m` to.
///
/// We want a vertex `P` of `outer` such that the segment `P–M` does not cross
/// any edge of `outer` or `hole` (mutual visibility). Among visible candidates
/// we pick the nearest. The building faces this tessellator meets (rectangular
/// caps with rectangular / round openings) always have a clear horizontal sight
/// line, so a visibility scan with a nearest-distance tie-break is sufficient;
/// if none tests visible (a pathological projection), we fall back to the
/// nearest vertex outright so the routine never panics.
fn best_bridge_vertex(outer: &[RingVertex], hole: &[RingVertex], m: RingVertex) -> usize {
    let mut best: Option<(usize, f64)> = None;
    let mut nearest: (usize, f64) = (0, f64::INFINITY);
    for (i, p) in outer.iter().enumerate() {
        let d = sq_dist(p.xy, m.xy);
        if d < nearest.1 {
            nearest = (i, d);
        }
        if bridge_visible(outer, hole, p.xy, m.xy) {
            match best {
                Some((_, bd)) if bd <= d => {}
                _ => best = Some((i, d)),
            }
        }
    }
    best.map(|(i, _)| i).unwrap_or(nearest.0)
}

/// Squared distance between two 2-D points.
fn sq_dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy
}

/// `true` if the segment `p–m` crosses no edge of `outer` or `hole` (excluding
/// edges that share the endpoint `p` or `m`).
fn bridge_visible(outer: &[RingVertex], hole: &[RingVertex], p: [f64; 2], m: [f64; 2]) -> bool {
    for ring in [outer, hole] {
        let n = ring.len();
        for i in 0..n {
            let a = ring[i].xy;
            let b = ring[(i + 1) % n].xy;
            // Skip edges incident to either bridge endpoint (sharing a point is
            // not a crossing).
            if same_pt(a, p) || same_pt(b, p) || same_pt(a, m) || same_pt(b, m) {
                continue;
            }
            if segments_cross(p, m, a, b) {
                return false;
            }
        }
    }
    true
}

/// Coordinate equality at the proper-crossing scale (well below `Tol::length`).
fn same_pt(a: [f64; 2], b: [f64; 2]) -> bool {
    sq_dist(a, b) <= 1e-18_f64
}

/// `true` if open segments `p1–p2` and `p3–p4` properly cross.
fn segments_cross(p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], p4: [f64; 2]) -> bool {
    let d1 = orient(p3, p4, p1);
    let d2 = orient(p3, p4, p2);
    let d3 = orient(p1, p2, p3);
    let d4 = orient(p1, p2, p4);
    (d1 * d2 < 0.0) && (d3 * d4 < 0.0)
}

/// Orientation determinant of `(a, b, c)` (positive = CCW turn).
fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Ear-clip a simple CCW polygon, emitting triangles through the builder.
///
/// Classic O(n²) ear clipping: repeatedly find a convex vertex whose triangle
/// contains no other vertex (an "ear"), emit it, and remove the ear tip. The
/// polygon is CCW (the caller guarantees it), so emitted triangles are CCW and
/// wind to the face's outward normal.
///
/// Degenerate handling: collinear (zero-area) vertices are valid ear tips and
/// are clipped away (the resulting zero-area triangles are dropped by the
/// builder). If no ear is found in a full pass (a numerically awkward ring) we
/// clip the least-bad vertex to guarantee termination rather than loop forever.
fn ear_clip(
    polygon: &[RingVertex],
    builder: &mut MeshBuilder,
    face_tag: u32,
) -> Result<(), TessError> {
    let n = polygon.len();
    if n < 3 {
        return Err(TessError::DegenerateFace);
    }
    // Work on indices into `polygon`.
    let mut idx: Vec<usize> = (0..n).collect();

    let mut guard = 0usize;
    let guard_max = n * n + 16;
    while idx.len() > 3 {
        guard += 1;
        if guard > guard_max {
            // Numerical fallback: emit a fan over whatever remains so we never
            // hang, then stop. The fan may not be perfectly clean for a wildly
            // non-convex remnant, but the building cases never reach here.
            emit_fan(polygon, &idx, builder, face_tag);
            return Ok(());
        }
        let m = idx.len();
        let mut clipped = false;
        for i in 0..m {
            let ip = idx[(i + m - 1) % m];
            let ic = idx[i];
            let in_ = idx[(i + 1) % m];
            let a = polygon[ip].xy;
            let b = polygon[ic].xy;
            let c = polygon[in_].xy;
            // Convex (CCW) corner?
            if orient(a, b, c) <= 0.0 {
                continue;
            }
            // No other vertex strictly inside the candidate ear? Vertices that
            // coincide with one of the ear's corners (the doubled bridge seam
            // endpoints) are skipped — they sit on the boundary and must not
            // block the ear, which is exactly the robustness fix the classic
            // hole-bridge pinch needs.
            let mut contains = false;
            for &k in &idx {
                if k == ip || k == ic || k == in_ {
                    continue;
                }
                let q = polygon[k].xy;
                if same_pt(q, a) || same_pt(q, b) || same_pt(q, c) {
                    continue;
                }
                if point_in_triangle_strict(q, a, b, c) {
                    contains = true;
                    break;
                }
            }
            if contains {
                continue;
            }
            // It is an ear: emit and clip.
            emit_tri(builder, polygon, [ip, ic, in_], face_tag);
            idx.remove(i);
            clipped = true;
            break;
        }
        if !clipped {
            // No ear found this pass: clip the first convex-ish vertex to make
            // progress (handles a sliver the strict test rejected).
            let i = pick_fallback_ear(polygon, &idx);
            let m = idx.len();
            let ip = idx[(i + m - 1) % m];
            let ic = idx[i];
            let in_ = idx[(i + 1) % m];
            emit_tri(builder, polygon, [ip, ic, in_], face_tag);
            idx.remove(i);
        }
    }
    // Final triangle.
    if idx.len() == 3 {
        emit_tri(builder, polygon, [idx[0], idx[1], idx[2]], face_tag);
    }
    Ok(())
}

/// Intern the three polygon vertices `[i, j, k]` and emit the triangle.
fn emit_tri(builder: &mut MeshBuilder, polygon: &[RingVertex], tri: [usize; 3], face_tag: u32) {
    let a = builder.vertex(polygon[tri[0]].p);
    let b = builder.vertex(polygon[tri[1]].p);
    let c = builder.vertex(polygon[tri[2]].p);
    builder.triangle(a, b, c, face_tag);
}

/// Emit a triangle fan over the remaining indices (numerical fallback).
fn emit_fan(polygon: &[RingVertex], idx: &[usize], builder: &mut MeshBuilder, face_tag: u32) {
    if idx.len() < 3 {
        return;
    }
    let p0 = polygon[idx[0]].p;
    let i0 = builder.vertex(p0);
    for w in 1..idx.len() - 1 {
        let i1 = builder.vertex(polygon[idx[w]].p);
        let i2 = builder.vertex(polygon[idx[w + 1]].p);
        builder.triangle(i0, i1, i2, face_tag);
    }
}

/// Pick a fallback ear when the strict scan found none: the most convex vertex.
fn pick_fallback_ear(polygon: &[RingVertex], idx: &[usize]) -> usize {
    let m = idx.len();
    let mut best = 0usize;
    let mut best_turn = f64::NEG_INFINITY;
    for i in 0..m {
        let a = polygon[idx[(i + m - 1) % m]].xy;
        let b = polygon[idx[i]].xy;
        let c = polygon[idx[(i + 1) % m]].xy;
        let turn = orient(a, b, c);
        if turn > best_turn {
            best_turn = turn;
            best = i;
        }
    }
    best
}

/// `true` if `p` lies **strictly** inside triangle `(a, b, c)`.
///
/// A point on an edge or at a vertex is *not* inside: a vertex on the ear's
/// boundary (a coincident bridge-seam endpoint, or a collinear neighbour) must
/// not block the ear, otherwise the doubled-seam pinch stalls the clipper.
fn point_in_triangle_strict(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    let d1 = orient(a, b, p);
    let d2 = orient(b, c, p);
    let d3 = orient(c, a, p);
    let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    // Strictly inside ⇒ all three determinants share one strict sign.
    !(has_neg && has_pos) && d1 != 0.0 && d2 != 0.0 && d3 != 0.0
}
