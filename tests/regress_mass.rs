//! Regression tests for the mass module (volume + formwork).
//!
//! Two bugs, both fixed in `src/mass/`:
//!
//! 1. `signed_volume` silently dropped a steep-cut cylinder side patch whose
//!    outer loop is `[Circle, Line, Circle, Ellipse]` (the `(2,1)` rim case),
//!    returning a far-too-small volume that disagreed with `signed_volume_checked`
//!    (`DESIGN.md` §6-4 "no silent zero"). `src/mass/volume.rs`.
//! 2. `formwork_area` dumped a whole cylinder side face into `side` without
//!    checking the axis orientation, so a *horizontal* round member lost its
//!    soffit (`bottom`). `src/mass/properties.rs`.

use std::f64::consts::PI;

use archi_kernel::boolean::{cut, CutResult, KeepSide};
use archi_kernel::build::extrude;
use archi_kernel::csg::Profile2d;
use archi_kernel::geom::{CurveGeom, SurfaceGeom};
use archi_kernel::mass::{formwork_area, signed_volume, signed_volume_checked, VolumeError};
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::{Line3, Plane};
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::ValidateLevel;

fn z_axis() -> Line3 {
    Line3::new(Point3::origin(), Vec3::Z).expect("z axis")
}

fn plane(p: Point3, n: Vec3) -> Plane {
    Plane::new(p, n).expect("plane")
}

/// A faithful numerical evaluation of the divergence-theorem volume of the
/// *stored* B-rep: `V = (1/3) ∮ x·n̂ dA`, by densely sampling every face's
/// boundary loop (curve edges sampled along their stored parameter range) and
/// summing the oriented shoelace area `½ Σ qᵢ × qᵢ₊₁` projected on each face's
/// outward normal, times the plane offset.
///
/// This is the reference the closed-form [`signed_volume`] must match *for the
/// geometry the B-rep actually encodes* (independent of how that geometry was
/// produced), so it pins the closed form without depending on any upstream
/// cut/parametrisation behaviour.
fn faithful_brep_volume(brep: &archi_kernel::brep::Brep) -> f64 {
    let mut three_v = 0.0_f64;
    for solid_id in &brep.solids {
        let s = brep.topo.solids.get(*solid_id).unwrap();
        for sh_id in &s.shells {
            let sh = brep.topo.shells.get(*sh_id).unwrap();
            for f_id in &sh.faces {
                let f = brep.topo.faces.get(*f_id).unwrap();
                let n_out = match brep.geom.surface(f.surface) {
                    Some(SurfaceGeom::Plane(p)) => {
                        let n = p.normal().as_vec();
                        match f.sense {
                            archi_kernel::topo::Sense::Reversed => -n,
                            _ => n,
                        }
                    }
                    // Cylinder faces: sample the patch directly below.
                    _ => {
                        three_v += faithful_cylinder_contribution(brep, f);
                        continue;
                    }
                };
                let mut loops = vec![f.outer];
                loops.extend(f.inners.iter().copied());
                for loop_id in loops {
                    let pts = sample_loop(brep, loop_id);
                    if pts.len() < 3 {
                        continue;
                    }
                    let o = Point3::origin();
                    let mut acc = Vec3::ZERO;
                    for i in 0..pts.len() {
                        let a = pts[i] - o;
                        let b = pts[(i + 1) % pts.len()] - o;
                        acc = acc + a.cross(b);
                    }
                    let area = 0.5 * acc.dot(n_out);
                    let d = (pts[0] - o).dot(n_out);
                    three_v += d * area;
                }
            }
        }
    }
    three_v / 3.0
}

/// Densely sample a face's boundary loop into points, walking each half-edge's
/// curve over its stored parameter range (the same range the kernel tessellates).
fn sample_loop(
    brep: &archi_kernel::brep::Brep,
    loop_id: archi_kernel::topo::arena::Id<archi_kernel::topo::Loop>,
) -> Vec<Point3> {
    let lp = brep.topo.loops.get(loop_id).unwrap();
    let mut pts = Vec::new();
    for &he_id in &lp.half_edges {
        let he = brep.topo.half_edges.get(he_id).unwrap();
        let [a, b] = he.boundary;
        match brep.geom.curve(he.curve) {
            Some(CurveGeom::Line(_)) => {
                let v = brep.topo.vertices.get(he.start).unwrap();
                pts.push(brep.geom.point(v.point).and_then(|g| g.as_point()).unwrap());
            }
            Some(c @ (CurveGeom::Circle(_) | CurveGeom::Ellipse(_))) => {
                let steps = 4000usize;
                for s in 0..steps {
                    let t = a + (b - a) * (s as f64) / (steps as f64);
                    let p = match c {
                        CurveGeom::Circle(circ) => circ.point_at(t),
                        CurveGeom::Ellipse(e) => e.point_at(t),
                        _ => unreachable!(),
                    };
                    pts.push(p);
                }
            }
            _ => {}
        }
    }
    pts
}

/// Faithful `∫ x·n̂ dA` over a cylinder patch face, by sampling the surface on a
/// `(φ, z)` grid bounded below by the bottom rim and above by each upper rim.
/// Used only as an independent reference for the closed-form cylinder integral.
fn faithful_cylinder_contribution(
    brep: &archi_kernel::brep::Brep,
    face: &archi_kernel::topo::Face,
) -> f64 {
    // Integrate ∫ x·n̂ dA over the ruled cylinder patch by a (φ, z) sweep: the
    // patch is bounded below by the bottom rim (z = 0 in the rim frame) and above
    // by the first cap plane the ruling at angle φ reaches. With n̂(φ) ⟂ axis the
    // integrand is constant in z, so ∫₀^{z₁} x·n̂ dz = (base·n̂)·z₁.
    let SurfaceGeom::Cylinder(cyl) = brep.geom.surface(face.surface).unwrap() else {
        return 0.0;
    };
    let axis = cyl.axis().dir().as_vec();
    let r = cyl.radius();
    // Build (u, v) the same deterministic way the kernel's `plane_basis` does.
    let seed = if axis.x.abs() <= axis.y.abs() && axis.x.abs() <= axis.z.abs() {
        Vec3::X
    } else if axis.y.abs() <= axis.z.abs() {
        Vec3::Y
    } else {
        Vec3::Z
    };
    let u = axis.cross(seed);
    let u = u * (1.0 / u.dot(u).sqrt());
    let v = axis.cross(u);

    // The lowest circular rim (smallest axial projection) is the patch bottom.
    let lp = brep.topo.loops.get(face.outer).unwrap();
    let mut bottom: Option<(f64, Point3)> = None;
    for &he_id in &lp.half_edges {
        let he = brep.topo.half_edges.get(he_id).unwrap();
        if let Some(CurveGeom::Circle(c)) = brep.geom.curve(he.curve) {
            let proj = (c.center() - Point3::origin()).dot(axis);
            if bottom.map(|(p, _)| proj < p).unwrap_or(true) {
                bottom = Some((proj, c.center()));
            }
        }
    }
    let (_, bc) = bottom.unwrap();

    // Cap planes the rulings rise to: every non-bottom circle's horizontal plane
    // and every ellipse's cut plane.
    let mut caps: Vec<Plane> = Vec::new();
    for &he_id in &lp.half_edges {
        let he = brep.topo.half_edges.get(he_id).unwrap();
        match brep.geom.curve(he.curve) {
            Some(CurveGeom::Circle(c)) => {
                let proj = (c.center() - Point3::origin()).dot(axis);
                let base = (bc - Point3::origin()).dot(axis);
                if proj > base + 1e-9 {
                    caps.push(Plane::new(c.center(), axis).unwrap());
                }
            }
            Some(CurveGeom::Ellipse(e)) => {
                caps.push(Plane::new(e.center(), e.normal().as_vec()).unwrap());
            }
            _ => {}
        }
    }

    // Bottom arc angular range in (u, v): collect φ of the bottom rim vertices.
    let mut phis: Vec<f64> = Vec::new();
    for &he_id in &lp.half_edges {
        let he = brep.topo.half_edges.get(he_id).unwrap();
        if let Some(CurveGeom::Circle(c)) = brep.geom.curve(he.curve) {
            let proj = (c.center() - Point3::origin()).dot(axis);
            let base = (bc - Point3::origin()).dot(axis);
            if proj <= base + 1e-9 {
                phis.push(he.boundary[0]);
                phis.push(he.boundary[1]);
            }
        }
    }
    let phi_lo = phis.iter().cloned().fold(f64::INFINITY, f64::min);
    let phi_hi = phis.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let n = 20000usize;
    let mut three_v = 0.0_f64;
    let dphi = (phi_hi - phi_lo) / n as f64;
    for i in 0..n {
        let phi = phi_lo + dphi * (i as f64 + 0.5);
        let nrm = u * phi.cos() + v * phi.sin();
        let base_pt = bc + nrm * r;
        // z₁(φ) = min over caps of the axial crossing height above the bottom.
        let mut z1 = f64::INFINITY;
        for cap in &caps {
            let np = cap.normal().as_vec();
            let denom = np.dot(axis);
            if denom.abs() < 1e-12 {
                continue;
            }
            let z = np.dot(cap.point() - base_pt) / denom;
            if z >= -1e-9 && z < z1 {
                z1 = z;
            }
        }
        if !z1.is_finite() || z1 <= 0.0 {
            continue;
        }
        // ∫_0^{z1} x·n̂ dz with x = base_pt + z·axis, n̂ = nrm (axis·n̂ = 0):
        // x·n̂ = base_pt·n̂ (constant in z) ⇒ contribution = (base_pt·n̂) · z1.
        let xn = (base_pt - Point3::origin()).dot(nrm);
        three_v += xn * z1 * r * dphi;
    }
    three_v
}

/// Bug 1: a steep oblique cut digs past the bottom rim, so the cylinder side
/// face's outer loop becomes `[Circle, Line, Circle, Ellipse]` (the `(2,1)` rim
/// case). The fixed `signed_volume` integrates this patch in closed form, so it
/// matches the Monte-Carlo truth and agrees with `signed_volume_checked` (no
/// `UnsupportedCylinderFace`).
#[test]
fn steep_cut_through_bottom_rim_clce_volume() {
    let tol = Tol::default();
    let radius = 0.3_f64;
    let length = 3.0_f64;
    let profile = Profile2d::circle(radius).expect("circle");
    let brep = extrude(&profile, &z_axis(), length, &tol).expect("extrude cyl");
    let solid = brep.solids[0];

    // Steep oblique plane through (0,0,1.5) with normal (0,8,1): the cut digs down
    // to the bottom rim, so the cylinder side patch becomes [C, L, C, E].
    let cut_plane = plane(Point3::new(0.0, 0.0, 1.5), Vec3::new(0.0, 8.0, 1.0));
    let below = cut(&brep, solid, &cut_plane, KeepSide::Below, &tol).expect("below cut");
    let CutResult::Cut { brep: bb, .. } = below else {
        panic!("expected a Cut result");
    };
    bb.validate(&tol, ValidateLevel::Full)
        .expect("steep cut valid");

    // Confirm at least one cylinder side face's outer loop really is the
    // [C, L, C, E] (2,1) rim case (the steep cut clips the bottom circle with an
    // axial line into two arcs; the seam splits the wall into half-patches, so
    // there may be more than one such face). This guards that the test exercises
    // the (2,1) path rather than some easier configuration.
    let mut clce_face_count = 0usize;
    for solid_id in &bb.solids {
        let s = bb.topo.solids.get(*solid_id).unwrap();
        for sh_id in &s.shells {
            let sh = bb.topo.shells.get(*sh_id).unwrap();
            for f_id in &sh.faces {
                let f = bb.topo.faces.get(*f_id).unwrap();
                if let Some(SurfaceGeom::Cylinder(_)) = bb.geom.surface(f.surface) {
                    let lp = bb.topo.loops.get(f.outer).unwrap();
                    let kinds: Vec<&str> = lp
                        .half_edges
                        .iter()
                        .map(|&he_id| {
                            let he = bb.topo.half_edges.get(he_id).unwrap();
                            match bb.geom.curve(he.curve) {
                                Some(CurveGeom::Circle(_)) => "C",
                                Some(CurveGeom::Ellipse(_)) => "E",
                                Some(CurveGeom::Line(_)) => "L",
                                _ => "?",
                            }
                        })
                        .collect();
                    if kinds == vec!["C", "L", "C", "E"] {
                        clce_face_count += 1;
                    }
                }
            }
        }
    }
    assert!(
        clce_face_count >= 1,
        "expected at least one [C, L, C, E] cylinder face, found {clce_face_count}"
    );

    let v = signed_volume(&bb);
    let checked = signed_volume_checked(&bb);

    // (a) No silent zero: the lenient `signed_volume` no longer drops the
    // `[C, L, C, E]` cylinder patch, so it agrees with the checked variant — they
    // disagreed (0.1245 vs Err(UnsupportedCylinderFace)) before the fix.
    assert_eq!(
        checked,
        Ok(v),
        "checked variant must match lenient and not error (no silent zero)"
    );

    // (b) The closed form is exact for the geometry the B-rep encodes: it matches
    // a faithful, densely-sampled divergence-theorem integral of the same stored
    // faces (planar loops sampled along their stored curve params, cylinder
    // patches swept on a (φ, z) grid). This pins the cylinder (2,1) integral and
    // the new ellipse-arc planar-segment correction without depending on any
    // upstream cut/parametrisation behaviour.
    let faithful = faithful_brep_volume(&bb);
    assert!(
        (v - faithful).abs() < 5e-3_f64,
        "signed_volume {v} disagrees with the faithful sampled integral {faithful} \
         (diff {})",
        (v - faithful).abs()
    );

    // (c) The cylinder side really contributes the bulk that used to be dropped:
    // the recovered value is far from the old silently-truncated 0.1245.
    assert!(
        v > 0.2_f64,
        "signed_volume {v} should recover the dropped cylinder patch (was ~0.1245)"
    );

    // Touch the error type so it stays part of the public surface under test.
    let _ = VolumeError::UnsupportedCylinderFace;
}

/// Bug 2: a *horizontal* round member's side face is split into `bottom` (the
/// lower-half soffit, ≈ πrL) and "no formwork" (the upper half); the two vertical
/// end caps remain `side`. Before the fix the whole lateral surface plus both
/// caps were dumped into `side`, with `bottom == 0`.
#[test]
fn horizontal_round_member_formwork_split() {
    let tol = Tol::default();
    let r = 0.3_f64;
    let len = 2.0_f64;
    let x_axis = Line3::new(Point3::origin(), Vec3::X).expect("axis");
    let brep = extrude(&Profile2d::circle(r).expect("circle"), &x_axis, len, &tol)
        .expect("horizontal round member");

    let fw = formwork_area(&brep, &tol).expect("formwork");

    let caps = 2.0 * PI * r * r; // 2πr², the two vertical end caps
    let correct_bottom = PI * r * len; // lower half of the lateral wall (soffit)
    let correct_side = caps; // only the end caps are vertical

    assert!(
        (fw.bottom - correct_bottom).abs() <= 1e-6_f64,
        "bottom = {}, expected soffit {correct_bottom}",
        fw.bottom
    );
    assert!(
        (fw.side - correct_side).abs() <= 1e-6_f64,
        "side = {}, expected just the end caps {correct_side}",
        fw.side
    );
}

/// Control: a *vertical* round column's whole wall is `side` and there is no
/// `bottom` (no downward-facing curved face). This guards the fix from
/// over-classifying a genuine column wall as soffit.
#[test]
fn vertical_round_column_formwork_all_side() {
    let tol = Tol::default();
    let r = 0.3_f64;
    let len = 2.0_f64;
    let brep = extrude(&Profile2d::circle(r).expect("circle"), &z_axis(), len, &tol)
        .expect("vertical round column");

    let fw = formwork_area(&brep, &tol).expect("formwork");

    // The lateral wall is vertical ⇒ all side. The end caps are horizontal: the
    // bottom cap faces down (a footing soffit ⇒ bottom), the top cap faces up
    // (no formwork). So bottom == one cap area, side == lateral wall.
    let lateral = 2.0 * PI * r * len;
    let cap = PI * r * r;
    assert!(
        (fw.side - lateral).abs() <= 1e-6_f64,
        "side = {}, expected the lateral wall {lateral}",
        fw.side
    );
    assert!(
        (fw.bottom - cap).abs() <= 1e-6_f64,
        "bottom = {}, expected the bottom cap {cap}",
        fw.bottom
    );
}
