//! Phase 3c — prismatic circular-void / sleeve tests (`DESIGN.md` §10 Phase 3c,
//! §4.2). The commercial-must family: a wall pierced by a circular sleeve, a
//! blind circular pocket, multiple sleeves, a circular column cut by a tool, and
//! a horizontal section through a circular void (arc loop). Every literal carries
//! an `f64` annotation and an explicit tolerance (`DESIGN.md` §12).

use std::f64::consts::PI;

use archi_kernel::boolean::prismatic::{self, ExtrudeLeaf, PrismError};
use archi_kernel::csg::Profile2d;
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::Plane;
use archi_kernel::section::{section, SectionEdge};
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::ValidateLevel;

const VOL_EPS: f64 = 1e-9;

/// A wall: box extruded along `+x`, 3.0 long (x), 0.2 thick (y), 2.0 tall (z),
/// base corner at the world origin.
fn wall() -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(1.0_f64, 0.1_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 0.1_f64, 1.0_f64),
        axis: Vec3::X,
        length: 3.0_f64,
    }
}

/// A circular sleeve through the wall thickness (along `y`), centred at `(x, z)`,
/// radius `r`, length `len` (≥ 0.2 to pass through).
fn sleeve(x: f64, z: f64, r: f64, len: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::circle(r).expect("circle"),
        origin: Point3::new(x, -0.1_f64, z),
        axis: Vec3::Y,
        length: len,
    }
}

// ── milestone: wall − circular sleeve (through) ───────────────────────────────

#[test]
fn wall_minus_circular_sleeve_through() {
    let tol = Tol::default();
    let r = 0.05_f64;
    let result = prismatic::difference(&wall(), &sleeve(1.5_f64, 1.0_f64, r, 0.4_f64), &tol)
        .expect("circular sleeve difference");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight wall with circular hole");
    assert_eq!(result.solids.len(), 1usize);
    // V = V_wall − πr²·t, with t = wall thickness 0.2.
    let expected = 3.0_f64 * 0.2_f64 * 2.0_f64 - PI * r * r * 0.2_f64;
    assert!(
        (result.signed_volume() - expected).abs() <= VOL_EPS,
        "volume = {} (expected {expected})",
        result.signed_volume()
    );
}

#[test]
fn wall_minus_blind_circular_pocket() {
    let tol = Tol::default();
    let r = 0.05_f64;
    // A pocket only 0.1 deep into the 0.2-thick wall.
    let pocket = ExtrudeLeaf {
        profile: Profile2d::circle(r).expect("circle"),
        origin: Point3::new(1.5_f64, 0.0_f64, 1.0_f64),
        axis: Vec3::Y,
        length: 0.1_f64,
    };
    let result = prismatic::difference(&wall(), &pocket, &tol).expect("blind pocket");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight blind pocket");
    assert_eq!(result.solids.len(), 1usize, "a blind pocket does not split");
    let expected = 1.2_f64 - PI * r * r * 0.1_f64;
    assert!(
        (result.signed_volume() - expected).abs() <= VOL_EPS,
        "volume = {} (expected {expected})",
        result.signed_volume()
    );
}

#[test]
fn wall_minus_two_separate_sleeves() {
    let tol = Tol::default();
    let r = 0.05_f64;
    let result = prismatic::opening_subtraction(
        &wall(),
        &[
            sleeve(0.8_f64, 1.0_f64, r, 0.4_f64),
            sleeve(2.2_f64, 1.0_f64, r, 0.4_f64),
        ],
        &tol,
    )
    .expect("two sleeves");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight two sleeves");
    let expected = 1.2_f64 - 2.0_f64 * PI * r * r * 0.2_f64;
    assert!(
        (result.signed_volume() - expected).abs() <= VOL_EPS,
        "volume = {} (expected {expected})",
        result.signed_volume()
    );
}

#[test]
fn wall_minus_overlapping_sleeves_fuse() {
    let tol = Tol::default();
    let r = 0.1_f64;
    // Two circular sleeves whose discs overlap fuse into one merged void.
    let result = prismatic::opening_subtraction(
        &wall(),
        &[
            sleeve(1.45_f64, 1.0_f64, r, 0.4_f64),
            sleeve(1.55_f64, 1.0_f64, r, 0.4_f64),
        ],
        &tol,
    )
    .expect("overlapping sleeves");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight fused sleeves");
    // The fused void removes less than two full discs (the overlap is counted
    // once), so the volume is strictly above wall − 2·disc·t.
    let two_discs = 1.2_f64 - 2.0_f64 * PI * r * r * 0.2_f64;
    assert!(
        result.signed_volume() > two_discs,
        "overlap must be counted once: {} vs {two_discs}",
        result.signed_volume()
    );
    assert!(result.signed_volume() < 1.2_f64);
}

#[test]
fn wall_minus_mixed_circular_and_rect_openings() {
    let tol = Tol::default();
    let r = 0.05_f64;
    let rect_open = ExtrudeLeaf {
        profile: Profile2d::rect(0.05_f64, 0.2_f64).expect("rect"),
        origin: Point3::new(2.2_f64, -0.1_f64, 1.0_f64),
        axis: Vec3::Y,
        length: 0.4_f64,
    };
    let result = prismatic::opening_subtraction(
        &wall(),
        &[sleeve(0.8_f64, 1.0_f64, r, 0.4_f64), rect_open],
        &tol,
    )
    .expect("mixed openings");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight mixed openings");
    // Removed: one circular sleeve (πr²·0.2) + one rect opening (0.1·0.4·0.2).
    let expected = 1.2_f64 - PI * r * r * 0.2_f64 - 0.1_f64 * 0.4_f64 * 0.2_f64;
    assert!(
        (result.signed_volume() - expected).abs() <= VOL_EPS,
        "volume = {} (expected {expected})",
        result.signed_volume()
    );
}

#[test]
fn circular_column_minus_rect_tool_common_axis() {
    let tol = Tol::default();
    // A round column along z, cut by a rectangular tool sweeping along z (common
    // direction = column axis z). The result is a circular column with a flat.
    let r = 0.3_f64;
    let column = ExtrudeLeaf {
        profile: Profile2d::circle(r).expect("circle"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 1.0_f64,
    };
    // A box tool removing the right side (x > 0.15), full z height.
    let tool = ExtrudeLeaf {
        profile: Profile2d::rect(0.5_f64, 0.5_f64).expect("rect"),
        origin: Point3::new(0.65_f64, 0.0_f64, -0.5_f64),
        axis: Vec3::Z,
        length: 2.0_f64,
    };
    let result = prismatic::difference(&column, &tool, &tol).expect("column minus tool");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight cut column");
    // The result is a circular segment column: less than the full disc, > 0.
    assert!(result.signed_volume() > 0.0_f64);
    assert!(result.signed_volume() < PI * r * r * 1.0_f64);
}

#[test]
fn horizontal_section_of_circular_void_has_arc_loop() {
    let tol = Tol::default();
    // A square column with a vertical circular void (both along z). A horizontal
    // section meets the cylindrical void in arcs (`SectionEdge::Arc`).
    let column = ExtrudeLeaf {
        profile: Profile2d::rect(0.5_f64, 0.5_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::Z,
        length: 2.0_f64,
    };
    let void = ExtrudeLeaf {
        profile: Profile2d::circle(0.2_f64).expect("circle"),
        origin: Point3::new(0.0_f64, 0.0_f64, -0.5_f64),
        axis: Vec3::Z,
        length: 3.0_f64,
    };
    let result = prismatic::difference(&column, &void, &tol).expect("column with void");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight voided column");

    let plane = Plane::new(Point3::new(0.0_f64, 0.0_f64, 1.0_f64), Vec3::Z).expect("plane");
    let result_sec = section(&result, result.solids[0], &plane, &tol).expect("section");
    let has_arc = result_sec.profiles.iter().any(|o| {
        o.outer
            .edges
            .iter()
            .any(|e| matches!(e, SectionEdge::Arc { .. }))
            || o.holes
                .iter()
                .any(|h| h.edges.iter().any(|e| matches!(e, SectionEdge::Arc { .. })))
    });
    assert!(
        has_arc,
        "a horizontal section through a circular void must yield an Arc loop"
    );
}

// ── regression: arc straddling a rectangle edge (0/2π winding seam) ──────────
//
// A round tool whose centre sits on (or near) a rectangle edge so that the
// surviving circle∩rect arc crosses the angular `0`/`2π` seam of the circle's
// parameterisation. Earlier this produced a non-watertight result
// (`MissingSibling` on the rim arcs) because the cylinder-wall and cap arcs were
// parameterised with inconsistent `0`/`2π` normalisation. All four placements
// (short-edge midpoint, long edge, corner, off-centre) must be watertight.

/// Numerically integrate the area of `{(x,y): in rect ∧ in disc}` for a disc of
/// centre `(cx, cy)` radius `r` clipped to `x∈[x0,x1] × y∈[y0,y1]`. A dense grid
/// is independent of the kernel, so it is a trustworthy oracle for the expected
/// removed volume.
fn disc_rect_overlap_area(cx: f64, cy: f64, r: f64, x0: f64, x1: f64, y0: f64, y1: f64) -> f64 {
    let n = 4000_usize;
    let dx = (x1 - x0) / n as f64;
    let dy = (y1 - y0) / n as f64;
    let cell = dx * dy;
    let r2 = r * r;
    let mut area = 0.0_f64;
    for i in 0..n {
        let x = x0 + (i as f64 + 0.5) * dx;
        let ddx = x - cx;
        for j in 0..n {
            let y = y0 + (j as f64 + 0.5) * dy;
            let ddy = y - cy;
            if ddx * ddx + ddy * ddy <= r2 {
                area += cell;
            }
        }
    }
    area
}

/// `base − tool` where `base` is a rectangular strip extruded along `+z` with the
/// world footprint `x∈[bx0,bx1] × y∈[by0,by1]` and `tool` is a cylinder of radius
/// `r` centred at world `(cx, cy)`, extruded along `+z` through the full height.
///
/// The shared prism frame's in-plane basis is `plane_basis(+z)`, which is a
/// deterministic 90°-rotated basis (not the identity), so the rectangle's
/// half-extents are supplied here in world terms and mapped onto the profile's
/// `(half_w, half_h)` accordingly — `half_w` spans world-Y and `half_h` spans
/// world-X. Building from a world AABB keeps the test readable and independent of
/// that internal basis choice.
///
/// Asserts the result is a single watertight solid whose volume equals the strip
/// volume minus the removed disc∩strip column.
#[allow(clippy::too_many_arguments)]
fn assert_strip_minus_cylinder(bx0: f64, bx1: f64, by0: f64, by1: f64, cx: f64, cy: f64, r: f64) {
    let tol = Tol::default();
    let h = 1.0_f64;
    // half_w → world-Y extent, half_h → world-X extent (rotated frame basis).
    let half_w = 0.5 * (by1 - by0);
    let half_h = 0.5 * (bx1 - bx0);
    let base = ExtrudeLeaf {
        profile: Profile2d::rect(half_w, half_h).expect("rect"),
        origin: Point3::new(0.5 * (bx0 + bx1), 0.5 * (by0 + by1), 0.0_f64),
        axis: Vec3::Z,
        length: h,
    };
    let tool = ExtrudeLeaf {
        profile: Profile2d::circle(r).expect("circle"),
        origin: Point3::new(cx, cy, -0.5_f64),
        axis: Vec3::Z,
        length: h + 1.0_f64,
    };
    let result = prismatic::difference(&base, &tool, &tol).expect("strip minus cylinder");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight strip with edge-straddling circular cut");
    assert_eq!(
        result.solids.len(),
        1usize,
        "the cut must not split the strip"
    );
    let removed = disc_rect_overlap_area(cx, cy, r, bx0, bx1, by0, by1) * h;
    let strip_vol = (bx1 - bx0) * (by1 - by0) * h;
    let expected = strip_vol - removed;
    // The oracle is a grid quadrature; allow a coarse tolerance for it.
    let area_eps = 5e-4_f64;
    assert!(
        (result.signed_volume() - expected).abs() <= area_eps,
        "volume = {} (expected ≈ {expected})",
        result.signed_volume()
    );
}

#[test]
fn strip_minus_cylinder_short_edge_midpoint() {
    // The exact reported failure: circle centred on the midpoint of the strip's
    // short edge (here the bottom edge y = 0), radius > the half-width, so the
    // surviving arc crosses both side edges and the circle's 0/2π angular seam.
    // Before the seam-normalisation fix this raised an `InvalidResult` with
    // `MissingSibling { boundary: [-0.608…, 0.608…] }` paired against
    // `[6.891…, 5.674…]` — the two seam-straddling halves wound up with
    // inconsistent 0/2π parameterisations.
    assert_strip_minus_cylinder(13.8, 14.2, 0.0, 6.0, 14.0, 0.0, 0.35);
}

#[test]
fn strip_minus_cylinder_short_edge_top() {
    // The mirror placement on the opposite short edge (y = 6).
    assert_strip_minus_cylinder(13.8, 14.2, 0.0, 6.0, 14.0, 6.0, 0.35);
}

#[test]
fn strip_minus_cylinder_long_edge() {
    // Circle straddling one long edge (x = 14.2) at mid-height.
    assert_strip_minus_cylinder(13.8, 14.2, 0.0, 6.0, 14.2, 3.0, 0.15);
}

#[test]
fn strip_minus_cylinder_corner() {
    // Circle centred on a corner of the strip — the arc spans a quarter-turn that
    // can also land across the seam depending on the corner.
    assert_strip_minus_cylinder(13.8, 14.2, 0.0, 6.0, 14.2, 6.0, 0.15);
}

#[test]
fn strip_minus_cylinder_short_edge_offset() {
    // Circle on the short edge but off-centre, so the two side-arc pieces have
    // unequal sweep — exercises asymmetric seam crossing.
    assert_strip_minus_cylinder(13.8, 14.2, 0.0, 6.0, 13.95, 0.0, 0.3);
}

#[test]
fn two_circles_with_different_axes_no_common_direction() {
    let tol = Tol::default();
    let cz = ExtrudeLeaf {
        profile: Profile2d::circle(0.1_f64).expect("circle"),
        origin: Point3::new(0.0_f64, 0.0_f64, -1.0_f64),
        axis: Vec3::Z,
        length: 2.0_f64,
    };
    let cx = ExtrudeLeaf {
        profile: Profile2d::circle(0.1_f64).expect("circle"),
        origin: Point3::new(-1.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::X,
        length: 2.0_f64,
    };
    assert!(matches!(
        prismatic::difference(&cz, &cx, &tol),
        Err(PrismError::NoCommonDirection)
    ));
}
