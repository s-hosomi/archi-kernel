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
    let loops = section(&result, result.solids[0], &plane, &tol).expect("section");
    let has_arc = loops.outlines.iter().any(|o| {
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
