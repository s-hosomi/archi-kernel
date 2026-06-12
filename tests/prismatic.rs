//! Prismatic (2.5-D) boolean tests (`DESIGN.md` §10 Phase 3b milestone).
//!
//! Covers the milestone (a rectangular beam notched by a crossing column, with
//! a section), blind and through cuts, opening subtraction (separate, fused,
//! edge), solid splitting, an oblique common direction, the unsupported / limit
//! fail-safes, and the `V(A∪B) = V(A)+V(B)−V(A∩B)` volume identity. Every
//! literal carries an `f64` annotation and an explicit tolerance (`DESIGN.md`
//! §12).

use archi_kernel::boolean::prismatic::{self, ExtrudeLeaf, PrismError};
use archi_kernel::csg::{
    CsgNode, EvalError, Member, Opening, OpeningId, Profile2d, UnsupportedReason,
};
use archi_kernel::math::{Point3, Vec3};
use archi_kernel::primitives::Plane;
use archi_kernel::section::section;
use archi_kernel::tolerance::Tol;
use archi_kernel::topo::ValidateLevel;

const VOL_EPS: f64 = 1e-9;

/// A true axis-aligned box `[x0,x0+wx] × [y0,y0+wy] × [z0,z0+wz]` extruded along
/// `+z`. The `+z` profile axes are `(u, v) = (Y, −X)`, so the profile half-width
/// is the `Y` half and the half-height is the `X` half.
fn box_z(x0: f64, y0: f64, z0: f64, wx: f64, wy: f64, wz: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(wy / 2.0, wx / 2.0).expect("valid rect"),
        origin: Point3::new(x0 + wx / 2.0, y0 + wy / 2.0, z0),
        axis: Vec3::Z,
        length: wz,
    }
}

/// A true axis-aligned box extruded along `+x`. The `+x` profile axes are
/// `(u, v) = (Z, −Y)`, so the profile half-width is the `Z` half and the
/// half-height is the `Y` half.
fn box_x(x0: f64, y0: f64, z0: f64, wx: f64, wy: f64, wz: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(wz / 2.0, wy / 2.0).expect("valid rect"),
        origin: Point3::new(x0, y0 + wy / 2.0, z0 + wz / 2.0),
        axis: Vec3::X,
        length: wx,
    }
}

/// A true axis-aligned box extruded along `+y`. The `+y` profile axes are
/// `(u, v) = (−Z, −X)`, so the profile half-width is the `Z` half and the
/// half-height is the `X` half.
fn box_y(x0: f64, y0: f64, z0: f64, wx: f64, wy: f64, wz: f64) -> ExtrudeLeaf {
    ExtrudeLeaf {
        profile: Profile2d::rect(wz / 2.0, wx / 2.0).expect("valid rect"),
        origin: Point3::new(x0 + wx / 2.0, y0, z0 + wz / 2.0),
        axis: Vec3::Y,
        length: wy,
    }
}

/// A square beam of side `s` centred on the x axis, from `x0`, length `len`,
/// extruded along `+x`.
fn beam_x(x0: f64, len: f64, s: f64) -> ExtrudeLeaf {
    box_x(x0, -s / 2.0, -s / 2.0, len, s, s)
}

// ── milestone: crossing beam − column → notched beam + section ────────────────

#[test]
fn beam_minus_crossing_column_notch_and_section() {
    let tol = Tol::default();
    // Beam: 0.3 × 0.3 square section, x ∈ [0, 2] (centred y, z ∈ [−0.15, 0.15]).
    let beam = beam_x(0.0_f64, 2.0_f64, 0.3_f64);
    // Column along z passing through the beam thickness in y (y ∈ [−0.5, 0.5]
    // covers the beam fully) but only the top of the section in z (z ∈ [0, 1]),
    // over x ∈ [0.9, 1.1]. It removes a top notch, leaving one solid with a
    // U-shaped section at x = 1.
    let column = box_z(0.9_f64, -0.5_f64, 0.0_f64, 0.2_f64, 1.0_f64, 1.0_f64);

    let result = prismatic::difference(&beam, &column, &tol).expect("difference");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight notched beam");
    assert_eq!(result.solids.len(), 1usize, "a top notch leaves one solid");

    // Beam volume 0.3·0.3·2 = 0.18; removed = section y∈[−0.15,0.15] (0.3),
    // z∈[0,0.15] (0.15) over x∈[0.9,1.1] (0.2) = 0.3·0.15·0.2 = 0.009; result
    // = 0.171.
    assert!(
        (result.signed_volume() - 0.171_f64).abs() <= VOL_EPS,
        "volume = {}",
        result.signed_volume()
    );

    // A section across the notch (plane x = 1) yields a section loop.
    let plane = Plane::new(Point3::new(1.0_f64, 0.0_f64, 0.0_f64), Vec3::X).expect("plane");
    let loops = section(&result, result.solids[0], &plane, &tol).expect("section");
    assert!(
        loops.loop_count() >= 1,
        "the notch must show a section loop"
    );
}

#[test]
fn beam_minus_blind_column_is_a_pocket() {
    let tol = Tol::default();
    let beam = beam_x(0.0_f64, 2.0_f64, 0.3_f64);
    // A blind pocket: the cutter reaches neither the far y face nor the bottom
    // z face. y ∈ [−0.05, 0.05] (inside the beam's ±0.15), z ∈ [0, 1] (open at
    // the top only), x ∈ [0.9, 1.1].
    let column = box_z(0.9_f64, -0.05_f64, 0.0_f64, 0.2_f64, 0.1_f64, 1.0_f64);

    let result = prismatic::difference(&beam, &column, &tol).expect("difference");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight pocket");
    assert_eq!(result.solids.len(), 1usize, "a blind pocket does not split");

    // Removed = x∈[0.9,1.1] (0.2), y∈[−0.05,0.05] (0.1), z∈[0,0.15] (0.15)
    // = 0.2·0.1·0.15 = 0.003; result = 0.18 − 0.003 = 0.177.
    assert!(
        (result.signed_volume() - 0.177_f64).abs() <= VOL_EPS,
        "volume = {}",
        result.signed_volume()
    );
}

#[test]
fn full_height_cut_splits_into_two_solids() {
    let tol = Tol::default();
    let beam = beam_x(0.0_f64, 2.0_f64, 0.3_f64);
    // A column covering the beam's whole y,z section over x ∈ [0.9, 1.1]:
    // y ∈ [−0.5, 0.5] and z ∈ [−1, 1] both cover the beam fully, so the band is
    // removed entirely and the beam parts into two pieces.
    let column = box_z(0.9_f64, -0.5_f64, -1.0_f64, 0.2_f64, 1.0_f64, 2.0_f64);

    let result = prismatic::difference(&beam, &column, &tol).expect("difference");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight pieces");
    assert_eq!(
        result.solids.len(),
        2usize,
        "a full-section cut splits the beam"
    );

    // Removed band x∈[0.9,1.1] of the full 0.3×0.3 section (0.09) ⇒ 0.018;
    // result 0.162, in two watertight pieces.
    assert!((result.signed_volume() - 0.162_f64).abs() <= VOL_EPS);
}

// ── opening subtraction (fused) ───────────────────────────────────────────────

/// A wall: 3.0 long (x), 2.0 tall (z), 0.2 thick (y), base corner at the world
/// origin, extruded along `+x`.
fn wall() -> CsgNode {
    leaf_node(box_x(0.0_f64, 0.0_f64, 0.0_f64, 3.0_f64, 0.2_f64, 2.0_f64))
}

/// A window opening punched through the wall thickness (along `y`), spanning
/// `x ∈ [x0, x0+wx]`, `z ∈ [z0, z0+wz]`, and the full wall thickness in `y`.
fn window(x0: f64, z0: f64, wx: f64, wz: f64) -> CsgNode {
    // y ∈ [−0.1, 0.3] comfortably covers the wall's y ∈ [0, 0.2].
    leaf_node(box_y(x0, -0.1_f64, z0, wx, 0.4_f64, wz))
}

/// Wrap an [`ExtrudeLeaf`] as a [`CsgNode::Extrude`].
fn leaf_node(leaf: ExtrudeLeaf) -> CsgNode {
    CsgNode::Extrude {
        profile: leaf.profile,
        origin: leaf.origin,
        axis: leaf.axis,
        length: leaf.length,
    }
}

#[test]
fn wall_with_two_separate_openings_validates() {
    let tol = Tol::default();
    let node = CsgNode::OpeningSubtraction {
        base: Box::new(wall()),
        openings: vec![
            // Two 0.4(x) × 0.4(z) windows, well separated along x.
            (
                OpeningId(1u64),
                Opening {
                    shape: window(0.5_f64, 0.8_f64, 0.4_f64, 0.4_f64),
                },
            ),
            (
                OpeningId(2u64),
                Opening {
                    shape: window(2.1_f64, 0.8_f64, 0.4_f64, 0.4_f64),
                },
            ),
        ],
    };
    let mut member = Member::new(node);
    let brep = member
        .brep(&tol)
        .expect("evaluate wall + 2 windows")
        .clone();
    brep.validate(&tol, ValidateLevel::Full)
        .expect("watertight wall with two through openings");
    assert_eq!(brep.solids.len(), 1usize);

    // Wall volume = 3·0.2·2 = 1.2. Each window removes the full thickness 0.2
    // over 0.4(x) × 0.4(z) = 0.2·0.16 = 0.032; two ⇒ 0.064; result = 1.136.
    assert!(
        (brep.signed_volume() - (1.2_f64 - 0.064_f64)).abs() <= VOL_EPS,
        "volume = {}",
        brep.signed_volume()
    );
}

#[test]
fn overlapping_openings_fuse_into_one_hole() {
    let tol = Tol::default();
    // A single [1.2,1.8]×z window vs two overlapping windows [1.2,1.6] and
    // [1.4,1.8] whose union is exactly [1.2,1.8]: identical removed volume.
    let separate = CsgNode::OpeningSubtraction {
        base: Box::new(wall()),
        openings: vec![(
            OpeningId(1u64),
            Opening {
                shape: window(1.2_f64, 0.8_f64, 0.6_f64, 0.4_f64),
            },
        )],
    };
    let fused = CsgNode::OpeningSubtraction {
        base: Box::new(wall()),
        openings: vec![
            (
                OpeningId(1u64),
                Opening {
                    shape: window(1.2_f64, 0.8_f64, 0.4_f64, 0.4_f64),
                },
            ),
            (
                OpeningId(2u64),
                Opening {
                    shape: window(1.4_f64, 0.8_f64, 0.4_f64, 0.4_f64),
                },
            ),
        ],
    };

    let mut m_sep = Member::new(separate);
    let mut m_fused = Member::new(fused);
    let b_sep = m_sep.brep(&tol).expect("single hole").clone();
    let b_fused = m_fused.brep(&tol).expect("fused hole").clone();

    b_fused
        .validate(&tol, ValidateLevel::Full)
        .expect("fused result watertight");
    assert!(
        (b_sep.signed_volume() - b_fused.signed_volume()).abs() <= VOL_EPS,
        "fused {} vs single {}",
        b_fused.signed_volume(),
        b_sep.signed_volume()
    );
    assert_eq!(
        b_fused.solids.len(),
        1usize,
        "fused opening leaves one solid"
    );
}

#[test]
fn edge_opening_makes_a_notch() {
    let tol = Tol::default();
    // An opening straddling the wall's top edge (z = 2) removes a C-notch rather
    // than a closed hole; the result stays one watertight solid.
    let node = CsgNode::OpeningSubtraction {
        base: Box::new(wall()),
        openings: vec![(
            OpeningId(1u64),
            // x ∈ [1.2,1.8], z ∈ [1.7,2.3] — the top half pokes past z = 2.
            Opening {
                shape: window(1.2_f64, 1.7_f64, 0.6_f64, 0.6_f64),
            },
        )],
    };
    let mut member = Member::new(node);
    let brep = member.brep(&tol).expect("edge opening").clone();
    brep.validate(&tol, ValidateLevel::Full)
        .expect("watertight C-notch");
    assert_eq!(brep.solids.len(), 1usize);

    // Removed: x∈[1.2,1.8] (0.6), z∈[1.7,2.0] (0.3, clipped to the top), full
    // thickness 0.2 ⇒ 0.6·0.3·0.2 = 0.036; result = 1.2 − 0.036 = 1.164.
    assert!(
        (brep.signed_volume() - 1.164_f64).abs() <= VOL_EPS,
        "volume = {}",
        brep.signed_volume()
    );
}

// ── oblique common direction ──────────────────────────────────────────────────

#[test]
fn oblique_common_direction_works() {
    let tol = Tol::default();
    // Beam extruded along d = (1,1,0): a common prismatic direction that is not
    // an axis. A z-column crossing it must still reduce to 2.5-D.
    let dir = Vec3::new(1.0_f64, 1.0_f64, 0.0_f64);
    let beam = ExtrudeLeaf {
        profile: Profile2d::rect(0.15_f64, 0.15_f64).expect("rect"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: dir,
        length: 2.0_f64,
    };
    let column = box_z(0.5_f64, -0.5_f64, -1.0_f64, 0.4_f64, 1.0_f64, 2.0_f64);

    let result = prismatic::difference(&beam, &column, &tol).expect("oblique difference");
    result
        .validate(&tol, ValidateLevel::Full)
        .expect("watertight oblique result");
    assert!(result.signed_volume() > 0.0_f64);
}

// ── unsupported / fail-safe paths ─────────────────────────────────────────────

#[test]
fn orthogonal_h_sections_are_unsupported() {
    let tol = Tol::default();
    let h1 = ExtrudeLeaf {
        profile: Profile2d::h_section(0.1_f64, 0.2_f64, 0.01_f64, 0.02_f64).expect("H"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::X,
        length: 2.0_f64,
    };
    let h2 = ExtrudeLeaf {
        profile: Profile2d::h_section(0.1_f64, 0.2_f64, 0.01_f64, 0.02_f64).expect("H"),
        origin: Point3::new(1.0_f64, 0.0_f64, -1.0_f64),
        axis: Vec3::Z,
        length: 2.0_f64,
    };
    assert!(matches!(
        prismatic::difference(&h1, &h2, &tol),
        Err(PrismError::NoCommonDirection)
    ));
}

#[test]
fn circular_operand_is_unsupported() {
    let tol = Tol::default();
    let beam = beam_x(0.0_f64, 2.0_f64, 0.3_f64);
    let round_column = ExtrudeLeaf {
        profile: Profile2d::circle(0.1_f64).expect("circle"),
        origin: Point3::new(1.0_f64, 0.0_f64, -1.0_f64),
        axis: Vec3::Z,
        length: 2.0_f64,
    };
    assert!(matches!(
        prismatic::difference(&beam, &round_column, &tol),
        Err(PrismError::CircularInvolved { .. })
    ));
}

#[test]
fn unsupported_member_surfaces_machine_readable_reason() {
    let tol = Tol::default();
    let h1 = CsgNode::Extrude {
        profile: Profile2d::h_section(0.1_f64, 0.2_f64, 0.01_f64, 0.02_f64).expect("H"),
        origin: Point3::new(0.0_f64, 0.0_f64, 0.0_f64),
        axis: Vec3::X,
        length: 2.0_f64,
    };
    let h2 = CsgNode::Extrude {
        profile: Profile2d::h_section(0.1_f64, 0.2_f64, 0.01_f64, 0.02_f64).expect("H"),
        origin: Point3::new(1.0_f64, 0.0_f64, -1.0_f64),
        axis: Vec3::Z,
        length: 2.0_f64,
    };
    let node = CsgNode::Difference {
        positive: Box::new(h1),
        negative: Box::new(h2),
    };
    let mut member = Member::new(node);
    match member.brep(&tol) {
        Err(EvalError::Unsupported3dBoolean {
            reason: UnsupportedReason::NoCommonDirection,
        }) => {}
        other => panic!(
            "expected NoCommonDirection reason, got {:?}",
            other.map(|_| ())
        ),
    }
}

#[test]
fn complexity_budget_isolates_a_member() {
    let tol = Tol::default();
    let beam = beam_x(0.0_f64, 2.0_f64, 0.3_f64);
    let column = box_z(0.9_f64, -0.5_f64, -1.0_f64, 0.2_f64, 1.0_f64, 2.0_f64);
    // A budget of 1 is below any non-trivial build's measure.
    match prismatic::difference_with_budget(&beam, &column, &tol, 1usize) {
        Err(PrismError::ComplexityLimit { measure, budget }) => {
            assert_eq!(budget, 1usize);
            assert!(measure > 1usize);
        }
        other => panic!("expected ComplexityLimit, got {other:?}"),
    }
}

// ── volume identity ───────────────────────────────────────────────────────────

#[test]
fn volume_identity_union_intersection() {
    let tol = Tol::default();
    // Two overlapping unit cubes offset by (0.3, 0.3, 0.3).
    let a = box_z(0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 1.0_f64, 1.0_f64);
    let b = box_z(0.3_f64, 0.3_f64, 0.3_f64, 1.0_f64, 1.0_f64, 1.0_f64);

    let v_union = prismatic::union_pair(&a, &b, &tol)
        .expect("union")
        .signed_volume();
    let v_inter = prismatic::intersection(&a, &b, &tol)
        .expect("intersection")
        .signed_volume();
    let v_a = 1.0_f64;
    let v_b = 1.0_f64;

    // V(A∪B) = V(A) + V(B) − V(A∩B).
    assert!(
        (v_union - (v_a + v_b - v_inter)).abs() <= 1e-9_f64,
        "union {} vs {} (= {}+{}−{})",
        v_union,
        v_a + v_b - v_inter,
        v_a,
        v_b,
        v_inter
    );
    // The intersection is the [0.3,1]³ box of volume 0.7³ = 0.343.
    assert!(
        (v_inter - 0.343_f64).abs() <= 1e-9_f64,
        "v_inter = {v_inter}"
    );
}
