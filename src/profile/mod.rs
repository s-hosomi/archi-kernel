//! 2-D cross-section outlines in local profile coordinates.
//!
//! The CSG vocabulary keeps [`Profile2d`](crate::csg::Profile2d) abstract (just
//! the shape parameters and positive-dimension validation). The extrusion
//! builder, however, needs the concrete 2-D outline in the profile's local
//! `(x, y)` plane — centred on the origin and wound counter-clockwise — so that
//! it can be lifted into 3-D along the extrusion axis.
//!
//! This module provides that outline as a [`ProfileGeom`]:
//!
//! * polygonal sections (rectangle, H-section) become an explicit CCW vertex
//!   ring;
//! * the circular section stays an analytic circle — it is *not* approximated
//!   by a polygon, so the extruded side surface is a true [`Cylinder`].
//!
//! All coordinates are in metres (SI), and the winding is CCW as seen looking
//! down the positive local `z` axis (the extrusion direction). With the
//! right-handed basis derived from the axis (`u × v = dir`), a CCW outline
//! extruded along `+dir` yields outward-facing side normals.
//!
//! [`Cylinder`]: crate::primitives::Cylinder

use crate::csg::Profile2d;
use crate::error::KernelError;

/// The concrete 2-D outline of a [`Profile2d`](crate::csg::Profile2d), in local
/// profile coordinates centred on the origin.
///
/// The winding of polygonal outlines is counter-clockwise. The circular outline
/// is kept analytic; no polygonal approximation is made.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum ProfileGeom {
    /// A closed polygon given by its vertices in CCW order. The ring is *not*
    /// repeated (the last vertex connects back to the first implicitly).
    Polygon(Vec<[f64; 2]>),
    /// An analytic circle centred on the local origin.
    Circle {
        /// The radius, in metres.
        radius: f64,
    },
}

impl Profile2d {
    /// Build the 2-D outline of this profile in local profile coordinates.
    ///
    /// The outline is centred on the local origin and, for polygonal sections,
    /// wound counter-clockwise. Returns [`KernelError`] for degenerate
    /// H-sections (a web wider than the flange, or flanges that meet or overlap
    /// through the web).
    ///
    /// # Errors
    ///
    /// * [`KernelError::NonPositiveDimension`] if a derived dimension (e.g. the
    ///   clear web height between the flanges) is not strictly positive.
    pub fn outline(&self) -> Result<ProfileGeom, KernelError> {
        match *self {
            Profile2d::Rect { half_w, half_h } => Ok(rect_outline(half_w, half_h)),
            Profile2d::HSection {
                half_w,
                half_h,
                web,
                flange,
            } => h_section_outline(half_w, half_h, web, flange),
            Profile2d::Circle { radius } => Ok(ProfileGeom::Circle { radius }),
        }
    }
}

/// CCW rectangle centred on the origin.
///
/// Corners in CCW order: `(-w, -h) → (w, -h) → (w, h) → (-w, h)`.
fn rect_outline(half_w: f64, half_h: f64) -> ProfileGeom {
    ProfileGeom::Polygon(vec![
        [-half_w, -half_h],
        [half_w, -half_h],
        [half_w, half_h],
        [-half_w, half_h],
    ])
}

/// CCW 12-vertex H-section centred on the origin.
///
/// The section corresponds to an ST-Bridge H-shape: overall depth `h = 2·half_h`
/// (along local `y`), flange width `b = 2·half_w` (along local `x`), web
/// thickness `tw = web` and flange thickness `tf = flange`.
///
/// The 12 corners, walked counter-clockwise starting at the bottom-right of the
/// bottom flange:
///
/// ```text
///   12┌──────────┐11        y
///     │          │          │
///   1 └──┐    ┌──┘10         └── x
///        │    │
///        │    │  (web)
///        │    │
///   2 ┌──┘    └──┐9
///     │          │
///   3 └──────────┘ ...
/// ```
///
/// Degeneracy checks (`DESIGN.md` §6-1): the web must be strictly narrower
/// than the flange (`tw < b`) and the two flanges must not meet or overlap
/// (`2·tf < h`, i.e. the clear web height is strictly positive).
fn h_section_outline(
    half_w: f64,
    half_h: f64,
    web: f64,
    flange: f64,
) -> Result<ProfileGeom, KernelError> {
    let b = 2.0 * half_w; // overall flange width
    let h = 2.0 * half_h; // overall depth

    // Web must be strictly narrower than the flange; equality produces
    // zero-length edges (degenerate polygon) that break downstream 2-D boolean.
    if web >= b {
        return Err(KernelError::NonPositiveDimension {
            name: "flange_width_minus_web",
            value: b - web,
        });
    }
    // The flanges must leave a strictly positive clear web height between them.
    let clear_web = h - 2.0 * flange;
    if clear_web <= 0.0 {
        return Err(KernelError::NonPositiveDimension {
            name: "clear_web_height",
            value: clear_web,
        });
    }

    let hw = half_w; // = b / 2
    let hh = half_h; // = h / 2
    let tw2 = web / 2.0; // half web thickness
    let yi = hh - flange; // inner face of a flange (|y| of the web-flange junction)

    // CCW from the bottom-left corner of the bottom flange.
    let ring = vec![
        [-hw, -hh],  //  1 bottom flange, bottom-left
        [hw, -hh],   //  2 bottom flange, bottom-right
        [hw, -yi],   //  3 bottom flange, top-right
        [tw2, -yi],  //  4 web, bottom-right
        [tw2, yi],   //  5 web, top-right
        [hw, yi],    //  6 top flange, bottom-right
        [hw, hh],    //  7 top flange, top-right
        [-hw, hh],   //  8 top flange, top-left
        [-hw, yi],   //  9 top flange, bottom-left
        [-tw2, yi],  // 10 web, top-left
        [-tw2, -yi], // 11 web, bottom-left
        [-hw, -yi],  // 12 bottom flange, top-left
    ];
    Ok(ProfileGeom::Polygon(ring))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rect_is_ccw_four_vertices() {
        let p = Profile2d::rect(0.15_f64, 0.3_f64).expect("valid rect");
        match p.outline().expect("outline") {
            ProfileGeom::Polygon(v) => {
                assert_eq!(v.len(), 4usize);
                // Shoelace area positive ⇒ CCW.
                assert!(signed_area(&v) > 0.0_f64);
            }
            ProfileGeom::Circle { .. } => panic!("rect must be a polygon"),
        }
    }

    #[test]
    fn h_section_is_ccw_twelve_vertices() {
        // b = 0.2, h = 0.4, tw = 0.01, tf = 0.02.
        let p = Profile2d::h_section(0.1_f64, 0.2_f64, 0.01_f64, 0.02_f64).expect("valid H");
        match p.outline().expect("outline") {
            ProfileGeom::Polygon(v) => {
                assert_eq!(v.len(), 12usize);
                assert!(signed_area(&v) > 0.0_f64);
            }
            ProfileGeom::Circle { .. } => panic!("H must be a polygon"),
        }
    }

    #[test]
    fn h_section_area_matches_formula() {
        // Area = 2·b·tf + (h − 2·tf)·tw.
        let half_w = 0.1_f64;
        let half_h = 0.2_f64;
        let web = 0.01_f64;
        let flange = 0.02_f64;
        let b = 2.0_f64 * half_w;
        let h = 2.0_f64 * half_h;
        let expected = 2.0_f64 * b * flange + (h - 2.0_f64 * flange) * web;
        let p = Profile2d::h_section(half_w, half_h, web, flange).expect("valid H");
        let area = match p.outline().expect("outline") {
            ProfileGeom::Polygon(v) => signed_area(&v),
            ProfileGeom::Circle { .. } => panic!("H must be a polygon"),
        };
        assert!((area - expected).abs() < 1e-12_f64, "area = {area}");
    }

    #[test]
    fn h_section_rejects_wide_web() {
        // web (0.3) > b (0.2).
        let p = Profile2d::h_section(0.1_f64, 0.2_f64, 0.3_f64, 0.02_f64).expect("ctor ok");
        assert!(matches!(
            p.outline(),
            Err(KernelError::NonPositiveDimension { .. })
        ));
    }

    #[test]
    fn h_section_rejects_thick_flanges() {
        // 2·tf (0.4) ≥ h (0.4).
        let p = Profile2d::h_section(0.1_f64, 0.2_f64, 0.01_f64, 0.2_f64).expect("ctor ok");
        assert!(matches!(
            p.outline(),
            Err(KernelError::NonPositiveDimension { .. })
        ));
    }

    #[test]
    fn circle_stays_analytic() {
        let p = Profile2d::circle(0.3_f64).expect("valid circle");
        assert_eq!(p.outline(), Ok(ProfileGeom::Circle { radius: 0.3_f64 }));
    }

    /// Shoelace signed area, positive for CCW winding.
    fn signed_area(v: &[[f64; 2]]) -> f64 {
        let n = v.len();
        let mut a = 0.0_f64;
        for i in 0..n {
            let p = v[i];
            let q = v[(i + 1) % n];
            a += p[0] * q[1] - q[0] * p[1];
        }
        a / 2.0_f64
    }
}
