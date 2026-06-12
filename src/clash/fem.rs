//! FEM / ST-Bridge adapter entry points (`DESIGN.md` §11,
//! `docs/research/06-domain.md` §1, §6).
//!
//! The ST-Bridge parser itself lives on the `fem-io` side (`DESIGN.md` §11): the
//! kernel is a pure library and does not parse `.stb`. What it *does* own is the
//! geometric entry the adapter feeds into — turning ST-Bridge's "two nodes + a
//! section reference" into a CSG [`Extrude`](crate::csg::CsgNode). This module is
//! that entry: [`member_from_axis`] maps a `(start, end, profile)` triple
//! one-to-one onto an extrusion along the member axis.
//!
//! # Units are the adapter's responsibility
//!
//! ST-Bridge node coordinates are in **millimetres**; the kernel is strictly SI
//! metres (`DESIGN.md` §8). The `× 10⁻³` conversion is **not** done here — it is
//! the adapter's job, applied to the node coordinates and section dimensions
//! *before* they reach the kernel. The doc example below shows the adapter
//! pattern explicitly.

use crate::csg::{CsgNode, Profile2d};
use crate::math::Point3;

/// A reason an axis member could not be built.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum AxisMemberError {
    /// `start` and `end` coincide, so the member has no axis direction / zero
    /// length. ST-Bridge gives every line member two distinct nodes, so this is
    /// malformed input.
    ZeroLength,
}

/// Build a member's CSG [`Extrude`](crate::csg::CsgNode) from its axis endpoints
/// and cross-section — the ST-Bridge "2 nodes + section" mapping
/// (`docs/research/06-domain.md` §1.2).
///
/// The extrusion runs from `start` to `end` along their difference, with length
/// `|end − start|`; `profile` is the cross-section placed at `start`. This is the
/// one-to-one image of a `StbColumn` / `StbGirder` (a node pair plus a section
/// id) in the kernel's CSG vocabulary. The resulting node is ready to wrap in a
/// [`Member`](crate::csg::Member) and insert into a [`Model`](crate::model::Model).
///
/// Coordinates are SI metres — the adapter must already have converted ST-Bridge
/// millimetres (module docs).
///
/// # Errors
///
/// [`AxisMemberError::ZeroLength`] if `start` and `end` coincide.
///
/// # Examples
///
/// Converting an ST-Bridge column (millimetre nodes, a rectangular RC section)
/// into a kernel member — note the `× 1e-3` is applied **here, in the adapter**,
/// never inside the kernel (`DESIGN.md` §8, §11):
///
/// ```
/// use archi_kernel::clash::member_from_axis;
/// use archi_kernel::csg::{Member, Profile2d};
/// use archi_kernel::math::Point3;
/// use archi_kernel::model::Model;
///
/// // Straight from a parsed <StbNode .../> pair, still in millimetres.
/// let bottom_mm = [0.0_f64, 0.0_f64, 0.0_f64];
/// let top_mm = [0.0_f64, 0.0_f64, 4000.0_f64]; // a 4 m storey
/// // RC column section 600 × 600 mm (width × depth).
/// let (width_mm, depth_mm) = (600.0_f64, 600.0_f64);
///
/// // Adapter responsibility: mm → m.
/// const MM: f64 = 1e-3;
/// let start = Point3::new(bottom_mm[0] * MM, bottom_mm[1] * MM, bottom_mm[2] * MM);
/// let end = Point3::new(top_mm[0] * MM, top_mm[1] * MM, top_mm[2] * MM);
/// let profile = Profile2d::rect(width_mm * MM / 2.0, depth_mm * MM / 2.0)
///     .expect("positive section");
///
/// let column = member_from_axis(profile, start, end).expect("distinct nodes");
///
/// // Bulk insertion and one-shot evaluation are on `Model`.
/// let mut model = Model::new();
/// model
///     .insert(archi_kernel::csg::StableId(33), Member::new(column))
///     .expect("fresh id");
/// let breps = model.evaluate_all(&Default::default());
/// assert!(breps[&archi_kernel::csg::StableId(33)].is_ok());
/// ```
pub fn member_from_axis(
    profile: Profile2d,
    start: Point3,
    end: Point3,
) -> Result<CsgNode, AxisMemberError> {
    let axis = end - start;
    let length = axis.norm();
    if axis.try_unit().is_none() || length == 0.0 {
        return Err(AxisMemberError::ZeroLength);
    }
    Ok(CsgNode::Extrude {
        profile,
        origin: start,
        axis,
        length,
    })
}

/// The conventional ST-Bridge millimetre→metre factor, named so an adapter can
/// reference the conversion the kernel expects it to have applied.
///
/// Applying it stays the adapter's responsibility (`DESIGN.md` §8, §11); this
/// constant is documentation in value form, not a conversion the kernel performs.
pub const ST_BRIDGE_MM_TO_M: f64 = 1e-3;
