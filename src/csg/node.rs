//! The CSG node vocabulary.

use crate::csg::ids::{OpeningId, StableId};
use crate::csg::profile::Profile2d;
use crate::math::{Point3, Vec3};

/// A node in a member's CSG tree.
///
/// The vocabulary follows `DESIGN.md` §5.1. Every variant is
/// `#[non_exhaustive]`-friendly via the enum attribute so new operations can be
/// added in a semver-compatible way.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum CsgNode {
    /// Extrude a 2-D profile along an axis for a given length.
    Extrude {
        /// The cross-section.
        profile: Profile2d,
        /// Where the extrusion starts (the bottom-cap centre). Members are
        /// placed in world coordinates so that inter-member booleans
        /// (priority clips, clash checks) see real positions.
        origin: Point3,
        /// The extrusion direction (need not be unit; length is separate).
        axis: Vec3,
        /// The extrusion length in metres.
        length: f64,
    },
    /// The union of several sub-nodes.
    Union(Vec<CsgNode>),
    /// `IfcRelVoidsElement`-equivalent semantic opening subtraction.
    ///
    /// Kept distinct from a general [`Difference`](CsgNode::Difference) so that
    /// formwork- and opening-area computations can complete by walking the tree
    /// (`DESIGN.md` §5.1).
    OpeningSubtraction {
        /// The base solid the openings are cut from.
        base: Box<CsgNode>,
        /// The openings, each with its stable id.
        openings: Vec<(OpeningId, Opening)>,
    },
    /// Priority-based deduction for quantity take-off (column → girder → beam →
    /// wall/slab). A deduction is *not* an opening, so it is not modelled with
    /// [`OpeningSubtraction`](CsgNode::OpeningSubtraction) (`DESIGN.md` §5.1).
    Clip {
        /// The base solid being clipped.
        base: Box<CsgNode>,
        /// The members (by stable id) that clip the base.
        clippers: Vec<StableId>,
        /// The rule deciding which member wins on overlap.
        rule: ClipRule,
    },
    /// A general boolean difference (joints, oblique notches).
    Difference {
        /// The solid to keep.
        positive: Box<CsgNode>,
        /// The solid to subtract.
        negative: Box<CsgNode>,
    },
}

/// A semantic opening (void) cut from a member.
///
/// For now an opening is described by its own CSG sub-tree (typically an
/// extrusion). The semantic distinction from a plain difference is preserved by
/// [`CsgNode::OpeningSubtraction`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Opening {
    /// The shape of the void.
    pub shape: CsgNode,
}

/// The rule deciding which member wins where clippers overlap the base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum ClipRule {
    /// The clipper with the higher priority wins (column over girder, …).
    Priority,
}
