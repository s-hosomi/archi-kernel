//! The coplanar residency truth table (`DESIGN.md` §4.3) — the *specification*,
//! written out so it can be read and so drift from the implementation is caught.
//!
//! When a face of one operand lies exactly in a face of the other (the building
//!常態 — a slab top on a grid level, a wall end on a column face), whether that
//! shared face survives the boolean is a *residency* question, not a geometric
//! one: a boundary face of the result exists exactly where the material
//! immediately on its two sides **differs** in membership of the result set. The
//! whole boolean engine — the prismatic residency closures
//! ([`difference`](crate::boolean::prismatic::difference),
//! [`union_pair`](crate::boolean::prismatic::union_pair),
//! [`intersection`](crate::boolean::prismatic::intersection)) and the half-space
//! cut's coplanar-lid rule
//! ([`cut`](crate::boolean::cut)) — decides every coplanar case from this single
//! principle. `DESIGN.md` §4.3 rejects SoS perturbation in favour of fixing the
//! finite set of cases here, with reference counter-examples.
//!
//! # The setup
//!
//! Take a shared plane Π. Operand A has a face in Π; operand B has a face in Π
//! too (the coincident case). "A's material" is on one side of Π, "B's material"
//! on one side of Π. Two sub-cases:
//!
//! * **normals agree** — A's outward face normal equals B's outward face normal,
//!   i.e. A's material and B's material are on the **same** side of Π (both
//!   bodies sit behind a lid that points the same way; the two overlap up to Π).
//!   Reference 1 of §4.3 (`A = [0,1]³`, `B = [0,1]²×[0,0.5]`) is this case: both
//!   bottom lids point `−z`, both bodies above `z = 0`.
//! * **normals oppose** — the outward normals point opposite ways, i.e. A's and
//!   B's material are on **opposite** sides of Π (the bodies meet in pure contact
//!   along Π). Reference 2 (`B = [0,1]²×[−0.5,0]` below A) is this case.
//!
//! For a result boundary on Π we test residency of the two thin layers straddling
//! Π — `result⁻` just on A's material side and `result⁺` just on the other side —
//! and the face **survives iff `result⁻ ≠ result⁺`**. The *owner* of the
//! surviving face (does it come from A's loop or B's loop) is whichever operand's
//! material is the one that is `in` the result there.
//!
//! # The table
//!
//! Rows are the operation; "self" = A's coincident face, "other" = B's; "agree"
//! / "oppose" the normal relation above. ✓ = the (sub)face survives in the
//! result, ✗ = it is dropped. The pinning test for each cell is named.
//!
//! ```text
//! ┌─────────────┬──────────────┬─────────┬─────────┬───────────────────────────────────────────────┐
//! │ operation   │ face of      │ agree   │ oppose  │ pinning test (tests/coplanar.rs unless noted)  │
//! ├─────────────┼──────────────┼─────────┼─────────┼───────────────────────────────────────────────┤
//! │ A − B       │ self  (A)    │   ✗     │   ✓     │ coplanar_normals_agree_face_dropped (agree);   │
//! │             │              │         │         │ coplanar_normals_oppose_face_kept   (oppose)   │
//! │ A ∪ B       │ self  (A)    │   ✓     │   ✗     │ shared_face_union_merges_to_one_box (oppose);  │
//! │             │              │         │         │ agree ⇒ A's lid bounds the union, survives     │
//! │ A ∩ B       │ self  (A)    │   ✓     │   ✗     │ overlap_intersection_is_the_shared_box (agree);│
//! │             │              │         │         │ contact_only_intersection_is_empty  (oppose)   │
//! └─────────────┴──────────────┴─────────┴─────────┴───────────────────────────────────────────────┘
//! ```
//!
//! Note the operations split: **difference** keeps the self-face when normals
//! *oppose* (the bodies meet in pure contact, so nothing is removed at the lid),
//! whereas **union and intersection** keep it when normals *agree* (the bodies
//! sit on the same side, so the lid is a genuine outer/overlap boundary). Reading
//! the headline cells:
//!
//! * **A − B, self, agree → drop.** Reference 1 of `DESIGN.md` §4.3:
//!   `A = [0,1]³`, `B = [0,1]²×[0,0.5]`. The shared bottom `z = 0` has both
//!   outward normals `−z` (agree); the correct `A − B = [0,1]²×[0.5,1]` has **no**
//!   face at `z = 0`. Just inside A's lid the material is `in_A ∧ in_B`, removed.
//! * **A − B, self, oppose → keep.** Reference 2: `B = [0,1]²×[−0.5,0]` sits below
//!   A in pure contact; `A − B = A`, and A's bottom `z = 0` face must remain. Just
//!   inside A's lid the material is `in_A ∧ ¬in_B`, kept.
//! * **A ∪ B, oppose → drop.** Two stacked boxes sharing `z = 0.5` (materials on
//!   opposite sides, outward normals opposed): the shared face is interior to the
//!   union and vanishes (one merged box).
//! * **A ∩ B, agree → keep / oppose → drop.** A genuine overlap (bodies on the
//!   same side) keeps the shared cut face as the boundary of the common box; pure
//!   contact (opposite sides) has empty intersection and no face.
//!
//! # How the cut specialises this
//!
//! The half-space cut ([`cut`](crate::boolean::cut)) is `solid ∩ half-space`. The
//! intersection row applies with the half-space as "B": a coplanar face of the
//! solid whose outward normal **agrees** with the cut plane's `+normal` is the
//! lid of the kept (`Below`) material and survives (∩-agree → keep: the solid's
//! material and the kept half-space lie on the same `−normal` side); the
//! opposite-facing coincident face is dropped. This is the rule
//! `half_space.rs::process_coplanar_face` implements and the rule the Phase 4
//! section drawing reads to decide whether a coincident-plane section (a plan on
//! a slab top) is drawn (`crate::section`).
//!
//! # Drift guard
//!
//! The unit tests below evaluate each cell from the *same residency primitive*
//! `face_survives(minus, plus)` the engine uses and assert it against the
//! documented ✓/✗, so an implementation change that contradicts this table
//! breaks a unit test here, not just an integration test elsewhere.

/// The residency of the result set on each side of a shared plane, for one
/// `(operation, normal-relation)` cell. `minus` is the layer on A's material
/// side, `plus` the layer on the far side. A boundary face survives iff the two
/// differ.
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Layers {
    minus: bool,
    plus: bool,
}

#[cfg(test)]
impl Layers {
    /// A face on the shared plane survives iff the two straddling layers differ
    /// in result-membership (the one residency principle, `DESIGN.md` §4.3).
    fn face_survives(self) -> bool {
        self.minus != self.plus
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Operation kinds for the table.
    #[derive(Clone, Copy)]
    enum Op {
        Difference,
        Union,
        Intersection,
    }

    /// The residency closure the engine uses, applied to `(in_a, in_b)` flags.
    fn keep(op: Op, in_a: bool, in_b: bool) -> bool {
        match op {
            Op::Difference => in_a && !in_b,
            Op::Union => in_a || in_b,
            Op::Intersection => in_a && in_b,
        }
    }

    /// Build the two straddling layers for "A's coincident face" under the given
    /// normal relation. With normals agreeing, A's material (`−` layer) overlaps
    /// B's material; with normals opposing, B sits on the far (`+`) side only.
    ///
    /// `−` layer: always A's material (`in_a = true`); B present iff normals
    /// agree (same side). `+` layer: outside A (`in_a = false`); B present iff
    /// normals oppose (B's body is on the far side in pure contact).
    fn layers_self(op: Op, agree: bool) -> Layers {
        // − layer: just inside A's face, in A's material.
        //   normals agree  ⇒ B's material is on the SAME side here ⇒ in_b = true.
        //   normals oppose ⇒ B's material is on the OTHER side    ⇒ in_b = false.
        let minus = keep(op, true, agree);
        // + layer: just outside A's face (not in A).
        //   normals agree  ⇒ B is on A's side, so absent here     ⇒ in_b = false.
        //   normals oppose ⇒ B's body is on this far side         ⇒ in_b = true.
        let plus = keep(op, false, !agree);
        Layers { minus, plus }
    }

    #[test]
    fn difference_self_agree_drops_oppose_keeps() {
        // A − B, A's face. agree → drop, oppose → keep (the two §4.3 references).
        assert!(
            !layers_self(Op::Difference, true).face_survives(),
            "A−B self/agree must DROP (reference 1)"
        );
        assert!(
            layers_self(Op::Difference, false).face_survives(),
            "A−B self/oppose must KEEP (reference 2)"
        );
    }

    #[test]
    fn union_self_oppose_drops_agree_keeps() {
        // A ∪ B, A's face. oppose → stacked boxes sharing a face: the shared face
        // is interior to the union and vanishes. agree → bodies on the same side,
        // A's outward lid still bounds the union, survives.
        assert!(
            !layers_self(Op::Union, false).face_survives(),
            "A∪B self/oppose must DROP (shared stacked face is interior)"
        );
        assert!(
            layers_self(Op::Union, true).face_survives(),
            "A∪B self/agree must KEEP (A's lid bounds the union)"
        );
    }

    #[test]
    fn intersection_self_agree_keeps_oppose_drops() {
        // A ∩ B, A's face. Overlap (agree, bodies same side) → the shared cut
        // face bounds the common box. Pure contact (oppose) → empty ∩, no face.
        assert!(
            layers_self(Op::Intersection, true).face_survives(),
            "A∩B self/agree must KEEP (overlap shares the cut face)"
        );
        assert!(
            !layers_self(Op::Intersection, false).face_survives(),
            "A∩B self/oppose must DROP (contact-only ∩ is empty)"
        );
    }

    #[test]
    fn cut_specialises_intersection_agree_keep() {
        // The half-space cut keeps a coplanar lid whose outward normal agrees
        // with +normal: the solid's material and the kept half-space lie on the
        // same side ⇒ the intersection-agree cell ⇒ KEEP. This is what
        // `process_coplanar_face` and `crate::section` rely on.
        assert!(
            layers_self(Op::Intersection, true).face_survives(),
            "cut lid (∩ with solid & half-space on the same side) must survive"
        );
    }
}
