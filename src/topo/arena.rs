//! A self-contained generational arena.
//!
//! Topology entities are stored in arenas and referenced by generational
//! handles ([`Id`]) rather than Rust references, which keeps the cyclic
//! half-edge structure free of borrow-checker conflicts (the community-standard
//! solution, see `docs/research/03-topology.md` §6).
//!
//! The "generation" component of each handle detects the classic *dangling*
//! bug: a slot can be removed and later re-used for a different entity; an
//! [`Id`] minted before the removal carries the old generation and so fails to
//! resolve, instead of silently aliasing the new occupant.
//!
//! The arena is held with a concrete element type per entity kind (no trait
//! abstraction) — abstraction layers are a hiding place for bugs and the
//! kernel keeps a zero-dependency, easy-to-verify implementation
//! (`DESIGN.md` §3.3).

use std::marker::PhantomData;

/// A generational handle into an [`Arena`].
///
/// The type parameter `T` ties the handle to the kind of entity it indexes, so
/// a `Id<Vertex>` cannot be used to index an `Arena<Face>`. Handles are cheap
/// `Copy` values carrying an index and a generation; equality and hashing
/// compare both fields (and the phantom type), so a stale handle never compares
/// equal to a live one occupying the same slot.
pub struct Id<T> {
    idx: u32,
    gen: u32,
    _marker: PhantomData<fn() -> T>,
}

impl<T> Id<T> {
    /// The slot index this handle refers to.
    #[inline]
    pub fn index(self) -> u32 {
        self.idx
    }

    /// The generation this handle was minted with.
    #[inline]
    pub fn generation(self) -> u32 {
        self.gen
    }
}

// Manual impls: deriving would add a spurious `T: Trait` bound. Handles are
// always `Copy`/`Eq`/`Hash` regardless of the entity type they index.
impl<T> Clone for Id<T> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Copy for Id<T> {}

impl<T> PartialEq for Id<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.idx == other.idx && self.gen == other.gen
    }
}

impl<T> Eq for Id<T> {}

impl<T> std::hash::Hash for Id<T> {
    #[inline]
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.idx.hash(state);
        self.gen.hash(state);
    }
}

impl<T> std::fmt::Debug for Id<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Id({}, gen {})", self.idx, self.gen)
    }
}

#[cfg(feature = "serde")]
impl<T> serde::Serialize for Id<T> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        (self.idx, self.gen).serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T> serde::Deserialize<'de> for Id<T> {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let (idx, gen) = <(u32, u32)>::deserialize(deserializer)?;
        Ok(Id {
            idx,
            gen,
            _marker: PhantomData,
        })
    }
}

/// One slot of an [`Arena`]: either a live value or a tombstone awaiting reuse.
///
/// The generation is stored on the slot, not the value, so it survives the
/// transition from `Live` to `Dead` and back.
#[derive(Debug, Clone)]
enum Slot<T> {
    Live { gen: u32, value: T },
    Dead { gen: u32 },
}

impl<T> Slot<T> {
    #[inline]
    fn generation(&self) -> u32 {
        match self {
            Slot::Live { gen, .. } | Slot::Dead { gen } => *gen,
        }
    }
}

/// A generational arena storing values of a single concrete type `T`.
///
/// Insertion returns an [`Id`]; resolution succeeds only while the handle's
/// generation matches the slot's, so handles into removed entities are
/// detected rather than silently reused.
#[derive(Debug, Clone)]
pub struct Arena<T> {
    slots: Vec<Slot<T>>,
    free: Vec<u32>,
}

impl<T> Arena<T> {
    /// Create an empty arena.
    pub fn new() -> Self {
        Self {
            slots: Vec::new(),
            free: Vec::new(),
        }
    }

    /// Insert `value`, reusing a freed slot if one is available.
    ///
    /// Returns a handle whose generation matches the chosen slot.
    pub fn insert(&mut self, value: T) -> Id<T> {
        if let Some(idx) = self.free.pop() {
            let gen = self.slots[idx as usize].generation();
            self.slots[idx as usize] = Slot::Live { gen, value };
            Id {
                idx,
                gen,
                _marker: PhantomData,
            }
        } else {
            let idx = self.slots.len() as u32;
            self.slots.push(Slot::Live { gen: 0, value });
            Id {
                idx,
                gen: 0,
                _marker: PhantomData,
            }
        }
    }

    /// Resolve `id` to a shared reference, or `None` if it is stale or out of
    /// range.
    pub fn get(&self, id: Id<T>) -> Option<&T> {
        match self.slots.get(id.idx as usize) {
            Some(Slot::Live { gen, value }) if *gen == id.gen => Some(value),
            _ => None,
        }
    }

    /// Resolve `id` to a mutable reference, or `None` if it is stale or out of
    /// range.
    pub fn get_mut(&mut self, id: Id<T>) -> Option<&mut T> {
        match self.slots.get_mut(id.idx as usize) {
            Some(Slot::Live { gen, value }) if *gen == id.gen => Some(value),
            _ => None,
        }
    }

    /// Remove the value referred to by `id`, returning it if it was live.
    ///
    /// The slot's generation is bumped so that the freed handle (and any copy
    /// of it) no longer resolves.
    ///
    /// If the generation counter would overflow `u32::MAX` (after ~4 billion
    /// remove+insert cycles on the same slot), the slot is **retired** rather
    /// than wrapped back to zero.  A retired slot stays `Dead` forever and is
    /// never returned to the free list, which ensures that the original gen-0
    /// handle cannot silently alias a future occupant.  In practice the
    /// retirement threshold is unreachable in any B-rep workload; the guard is
    /// present for correctness, not performance.
    pub fn remove(&mut self, id: Id<T>) -> Option<T> {
        let slot = self.slots.get_mut(id.idx as usize)?;
        match slot {
            Slot::Live { gen, .. } if *gen == id.gen => {
                // Use checked_add: if the generation would overflow u32::MAX we
                // retire the slot permanently (do not push back to `free`).
                // This prevents the gen-0 aliasing hazard while keeping the
                // public API panic-free.
                let next_gen = gen.checked_add(1);
                let dead_gen = next_gen.unwrap_or(*gen); // stay at MAX when retiring
                let old = std::mem::replace(slot, Slot::Dead { gen: dead_gen });
                if next_gen.is_some() {
                    // Normal case: slot can be recycled.
                    self.free.push(id.idx);
                }
                // If next_gen is None the slot is silently retired; it will
                // never be handed out again, preventing generation wraparound.
                match old {
                    Slot::Live { value, .. } => Some(value),
                    Slot::Dead { .. } => unreachable!(),
                }
            }
            _ => None,
        }
    }

    /// `true` if `id` currently resolves to a live value.
    #[inline]
    pub fn contains(&self, id: Id<T>) -> bool {
        self.get(id).is_some()
    }

    /// Number of live values in the arena.
    pub fn len(&self) -> usize {
        self.slots.len() - self.free.len()
    }

    /// `true` if the arena holds no live values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate over `(Id, &value)` for every live slot.
    pub fn iter(&self) -> impl Iterator<Item = (Id<T>, &T)> {
        self.slots
            .iter()
            .enumerate()
            .filter_map(|(idx, slot)| match slot {
                Slot::Live { gen, value } => Some((
                    Id {
                        idx: idx as u32,
                        gen: *gen,
                        _marker: PhantomData,
                    },
                    value,
                )),
                Slot::Dead { .. } => None,
            })
    }

    /// Iterate over the [`Id`] of every live slot.
    pub fn ids(&self) -> impl Iterator<Item = Id<T>> + '_ {
        self.iter().map(|(id, _)| id)
    }
}

impl<T> Default for Arena<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_and_get_roundtrips() {
        let mut arena: Arena<i32> = Arena::new();
        let a = arena.insert(10_i32);
        let b = arena.insert(20_i32);
        assert_eq!(arena.get(a), Some(&10_i32));
        assert_eq!(arena.get(b), Some(&20_i32));
        assert_eq!(arena.len(), 2);
    }

    #[test]
    fn removed_handle_does_not_resolve() {
        let mut arena: Arena<i32> = Arena::new();
        let a = arena.insert(10_i32);
        assert_eq!(arena.remove(a), Some(10_i32));
        assert_eq!(arena.get(a), None);
        assert!(!arena.contains(a));
    }

    #[test]
    fn stale_handle_detected_after_slot_reuse() {
        let mut arena: Arena<i32> = Arena::new();
        let a = arena.insert(10_i32);
        arena.remove(a);
        // Re-using the freed slot bumps the generation: the new handle differs.
        let b = arena.insert(99_i32);
        assert_eq!(a.index(), b.index());
        assert_ne!(a, b);
        // The stale handle must not alias the new occupant.
        assert_eq!(arena.get(a), None);
        assert_eq!(arena.get(b), Some(&99_i32));
    }

    #[test]
    fn get_mut_mutates_in_place() {
        let mut arena: Arena<i32> = Arena::new();
        let a = arena.insert(1_i32);
        *arena.get_mut(a).expect("live") += 41_i32;
        assert_eq!(arena.get(a), Some(&42_i32));
    }

    #[test]
    fn iter_visits_only_live_slots() {
        let mut arena: Arena<i32> = Arena::new();
        let a = arena.insert(1_i32);
        let _b = arena.insert(2_i32);
        let c = arena.insert(3_i32);
        arena.remove(a);
        let mut vals: Vec<i32> = arena.iter().map(|(_, v)| *v).collect();
        vals.sort_unstable();
        assert_eq!(vals, vec![2_i32, 3_i32]);
        assert!(arena.contains(c));
    }

    /// After ~2^32 remove+insert cycles on the same slot the generation counter
    /// would have wrapped back to 0, making the original gen-0 handle silently
    /// alias the new occupant.  The fix retires the slot when `checked_add`
    /// overflows (i.e. gen == u32::MAX): `remove()` does not push the slot back
    /// onto `free`, so it can never be reused again.
    ///
    /// We fast-forward the slot to `Live { gen: u32::MAX }` to avoid running
    /// ~4 billion real cycles; everything after that uses only the public API.
    #[test]
    fn slot_is_retired_at_generation_max_not_wrapped() {
        let mut arena: Arena<i32> = Arena::new();

        // Obtain the first handle (gen == 0, slot 0).
        let id0 = arena.insert(10_i32);
        assert_eq!(id0.index(), 0u32);
        assert_eq!(id0.generation(), 0u32);

        // Fast-forward: put slot 0 into the state it would have after
        // u32::MAX - 1 remove+insert cycles, i.e. Live { gen: u32::MAX }.
        arena.slots[0] = Slot::Live {
            gen: u32::MAX,
            value: 10_i32,
        };
        // Clear the free list so slot 0 is treated as currently live.
        arena.free.clear();
        let id_max = Id::<i32> {
            idx: 0u32,
            gen: u32::MAX,
            _marker: PhantomData,
        };

        // The 2^32-th remove: without the fix, wrapping_add(1) would produce
        // Dead { gen: 0 } and the slot would go back onto `free`.  With the
        // fix, checked_add overflows → slot is retired (stays Dead { gen:
        // u32::MAX }) and is NOT added to `free`.
        assert_eq!(arena.remove(id_max), Some(10_i32));

        // Slot must be Dead with gen == u32::MAX (retired, not wrapped to 0).
        assert!(
            matches!(arena.slots[0], Slot::Dead { gen: u32::MAX }),
            "retired slot must keep gen == u32::MAX, not wrap to 0: {:?}",
            arena.slots[0]
        );

        // The free list must NOT contain slot 0 — it has been permanently
        // retired and must never be handed out again.
        assert!(
            !arena.free.contains(&0u32),
            "retired slot must not be on the free list"
        );

        // insert() must allocate a fresh slot (index 1) rather than reusing
        // the retired slot 0.
        let id_new = arena.insert(777_i32);
        assert_eq!(
            id_new.index(),
            1u32,
            "insert after retirement must use a new slot, not the retired one"
        );
        assert_eq!(id_new.generation(), 0u32);

        // The original gen-0 handle must still be stale (returns None).
        assert_eq!(
            arena.get(id0),
            None,
            "stale gen-0 handle must not resolve after slot retirement"
        );
    }

    /// Verify the cycle-counting premise independently: each remove+insert on
    /// the same slot bumps the generation by exactly 1, confirming that 2^32
    /// cycles would be needed to cause a wraparound (which is now prevented).
    #[test]
    fn generation_increments_by_one_per_remove_insert_cycle() {
        let mut arena: Arena<i32> = Arena::new();
        let mut id = arena.insert(0_i32);
        assert_eq!(id.index(), 0u32);
        assert_eq!(id.generation(), 0u32);

        for expected_gen in 1u32..=1000u32 {
            arena.remove(id);
            id = arena.insert(0_i32);
            assert_eq!(id.index(), 0u32, "slot is reused every cycle");
            assert_eq!(
                id.generation(),
                expected_gen,
                "generation must increment by exactly 1 per cycle"
            );
        }
    }
}
