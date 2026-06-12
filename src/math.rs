//! Minimal 3-D vector mathematics for the archi-kernel.
//!
//! This module intentionally exposes only the operations needed by the kernel
//! itself — no BLAS-level generality. Keeping the surface small prevents the
//! public API from coupling to a third-party linear-algebra library.
//!
//! All computations are in SI metres (lengths) or radians (angles).

use std::ops::{Add, Deref, Mul, Neg, Sub};

// ── Point3 ──────────────────────────────────────────────────────────────────

/// A point in 3-D Euclidean space (coordinates in metres, SI).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Point3 {
    /// x coordinate.
    pub x: f64,
    /// y coordinate.
    pub y: f64,
    /// z coordinate.
    pub z: f64,
}

impl Point3 {
    /// Construct a point from its three coordinates.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The origin `(0, 0, 0)`.
    #[inline]
    pub fn origin() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }
}

// Point3 − Point3 → Vec3
impl Sub for Point3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, rhs: Self) -> Vec3 {
        Vec3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

// Point3 + Vec3 → Point3
impl Add<Vec3> for Point3 {
    type Output = Point3;
    #[inline]
    fn add(self, rhs: Vec3) -> Point3 {
        Point3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

// Point3 − Vec3 → Point3
impl Sub<Vec3> for Point3 {
    type Output = Point3;
    #[inline]
    fn sub(self, rhs: Vec3) -> Point3 {
        Point3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

// ── Vec3 ─────────────────────────────────────────────────────────────────────

/// A free vector in 3-D space.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl Vec3 {
    /// Construct a vector from its three components.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The zero vector.
    pub const ZERO: Vec3 = Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// The unit vector along the x axis.
    pub const X: Vec3 = Vec3 {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };

    /// The unit vector along the y axis.
    pub const Y: Vec3 = Vec3 {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    };

    /// The unit vector along the z axis.
    pub const Z: Vec3 = Vec3 {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };

    /// Dot product.
    #[inline]
    pub fn dot(self, rhs: Vec3) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    /// Cross product.
    #[inline]
    pub fn cross(self, rhs: Vec3) -> Vec3 {
        Vec3::new(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x,
        )
    }

    /// Euclidean norm (length).
    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_squared().sqrt()
    }

    /// Squared Euclidean norm.
    #[inline]
    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    /// Try to produce a unit vector in the same direction.
    ///
    /// Returns `None` if any component is non-finite or if the norm is zero.
    pub fn try_unit(self) -> Option<Unit3> {
        if !self.x.is_finite() || !self.y.is_finite() || !self.z.is_finite() {
            return None;
        }
        let n = self.norm();
        if n == 0.0 {
            return None;
        }
        Some(Unit3(Vec3::new(self.x / n, self.y / n, self.z / n)))
    }
}

// Vec3 + Vec3
impl Add for Vec3 {
    type Output = Vec3;
    #[inline]
    fn add(self, rhs: Vec3) -> Vec3 {
        Vec3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

// Vec3 − Vec3
impl Sub for Vec3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, rhs: Vec3) -> Vec3 {
        Vec3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

// Vec3 * f64
impl Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, s: f64) -> Vec3 {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }
}

// f64 * Vec3
impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}

// −Vec3
impl Neg for Vec3 {
    type Output = Vec3;
    #[inline]
    fn neg(self) -> Vec3 {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}

// ── Unit3 ────────────────────────────────────────────────────────────────────

/// A unit vector in 3-D space.
///
/// The invariant `|v| == 1` is guaranteed by construction: the only public
/// way to obtain a `Unit3` is via [`Vec3::try_unit`]. Internal code may use
/// [`Unit3::new_unchecked`] when the invariant is already established.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Unit3(Vec3);

impl Unit3 {
    /// The unit vector along the x axis.
    pub const X: Unit3 = Unit3(Vec3::X);

    /// The unit vector along the y axis.
    pub const Y: Unit3 = Unit3(Vec3::Y);

    /// The unit vector along the z axis.
    pub const Z: Unit3 = Unit3(Vec3::Z);

    /// Return the underlying [`Vec3`].
    #[inline]
    pub fn as_vec(self) -> Vec3 {
        self.0
    }

    /// Construct a `Unit3` without checking the norm.
    ///
    /// # Safety (logical)
    /// The caller must ensure that `v` is a unit vector. This function is
    /// intended only for internal use where the invariant is already
    /// established by prior computation.
    #[inline]
    pub(crate) fn new_unchecked(v: Vec3) -> Self {
        Self(v)
    }
}

impl Deref for Unit3 {
    type Target = Vec3;
    #[inline]
    fn deref(&self) -> &Vec3 {
        &self.0
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_subtraction_gives_vec() {
        let a = Point3::new(1.0_f64, 2.0_f64, 3.0_f64);
        let b = Point3::new(4.0_f64, 6.0_f64, 8.0_f64);
        let v = b - a;
        assert_eq!(v, Vec3::new(3.0_f64, 4.0_f64, 5.0_f64));
    }

    #[test]
    fn point_plus_vec() {
        let p = Point3::origin();
        let v = Vec3::new(1.0_f64, 2.0_f64, 3.0_f64);
        assert_eq!(p + v, Point3::new(1.0_f64, 2.0_f64, 3.0_f64));
    }

    #[test]
    fn point_minus_vec() {
        let p = Point3::new(1.0_f64, 2.0_f64, 3.0_f64);
        let v = Vec3::new(1.0_f64, 2.0_f64, 3.0_f64);
        assert_eq!(p - v, Point3::origin());
    }

    #[test]
    fn cross_product_orthogonal() {
        let x = Vec3::X;
        let y = Vec3::Y;
        let z = x.cross(y);
        assert!((z - Vec3::Z).norm() < 1e-15_f64);
    }

    #[test]
    fn try_unit_normalises() {
        let v = Vec3::new(3.0_f64, 4.0_f64, 0.0_f64);
        let u = v.try_unit().expect("non-zero vector must succeed");
        assert!((u.norm() - 1.0_f64).abs() < 1e-15_f64);
    }

    #[test]
    fn try_unit_zero_returns_none() {
        assert!(Vec3::ZERO.try_unit().is_none());
    }

    #[test]
    fn try_unit_nan_returns_none() {
        assert!(Vec3::new(f64::NAN, 0.0_f64, 0.0_f64).try_unit().is_none());
    }

    #[test]
    fn unit3_deref_gives_vec() {
        let u = Unit3::Z;
        assert_eq!(u.z, 1.0_f64);
        assert_eq!(u.norm(), 1.0_f64);
    }
}
