//! Typed 3-vector and 3-matrix arithmetic, mirroring the C `vectmath.h`
//! macros 1:1 (ADDV, MULVS, DOTVP, MULMV, DOTPMULMV, ...).
//!
//! `Vector` and `Matrix` are `#[repr(C)]` newtypes over fixed-size
//! arrays, so they keep the exact memory layout of the C `vector`/`matrix`
//! types (required by the byte-exact state-file format). They deref to
//! their backing arrays, so existing `v[k]` / `m[i][j]` indexing keeps
//! working, and they implement the usual arithmetic operators. Every C
//! macro also has a C-named method (`add`, `mul_scalar`, `dot`, ...).

#![allow(clippy::needless_range_loop)]
// The C-named methods `add`/`sub`/`mul` must stay 1:1 with vectmath.h even
// though they match the operator trait method names (which are also impl'd).
#![allow(clippy::should_implement_trait)]

use std::ops::{
    Add, AddAssign, Deref, DerefMut, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
};

pub type Real = f32;
pub const NDIM: usize = 3;

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector(pub [Real; NDIM]);

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Matrix(pub [[Real; NDIM]; NDIM]);

// ---------------------------------------------------------------------
// Array interop: ergonomic conversion and comparison with raw arrays.
// ---------------------------------------------------------------------

impl From<[Real; NDIM]> for Vector {
    fn from(v: [Real; NDIM]) -> Self {
        Self(v)
    }
}

impl From<Vector> for [Real; NDIM] {
    fn from(v: Vector) -> Self {
        v.0
    }
}

impl PartialEq<[Real; NDIM]> for Vector {
    fn eq(&self, other: &[Real; NDIM]) -> bool {
        self.0 == *other
    }
}

impl PartialEq<Vector> for [Real; NDIM] {
    fn eq(&self, other: &Vector) -> bool {
        *self == other.0
    }
}

impl From<[[Real; NDIM]; NDIM]> for Matrix {
    fn from(m: [[Real; NDIM]; NDIM]) -> Self {
        Self(m)
    }
}

impl From<Matrix> for [[Real; NDIM]; NDIM] {
    fn from(m: Matrix) -> Self {
        m.0
    }
}

impl PartialEq<[[Real; NDIM]; NDIM]> for Matrix {
    fn eq(&self, other: &[[Real; NDIM]; NDIM]) -> bool {
        self.0 == *other
    }
}

impl PartialEq<Matrix> for [[Real; NDIM]; NDIM] {
    fn eq(&self, other: &Matrix) -> bool {
        *self == other.0
    }
}

// ---------------------------------------------------------------------
// Deref to the backing arrays, so `v[k]` and `m[i][j]` keep working.
// ---------------------------------------------------------------------

impl Deref for Vector {
    type Target = [Real; NDIM];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Vector {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl AsRef<[Real]> for Vector {
    fn as_ref(&self) -> &[Real] {
        &self.0
    }
}

impl AsMut<[Real]> for Vector {
    fn as_mut(&mut self) -> &mut [Real] {
        &mut self.0
    }
}

impl Deref for Matrix {
    type Target = [[Real; NDIM]; NDIM];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Matrix {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<'a> IntoIterator for &'a Vector {
    type Item = &'a Real;
    type IntoIter = std::slice::Iter<'a, Real>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a Matrix {
    type Item = &'a [Real; NDIM];
    type IntoIter = std::slice::Iter<'a, [Real; NDIM]>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

// ---------------------------------------------------------------------
// Vector operators.
// ---------------------------------------------------------------------

impl Add for Vector {
    type Output = Self;

    fn add(self, w: Self) -> Self {
        let mut v = Self::default();
        for i in 0..NDIM {
            v[i] = self[i] + w[i];
        }
        v
    }
}

impl Sub for Vector {
    type Output = Self;

    fn sub(self, w: Self) -> Self {
        let mut v = Self::default();
        for i in 0..NDIM {
            v[i] = self[i] - w[i];
        }
        v
    }
}

impl Neg for Vector {
    type Output = Self;

    fn neg(self) -> Self {
        let mut v = Self::default();
        for i in 0..NDIM {
            v[i] = -self[i];
        }
        v
    }
}

impl Mul<Real> for Vector {
    type Output = Self;

    fn mul(self, s: Real) -> Self {
        let mut v = Self::default();
        for i in 0..NDIM {
            v[i] = self[i] * s;
        }
        v
    }
}

impl Mul<Vector> for Real {
    type Output = Vector;

    fn mul(self, v: Vector) -> Vector {
        v * self
    }
}

impl Div<Real> for Vector {
    type Output = Self;

    fn div(self, s: Real) -> Self {
        let mut v = Self::default();
        for i in 0..NDIM {
            v[i] = self[i] / s;
        }
        v
    }
}

impl AddAssign for Vector {
    fn add_assign(&mut self, w: Self) {
        *self = *self + w;
    }
}

impl SubAssign for Vector {
    fn sub_assign(&mut self, w: Self) {
        *self = *self - w;
    }
}

impl MulAssign<Real> for Vector {
    fn mul_assign(&mut self, s: Real) {
        *self = *self * s;
    }
}

impl DivAssign<Real> for Vector {
    fn div_assign(&mut self, s: Real) {
        *self = *self / s;
    }
}

impl std::iter::Sum for Vector {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |acc, v| acc + v)
    }
}

// ---------------------------------------------------------------------
// Matrix operators.
// ---------------------------------------------------------------------

impl Add for Matrix {
    type Output = Self;

    fn add(self, r: Self) -> Self {
        let mut p = Self::default();
        for i in 0..NDIM {
            for j in 0..NDIM {
                p[i][j] = self[i][j] + r[i][j];
            }
        }
        p
    }
}

impl Sub for Matrix {
    type Output = Self;

    fn sub(self, r: Self) -> Self {
        let mut p = Self::default();
        for i in 0..NDIM {
            for j in 0..NDIM {
                p[i][j] = self[i][j] - r[i][j];
            }
        }
        p
    }
}

impl Mul for Matrix {
    type Output = Self;

    fn mul(self, r: Self) -> Self {
        let mut p = Self::default();
        for i in 0..NDIM {
            for j in 0..NDIM {
                for k in 0..NDIM {
                    p[i][j] += self[i][k] * r[k][j];
                }
            }
        }
        p
    }
}

impl Mul<Vector> for Matrix {
    type Output = Vector;

    fn mul(self, u: Vector) -> Vector {
        let mut v = Vector::default();
        for i in 0..NDIM {
            for j in 0..NDIM {
                v[i] += self[i][j] * u[j];
            }
        }
        v
    }
}

impl Mul<Real> for Matrix {
    type Output = Self;

    fn mul(self, s: Real) -> Self {
        let mut p = Self::default();
        for i in 0..NDIM {
            for j in 0..NDIM {
                p[i][j] = self[i][j] * s;
            }
        }
        p
    }
}

impl AddAssign for Matrix {
    fn add_assign(&mut self, r: Self) {
        *self = *self + r;
    }
}

impl MulAssign<Real> for Matrix {
    fn mul_assign(&mut self, s: Real) {
        *self = *self * s;
    }
}

// ---------------------------------------------------------------------
// Vector methods (1:1 with vectmath.h).
// ---------------------------------------------------------------------

impl Vector {
    /// CLRV: all-zero vector.
    pub const fn zero() -> Self {
        Self([0.0; NDIM])
    }

    /// SETVS(s): vector with every component set to `s`.
    pub const fn ones() -> Self {
        Self([1.0; NDIM])
    }

    /// UNITV(j): unit vector along axis `j`.
    pub const fn unit(j: usize) -> Self {
        let mut v = Self::zero();
        if j < NDIM {
            v.0[j] = 1.0;
        }
        v
    }

    /// ADDV: element-wise sum.
    pub fn add(self, w: Self) -> Self {
        self + w
    }

    /// SUBV: element-wise difference.
    pub fn sub(self, w: Self) -> Self {
        self - w
    }

    /// MULVS: multiply by a scalar.
    pub fn mul_scalar(self, s: Real) -> Self {
        self * s
    }

    /// DIVVS: divide by a scalar.
    pub fn div_scalar(self, s: Real) -> Self {
        self / s
    }

    /// DOTVP: dot product.
    pub fn dot(self, u: Self) -> Real {
        let mut s = 0.0;
        for i in 0..NDIM {
            s += self[i] * u[i];
        }
        s
    }

    /// ABSV: length (`rsqrt` == `sqrt` in single precision).
    pub fn length(self) -> Real {
        self.dot(self).sqrt()
    }

    /// DISTV: distance to another vector.
    pub fn distance(self, v: Self) -> Real {
        let mut tmp = 0.0;
        for i in 0..NDIM {
            let d = self[i] - v[i];
            tmp += d * d;
        }
        tmp.sqrt()
    }

    /// CROSSVP: cross product.
    pub fn cross(self, u: Self) -> Self {
        Self([
            self[1] * u[2] - self[2] * u[1],
            self[2] * u[0] - self[0] * u[2],
            self[0] * u[1] - self[1] * u[0],
        ])
    }

    /// SETVS: set every component to `s`.
    pub fn set_scalar(&mut self, s: Real) {
        for i in 0..NDIM {
            self[i] = s;
        }
    }

    /// ADDVS: replace `self` with `u + s` (per component).
    pub fn add_scalar(&mut self, u: Self, s: Real) {
        for i in 0..NDIM {
            self[i] = u[i] + s;
        }
    }
}

// ---------------------------------------------------------------------
// Matrix methods (1:1 with vectmath.h).
// ---------------------------------------------------------------------

impl Matrix {
    /// CLRM: all-zero matrix.
    pub const fn zero() -> Self {
        Self([[0.0; NDIM]; NDIM])
    }

    /// SETMS(s): matrix with every entry set to `s`.
    pub const fn ones() -> Self {
        Self([[1.0; NDIM]; NDIM])
    }

    /// SETMI: identity matrix.
    pub const fn identity() -> Self {
        let mut m = Self::zero();
        let mut i = 0;
        while i < NDIM {
            m.0[i][i] = 1.0;
            i += 1;
        }
        m
    }

    /// ADDM: element-wise sum.
    pub fn add(self, r: Self) -> Self {
        self + r
    }

    /// SUBM: element-wise difference.
    pub fn sub(self, r: Self) -> Self {
        self - r
    }

    /// MULM: matrix product.
    pub fn mul(self, r: Self) -> Self {
        self * r
    }

    /// TRANM: transpose.
    pub fn transpose(self) -> Self {
        let mut p = Self::default();
        for i in 0..NDIM {
            for j in 0..NDIM {
                p[i][j] = self[j][i];
            }
        }
        p
    }

    /// MULMS: multiply by a scalar.
    pub fn mul_scalar(self, s: Real) -> Self {
        self * s
    }

    /// DIVMS: divide by a scalar.
    pub fn div_scalar(self, s: Real) -> Self {
        let mut p = Self::default();
        for i in 0..NDIM {
            for j in 0..NDIM {
                p[i][j] = self[i][j] / s;
            }
        }
        p
    }

    /// MULMV: matrix times vector.
    pub fn mul_vec(self, u: Vector) -> Vector {
        self * u
    }

    /// TRACEM: trace (sum of diagonal entries).
    pub fn trace(self) -> Real {
        let mut s = 0.0;
        for i in 0..NDIM {
            s += self[i][i];
        }
        s
    }

    /// SETMS: set every entry to `s`.
    pub fn set_scalar(&mut self, s: Real) {
        for i in 0..NDIM {
            for j in 0..NDIM {
                self[i][j] = s;
            }
        }
    }
}

// ---------------------------------------------------------------------
// Free functions (1:1 with vectmath.h macros that take an out-param).
// ---------------------------------------------------------------------

/// CLRV: zero a vector in place.
pub fn vector_zero(v: &mut Vector) {
    *v = Vector::zero();
}

/// ABSV: length of a vector.
pub fn vector_length(v: &Vector) -> Real {
    v.dot(*v).sqrt()
}

/// CLRM: zero a matrix in place.
pub fn matrix_zero(m: &mut Matrix) {
    *m = Matrix::zero();
}

/// SETMI: set a matrix to identity in place.
pub fn matrix_identity(m: &mut Matrix) {
    *m = Matrix::identity();
}

/// OUTVP: outer product `v ⊗ u`.
pub fn outer_product(v: &Vector, u: &Vector) -> Matrix {
    let mut p = Matrix::default();
    for i in 0..NDIM {
        for j in 0..NDIM {
            p[i][j] = v[i] * u[j];
        }
    }
    p
}

/// DOTPSUBV: subtract vectors, form dot product.
/// Returns `(s, v)` with `v = u - w` and `s = v·v`.
pub fn dot_sub(u: &Vector, w: &Vector) -> (Real, Vector) {
    let mut v = Vector::default();
    let mut s = 0.0;
    for i in 0..NDIM {
        v[i] = u[i] - w[i];
        s += v[i] * v[i];
    }
    (s, v)
}

/// DOTPMULMV: multiply matrix by vector, form dot product.
/// Returns `(s, v)` with `v = p·u` and `s = v·u`.
pub fn dot_mul_mat(p: &Matrix, u: &Vector) -> (Real, Vector) {
    let mut v = Vector::default();
    let mut s = 0.0;
    for i in 0..NDIM {
        v[i] = Vector::from(p[i]).dot(*u);
        s += v[i] * u[i];
    }
    (s, v)
}

/// ADDMULVS: `v += u * s`.
pub fn add_mul_scalar(v: &mut Vector, u: &Vector, s: Real) {
    for i in 0..NDIM {
        v[i] += u[i] * s;
    }
}

/// ADDMULVS2: `v += u * s + w * r`.
pub fn add_mul_scalar2(v: &mut Vector, u: &Vector, s: Real, w: &Vector, r: Real) {
    for i in 0..NDIM {
        v[i] += u[i] * s + w[i] * r;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layout_is_array_compatible() {
        assert_eq!(
            std::mem::size_of::<Vector>(),
            std::mem::size_of::<[Real; NDIM]>()
        );
        assert_eq!(
            std::mem::size_of::<Matrix>(),
            std::mem::size_of::<[[Real; NDIM]; NDIM]>()
        );
        assert_eq!(std::mem::size_of::<Matrix>(), 36);
    }

    #[test]
    fn operators_match_free_functions() {
        let u = Vector::from([1.0, 2.0, 3.0]);
        let w = Vector::from([4.0, 5.0, 6.0]);
        assert_eq!(u + w, Vector::from([5.0, 7.0, 9.0]));
        assert_eq!(u - w, Vector::from([-3.0, -3.0, -3.0]));
        assert_eq!(u * 2.0, Vector::from([2.0, 4.0, 6.0]));
        assert_eq!(2.0 * u, Vector::from([2.0, 4.0, 6.0]));
        assert_eq!(u / 2.0, Vector::from([0.5, 1.0, 1.5]));
        assert_eq!(-u, Vector::from([-1.0, -2.0, -3.0]));

        let m = Matrix::from([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(m * u, u);
        assert_eq!(m + m, m * 2.0);
        assert_eq!(m.trace(), 3.0);
    }
}
