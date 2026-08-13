use approx::assert_relative_eq;
use treecode::types::{Matrix, Real, Vector, NDIM};
use treecode::vecmath::{add_mul_scalar, add_mul_scalar2, dot_mul_mat, dot_sub, outer_product};

#[test]
fn test_vector_zero() {
    let v = Vector::zero();
    assert_eq!(v, [0.0, 0.0, 0.0]);
}

#[test]
fn test_vector_unit() {
    let v0 = Vector::unit(0);
    assert_eq!(v0, [1.0, 0.0, 0.0]);
    let v1 = Vector::unit(1);
    assert_eq!(v1, [0.0, 1.0, 0.0]);
    let v2 = Vector::unit(2);
    assert_eq!(v2, [0.0, 0.0, 1.0]);
    let v3 = Vector::unit(3);
    assert_eq!(v3, [0.0, 0.0, 0.0]);
}

#[test]
fn test_vector_add() {
    let u = Vector::from([1.0, 2.0, 3.0]);
    let w = Vector::from([4.0, 5.0, 6.0]);
    assert_eq!(u.add(w), [5.0, 7.0, 9.0]);
    assert_eq!(u + w, [5.0, 7.0, 9.0]);
}

#[test]
fn test_vector_sub() {
    let u = Vector::from([4.0, 5.0, 6.0]);
    let w = Vector::from([1.0, 2.0, 3.0]);
    assert_eq!(u.sub(w), [3.0, 3.0, 3.0]);
    assert_eq!(u - w, [3.0, 3.0, 3.0]);
}

#[test]
fn test_vector_mul_scalar() {
    let u = Vector::from([1.0, 2.0, 3.0]);
    assert_eq!(u.mul_scalar(2.0), [2.0, 4.0, 6.0]);
    assert_eq!(u.mul_scalar(0.0), [0.0, 0.0, 0.0]);
    assert_eq!(u * 2.0, [2.0, 4.0, 6.0]);
    assert_eq!(2.0 * u, [2.0, 4.0, 6.0]);
}

#[test]
fn test_vector_div_scalar() {
    let u = Vector::from([2.0, 4.0, 6.0]);
    assert_eq!(u.div_scalar(2.0), [1.0, 2.0, 3.0]);
    assert_eq!(u / 2.0, [1.0, 2.0, 3.0]);
}

#[test]
fn test_vector_dot() {
    let u = Vector::from([1.0, 2.0, 3.0]);
    let v = Vector::from([4.0, 5.0, 6.0]);
    assert_relative_eq!(u.dot(v), 32.0);
    assert_relative_eq!(u.dot(u), 14.0);
}

#[test]
fn test_vector_length() {
    let v = Vector::from([3.0, 4.0, 0.0]);
    assert_relative_eq!(v.length(), 5.0);
    let v = Vector::zero();
    assert_relative_eq!(v.length(), 0.0);
}

#[test]
fn test_vector_distance() {
    let u = Vector::from([1.0, 2.0, 3.0]);
    let v = Vector::from([4.0, 6.0, 3.0]);
    assert_relative_eq!(u.distance(v), 5.0);
}

#[test]
fn test_vector_cross() {
    let u = Vector::from([1.0, 0.0, 0.0]);
    let v = Vector::from([0.0, 1.0, 0.0]);
    assert_eq!(u.cross(v), [0.0, 0.0, 1.0]);
    let u = Vector::from([0.0, 1.0, 0.0]);
    let v = Vector::from([1.0, 0.0, 0.0]);
    assert_eq!(u.cross(v), [0.0, 0.0, -1.0]);
}

#[test]
fn test_vector_set_scalar() {
    let mut v = Vector::from([1.0, 2.0, 3.0]);
    v.set_scalar(5.0);
    assert_eq!(v, [5.0, 5.0, 5.0]);
}

#[test]
fn test_vector_add_scalar() {
    let mut v = Vector::from([1.0, 2.0, 3.0]);
    let u = Vector::from([4.0, 5.0, 6.0]);
    v.add_scalar(u, 10.0);
    assert_eq!(v, [14.0, 15.0, 16.0]);
}

#[test]
fn test_matrix_zero() {
    let m = Matrix::zero();
    for row in &m {
        for val in row {
            assert_relative_eq!(*val, 0.0);
        }
    }
}

#[test]
fn test_matrix_identity() {
    let m = Matrix::identity();
    for i in 0..NDIM {
        for j in 0..NDIM {
            if i == j {
                assert_relative_eq!(m[i][j], 1.0);
            } else {
                assert_relative_eq!(m[i][j], 0.0);
            }
        }
    }
}

#[test]
fn test_matrix_add() {
    let q = Matrix::from([[1.0; NDIM]; NDIM]);
    let r = Matrix::from([[2.0; NDIM]; NDIM]);
    let p = q.add(r);
    for row in &p {
        for val in row {
            assert_relative_eq!(*val, 3.0);
        }
    }
    assert_eq!(q + r, p);
}

#[test]
fn test_matrix_sub() {
    let q = Matrix::from([[3.0; NDIM]; NDIM]);
    let r = Matrix::from([[1.0; NDIM]; NDIM]);
    let p = q.sub(r);
    for row in &p {
        for val in row {
            assert_relative_eq!(*val, 2.0);
        }
    }
}

#[test]
fn test_matrix_mul() {
    let q = Matrix::identity();
    let r = Matrix::from([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    let p = q.mul(r);
    assert_eq!(p, r);
    assert_eq!(q * r, r);
}

#[test]
fn test_matrix_transpose() {
    let q = Matrix::from([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    let t = q.transpose();
    assert_eq!(t[0], [1.0, 4.0, 7.0]);
    assert_eq!(t[1], [2.0, 5.0, 8.0]);
    assert_eq!(t[2], [3.0, 6.0, 9.0]);
}

#[test]
fn test_matrix_mul_scalar() {
    let q = Matrix::from([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    let p = q.mul_scalar(2.0);
    assert_eq!(p[0], [2.0, 4.0, 6.0]);
    assert_eq!(p[1], [8.0, 10.0, 12.0]);
    assert_eq!(p[2], [14.0, 16.0, 18.0]);
}

#[test]
fn test_matrix_div_scalar() {
    let q = Matrix::from([[2.0, 4.0, 6.0], [8.0, 10.0, 12.0], [14.0, 16.0, 18.0]]);
    let p = q.div_scalar(2.0);
    assert_eq!(p[0], [1.0, 2.0, 3.0]);
    assert_eq!(p[1], [4.0, 5.0, 6.0]);
    assert_eq!(p[2], [7.0, 8.0, 9.0]);
}

#[test]
fn test_matrix_mul_vec() {
    let m = Matrix::identity();
    let u = Vector::from([1.0, 2.0, 3.0]);
    let v = m.mul_vec(u);
    assert_eq!(v, u);
    assert_eq!(m * u, u);
}

#[test]
fn test_matrix_set_scalar() {
    let mut m = Matrix::zero();
    m.set_scalar(7.0);
    for row in &m {
        for val in row {
            assert_relative_eq!(*val, 7.0);
        }
    }
}

#[test]
fn test_outer_product() {
    let u = Vector::from([1.0, 2.0, 3.0]);
    let v = Vector::from([4.0, 5.0, 6.0]);
    let p = outer_product(&u, &v);
    assert_eq!(p[0], [4.0, 5.0, 6.0]);
    assert_eq!(p[1], [8.0, 10.0, 12.0]);
    assert_eq!(p[2], [12.0, 15.0, 18.0]);
}

#[test]
fn test_matrix_trace() {
    let m = Matrix::from([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    assert_relative_eq!(m.trace(), 15.0);
}

#[test]
fn test_dot_sub() {
    let u = Vector::from([4.0, 5.0, 6.0]);
    let w = Vector::from([1.0, 2.0, 3.0]);
    let (s, v) = dot_sub(&u, &w);
    assert_eq!(v, [3.0, 3.0, 3.0]);
    assert_relative_eq!(s, 27.0);
}

#[test]
fn test_dot_mul_mat() {
    let p = Matrix::identity();
    let u = Vector::from([1.0, 2.0, 3.0]);
    let (s, v) = dot_mul_mat(&p, &u);
    assert_eq!(v, u);
    assert_relative_eq!(s, 14.0);
}

#[test]
fn test_add_mul_scalar() {
    let mut v = Vector::from([1.0, 2.0, 3.0]);
    let u = Vector::from([4.0, 5.0, 6.0]);
    add_mul_scalar(&mut v, &u, 2.0);
    assert_eq!(v, [9.0, 12.0, 15.0]);
}

#[test]
fn test_add_mul_scalar2() {
    let mut v = Vector::from([1.0, 2.0, 3.0]);
    let u = Vector::from([4.0, 5.0, 6.0]);
    let w = Vector::from([1.0, 1.0, 1.0]);
    add_mul_scalar2(&mut v, &u, 2.0, &w, 3.0);
    assert_eq!(v, [12.0, 15.0, 18.0]);
}

#[test]
fn test_from_into() {
    let a: [Real; NDIM] = [1.0, 2.0, 3.0];
    let v: Vector = a.into();
    let back: [Real; NDIM] = v.into();
    assert_eq!(a, back);

    let m: Matrix = [[1.0; NDIM]; NDIM].into();
    let back: [[Real; NDIM]; NDIM] = m.into();
    assert_eq!([[1.0; NDIM]; NDIM], back);
}

#[test]
fn test_indexing() {
    let mut v = Vector::from([1.0, 2.0, 3.0]);
    v[1] = 9.0;
    assert_eq!(v, [1.0, 9.0, 3.0]);
    assert_eq!(v[0], 1.0);

    let mut m = Matrix::zero();
    m[0][1] = 5.0;
    assert_eq!(m[0][1], 5.0);
}

#[test]
fn test_add_assign() {
    let mut v = Vector::from([1.0, 2.0, 3.0]);
    v += Vector::from([4.0, 5.0, 6.0]);
    assert_eq!(v, [5.0, 7.0, 9.0]);

    let mut m = Matrix::identity();
    m += m;
    assert_eq!(m, Matrix::identity().mul_scalar(2.0));
}

#[test]
fn test_ones() {
    assert_eq!(Vector::ones(), [1.0, 1.0, 1.0]);
    let m = Matrix::ones();
    assert_eq!(m[0], [1.0, 1.0, 1.0]);
    assert_eq!(m[2][2], 1.0);
}

#[test]
fn test_const_constructors() {
    const ZV: Vector = Vector::zero();
    const OV: Vector = Vector::ones();
    const UV: Vector = Vector::unit(1);
    const ZM: Matrix = Matrix::zero();
    const OM: Matrix = Matrix::ones();
    const IM: Matrix = Matrix::identity();

    assert_eq!(ZV, [0.0; 3]);
    assert_eq!(OV, [1.0; 3]);
    assert_eq!(UV, [0.0, 1.0, 0.0]);
    assert_eq!(ZM, [[0.0; 3]; 3]);
    assert_eq!(OM, [[1.0; 3]; 3]);
    assert_eq!(IM, Matrix::identity());
}
