#![allow(clippy::needless_range_loop)]

use approx::assert_relative_eq;
use treecode::vector::*;

#[test]
fn test_vector_zero() {
    let v = vector_zero();
    assert_eq!(v, [0.0, 0.0, 0.0]);
}

#[test]
fn test_vector_unit() {
    let v0 = vector_unit(0);
    assert_eq!(v0, [1.0, 0.0, 0.0]);
    let v1 = vector_unit(1);
    assert_eq!(v1, [0.0, 1.0, 0.0]);
    let v2 = vector_unit(2);
    assert_eq!(v2, [0.0, 0.0, 1.0]);
    let v3 = vector_unit(3);
    assert_eq!(v3, [0.0, 0.0, 0.0]);
}

#[test]
fn test_vector_add() {
    let u = [1.0, 2.0, 3.0];
    let w = [4.0, 5.0, 6.0];
    assert_eq!(vector_add(&u, &w), [5.0, 7.0, 9.0]);
}

#[test]
fn test_vector_sub() {
    let u = [4.0, 5.0, 6.0];
    let w = [1.0, 2.0, 3.0];
    assert_eq!(vector_sub(&u, &w), [3.0, 3.0, 3.0]);
}

#[test]
fn test_vector_mul_scalar() {
    let u = [1.0, 2.0, 3.0];
    assert_eq!(vector_mul_scalar(&u, 2.0), [2.0, 4.0, 6.0]);
    assert_eq!(vector_mul_scalar(&u, 0.0), [0.0, 0.0, 0.0]);
}

#[test]
fn test_vector_div_scalar() {
    let u = [2.0, 4.0, 6.0];
    assert_eq!(vector_div_scalar(&u, 2.0), [1.0, 2.0, 3.0]);
}

#[test]
fn test_vector_dot() {
    let u = [1.0, 2.0, 3.0];
    let v = [4.0, 5.0, 6.0];
    assert_relative_eq!(vector_dot(&u, &v), 32.0);
    assert_relative_eq!(vector_dot(&u, &u), 14.0);
}

#[test]
fn test_vector_length() {
    let v = [3.0, 4.0, 0.0];
    assert_relative_eq!(vector_length(&v), 5.0);
    let v = [0.0, 0.0, 0.0];
    assert_relative_eq!(vector_length(&v), 0.0);
}

#[test]
fn test_vector_distance() {
    let u = [1.0, 2.0, 3.0];
    let v = [4.0, 6.0, 3.0];
    assert_relative_eq!(vector_distance(&u, &v), 5.0);
}

#[test]
fn test_vector_cross() {
    let u = [1.0, 0.0, 0.0];
    let v = [0.0, 1.0, 0.0];
    assert_eq!(vector_cross(&u, &v), [0.0, 0.0, 1.0]);
    let u = [0.0, 1.0, 0.0];
    let v = [1.0, 0.0, 0.0];
    assert_eq!(vector_cross(&u, &v), [0.0, 0.0, -1.0]);
}

#[test]
fn test_matrix_zero() {
    let m = matrix_zero();
    for row in &m {
        for val in row {
            assert_relative_eq!(*val, 0.0);
        }
    }
}

#[test]
fn test_matrix_identity() {
    let m = matrix_identity();
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
    let q = [[1.0; NDIM]; NDIM];
    let r = [[2.0; NDIM]; NDIM];
    let p = matrix_add(&q, &r);
    for row in &p {
        for val in row {
            assert_relative_eq!(*val, 3.0);
        }
    }
}

#[test]
fn test_matrix_sub() {
    let q = [[3.0; NDIM]; NDIM];
    let r = [[1.0; NDIM]; NDIM];
    let p = matrix_sub(&q, &r);
    for row in &p {
        for val in row {
            assert_relative_eq!(*val, 2.0);
        }
    }
}

#[test]
fn test_matrix_mul() {
    let q = matrix_identity();
    let r = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    let p = matrix_mul(&q, &r);
    assert_eq!(p, r);
}

#[test]
fn test_matrix_transpose() {
    let q = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    let t = matrix_transpose(&q);
    assert_eq!(t[0], [1.0, 4.0, 7.0]);
    assert_eq!(t[1], [2.0, 5.0, 8.0]);
    assert_eq!(t[2], [3.0, 6.0, 9.0]);
}

#[test]
fn test_matrix_mul_scalar() {
    let q = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    let p = matrix_mul_scalar(&q, 2.0);
    assert_eq!(p[0], [2.0, 4.0, 6.0]);
    assert_eq!(p[1], [8.0, 10.0, 12.0]);
    assert_eq!(p[2], [14.0, 16.0, 18.0]);
}

#[test]
fn test_matrix_mul_vec() {
    let m = matrix_identity();
    let u = [1.0, 2.0, 3.0];
    let v = matrix_mul_vec(&m, &u);
    assert_eq!(v, u);
}

#[test]
fn test_outer_product() {
    let u = [1.0, 2.0, 3.0];
    let v = [4.0, 5.0, 6.0];
    let p = outer_product(&u, &v);
    assert_eq!(p[0], [4.0, 5.0, 6.0]);
    assert_eq!(p[1], [8.0, 10.0, 12.0]);
    assert_eq!(p[2], [12.0, 15.0, 18.0]);
}

#[test]
fn test_matrix_trace() {
    let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    assert_relative_eq!(matrix_trace(&m), 15.0);
}

#[test]
fn test_dot_sub() {
    let u = [4.0, 5.0, 6.0];
    let w = [1.0, 2.0, 3.0];
    let (s, v) = dot_sub(&u, &w);
    assert_eq!(v, [3.0, 3.0, 3.0]);
    assert_relative_eq!(s, 27.0);
}

#[test]
fn test_dot_mul_mat() {
    let p = matrix_identity();
    let u = [1.0, 2.0, 3.0];
    let (s, v) = dot_mul_mat(&p, &u);
    assert_eq!(v, u);
    assert_relative_eq!(s, 14.0);
}

#[test]
fn test_add_mul_scalar() {
    let mut v = [1.0, 2.0, 3.0];
    let u = [4.0, 5.0, 6.0];
    add_mul_scalar(&mut v, &u, 2.0);
    assert_eq!(v, [9.0, 12.0, 15.0]);
}

#[test]
fn test_add_mul_scalar2() {
    let mut v = [1.0, 2.0, 3.0];
    let u = [4.0, 5.0, 6.0];
    let w = [1.0, 1.0, 1.0];
    add_mul_scalar2(&mut v, &u, 2.0, &w, 3.0);
    assert_eq!(v, [12.0, 15.0, 18.0]);
}

#[test]
fn test_vector_set_scalar() {
    let mut v = [1.0, 2.0, 3.0];
    vector_set_scalar(&mut v, 5.0);
    assert_eq!(v, [5.0, 5.0, 5.0]);
}

#[test]
fn test_vector_add_scalar() {
    let mut v = [1.0, 2.0, 3.0];
    let u = [4.0, 5.0, 6.0];
    vector_add_scalar(&mut v, &u, 10.0);
    assert_eq!(v, [14.0, 15.0, 16.0]);
}

#[test]
fn test_matrix_set_scalar() {
    let mut m = [[0.0; NDIM]; NDIM];
    matrix_set_scalar(&mut m, 7.0);
    for row in &m {
        for val in row {
            assert_relative_eq!(*val, 7.0);
        }
    }
}
