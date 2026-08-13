#![allow(clippy::needless_range_loop)]

pub type Real = f32;
pub const NDIM: usize = 3;

pub type Vector = [Real; NDIM];
pub type Matrix = [[Real; NDIM]; NDIM];

pub fn vector_zero() -> Vector {
    [0.0; NDIM]
}

pub fn vector_unit(j: usize) -> Vector {
    let mut v = [0.0; NDIM];
    if j < NDIM {
        v[j] = 1.0;
    }
    v
}

pub fn vector_add(u: &Vector, w: &Vector) -> Vector {
    let mut v = [0.0; NDIM];
    for i in 0..NDIM {
        v[i] = u[i] + w[i];
    }
    v
}

pub fn vector_sub(u: &Vector, w: &Vector) -> Vector {
    let mut v = [0.0; NDIM];
    for i in 0..NDIM {
        v[i] = u[i] - w[i];
    }
    v
}

pub fn vector_mul_scalar(u: &Vector, s: Real) -> Vector {
    let mut v = [0.0; NDIM];
    for i in 0..NDIM {
        v[i] = u[i] * s;
    }
    v
}

pub fn vector_div_scalar(u: &Vector, s: Real) -> Vector {
    let mut v = [0.0; NDIM];
    for i in 0..NDIM {
        v[i] = u[i] / s;
    }
    v
}

pub fn vector_dot(v: &Vector, u: &Vector) -> Real {
    let mut s = 0.0;
    for i in 0..NDIM {
        s += v[i] * u[i];
    }
    s
}

pub fn vector_length(v: &Vector) -> Real {
    vector_dot(v, v).sqrt()
}

pub fn vector_distance(u: &Vector, v: &Vector) -> Real {
    let mut tmp = 0.0;
    for i in 0..NDIM {
        let d = u[i] - v[i];
        tmp += d * d;
    }
    tmp.sqrt()
}

pub fn vector_cross(v: &Vector, u: &Vector) -> Vector {
    [
        v[1] * u[2] - v[2] * u[1],
        v[2] * u[0] - v[0] * u[2],
        v[0] * u[1] - v[1] * u[0],
    ]
}

pub fn matrix_zero() -> Matrix {
    [[0.0; NDIM]; NDIM]
}

pub fn matrix_identity() -> Matrix {
    let mut m = [[0.0; NDIM]; NDIM];
    for i in 0..NDIM {
        m[i][i] = 1.0;
    }
    m
}

pub fn matrix_add(q: &Matrix, r: &Matrix) -> Matrix {
    let mut p = [[0.0; NDIM]; NDIM];
    for i in 0..NDIM {
        for j in 0..NDIM {
            p[i][j] = q[i][j] + r[i][j];
        }
    }
    p
}

pub fn matrix_sub(q: &Matrix, r: &Matrix) -> Matrix {
    let mut p = [[0.0; NDIM]; NDIM];
    for i in 0..NDIM {
        for j in 0..NDIM {
            p[i][j] = q[i][j] - r[i][j];
        }
    }
    p
}

pub fn matrix_mul(q: &Matrix, r: &Matrix) -> Matrix {
    let mut p = [[0.0; NDIM]; NDIM];
    for i in 0..NDIM {
        for j in 0..NDIM {
            for k in 0..NDIM {
                p[i][j] += q[i][k] * r[k][j];
            }
        }
    }
    p
}

pub fn matrix_transpose(q: &Matrix) -> Matrix {
    let mut p = [[0.0; NDIM]; NDIM];
    for i in 0..NDIM {
        for j in 0..NDIM {
            p[i][j] = q[j][i];
        }
    }
    p
}

pub fn matrix_mul_scalar(q: &Matrix, s: Real) -> Matrix {
    let mut p = [[0.0; NDIM]; NDIM];
    for i in 0..NDIM {
        for j in 0..NDIM {
            p[i][j] = q[i][j] * s;
        }
    }
    p
}

pub fn matrix_mul_vec(p: &Matrix, u: &Vector) -> Vector {
    let mut v = [0.0; NDIM];
    for i in 0..NDIM {
        for j in 0..NDIM {
            v[i] += p[i][j] * u[j];
        }
    }
    v
}

pub fn outer_product(v: &Vector, u: &Vector) -> Matrix {
    let mut p = [[0.0; NDIM]; NDIM];
    for i in 0..NDIM {
        for j in 0..NDIM {
            p[i][j] = v[i] * u[j];
        }
    }
    p
}

pub fn matrix_trace(p: &Matrix) -> Real {
    let mut s = 0.0;
    for i in 0..NDIM {
        s += p[i][i];
    }
    s
}

pub fn dot_sub(u: &Vector, w: &Vector) -> (Real, Vector) {
    let mut v = [0.0; NDIM];
    let mut s = 0.0;
    for i in 0..NDIM {
        v[i] = u[i] - w[i];
        s += v[i] * v[i];
    }
    (s, v)
}

pub fn dot_mul_mat(p: &Matrix, u: &Vector) -> (Real, Vector) {
    let mut v = [0.0; NDIM];
    let mut s = 0.0;
    for i in 0..NDIM {
        v[i] = vector_dot(&p[i], u);
        s += v[i] * u[i];
    }
    (s, v)
}

pub fn add_mul_scalar(v: &mut Vector, u: &Vector, s: Real) {
    for i in 0..NDIM {
        v[i] += u[i] * s;
    }
}

pub fn add_mul_scalar2(v: &mut Vector, u: &Vector, s: Real, w: &Vector, r: Real) {
    for i in 0..NDIM {
        v[i] += u[i] * s + w[i] * r;
    }
}

pub fn vector_set_scalar(v: &mut Vector, s: Real) {
    for i in 0..NDIM {
        v[i] = s;
    }
}

pub fn vector_add_scalar(v: &mut Vector, u: &Vector, s: Real) {
    for i in 0..NDIM {
        v[i] = u[i] + s;
    }
}

pub fn matrix_set_scalar(p: &mut Matrix, s: Real) {
    for i in 0..NDIM {
        for j in 0..NDIM {
            p[i][j] = s;
        }
    }
}
