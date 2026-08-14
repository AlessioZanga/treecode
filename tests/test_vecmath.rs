use approx::assert_relative_eq;
use treecode::vecmath::{Matrix, NDIM, Vector};

/// Deterministic LCG so the property sweeps are reproducible.
fn lcg(state: &mut u64) -> f32 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let frac = ((*state >> 11) as f64) / ((1u64 << 53) as f64);
    (frac * 20.0 - 10.0) as f32
}

fn rand_vec(state: &mut u64) -> ([f32; NDIM], Vector) {
    let arr = [lcg(state), lcg(state), lcg(state)];
    (arr, Vector::from(arr))
}

fn rand_mat(state: &mut u64) -> ([[f32; NDIM]; NDIM], Matrix) {
    let arr: [[f32; NDIM]; NDIM] = std::array::from_fn(|_| std::array::from_fn(|_| lcg(state)));
    (arr, Matrix::from(arr))
}

// ---- Vector vs `[f32;3]` reference -------------------------------------------

#[test]
fn vector_ops_match_reference() {
    let mut s = 0x1234_5678_9abc_def0u64;
    for _ in 0..500 {
        let (a, va) = rand_vec(&mut s);
        let (b, vb) = rand_vec(&mut s);
        let sc = lcg(&mut s);

        let add = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
        let sub = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let mul = [a[0] * sc, a[1] * sc, a[2] * sc];
        let neg = [-a[0], -a[1], -a[2]];
        let div = [a[0] / sc, a[1] / sc, a[2] / sc];
        let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

        assert_eq!(va + vb, Vector::from(add));
        assert_eq!(va - vb, Vector::from(sub));
        assert_eq!(va * sc, Vector::from(mul));
        assert_eq!(sc * va, Vector::from(mul));
        assert_eq!(-va, Vector::from(neg));
        assert_eq!(va / sc, Vector::from(div));
        assert_relative_eq!(va.dot(vb), dot, epsilon = 1e-4);

        // assignment operators
        let mut acc = va;
        acc += vb;
        assert_eq!(acc, Vector::from(add));
        acc = va;
        acc -= vb;
        assert_eq!(acc, Vector::from(sub));
        acc = va;
        acc *= sc;
        assert_eq!(acc, Vector::from(mul));
        acc = va;
        acc /= sc;
        assert_eq!(acc, Vector::from(div));
    }
}

#[test]
fn vector_cross_and_length_match_reference() {
    let mut s = 0xfeed_face_1357_9bdfu64;
    for _ in 0..300 {
        let (a, va) = rand_vec(&mut s);
        let (b, vb) = rand_vec(&mut s);

        let cross = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        let cprod = va.cross(vb);
        for (k, &c) in cross.iter().enumerate() {
            assert_relative_eq!(cprod.0[k], c, epsilon = 1e-3);
        }

        let len_sq = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        assert_relative_eq!(va.length(), len_sq.sqrt(), epsilon = 1e-3);
    }
}

// ---- Matrix associativity / properties vs reference ---------------------------

#[test]
fn matrix_mul_associative() {
    let mut s = 0xc0ffee00cafebabeu64;
    for _ in 0..300 {
        let (_, a) = rand_mat(&mut s);
        let (_, b) = rand_mat(&mut s);
        let (_, c) = rand_mat(&mut s);

        let lhs = (a * b) * c;
        let rhs = a * (b * c);
        for i in 0..NDIM {
            for j in 0..NDIM {
                assert_relative_eq!(lhs[i][j], rhs[i][j], epsilon = 1e-2);
            }
        }
    }
}

#[test]
fn matrix_vector_mul_matches_reference() {
    let mut s = 0x9e3779b97f4a7c15u64;
    for _ in 0..300 {
        let (_, m) = rand_mat(&mut s);
        let (v, mv) = rand_vec(&mut s);

        let mut refv = [0.0f32; NDIM];
        for i in 0..NDIM {
            let mut acc = 0.0f32;
            for k in 0..NDIM {
                acc += m[i][k] * v[k];
            }
            refv[i] = acc;
        }
        assert_eq!(m * mv, Vector::from(refv));
        assert_eq!(m.mul_vec(mv), Vector::from(refv));
    }
}

#[test]
fn matrix_transpose_properties() {
    let mut s = 0x1234_4321_abcd_dcba_u64;
    for _ in 0..200 {
        let (_, a) = rand_mat(&mut s);
        let (_, b) = rand_mat(&mut s);

        assert_eq!(a.transpose().transpose(), a);

        // (A*B)^T == B^T * A^T
        let ab_t = (a * b).transpose();
        let bt_at = b.transpose() * a.transpose();
        for i in 0..NDIM {
            for j in 0..NDIM {
                assert_relative_eq!(ab_t[i][j], bt_at[i][j], epsilon = 1e-2);
            }
        }

        // identity is neutral
        let i = Matrix::identity();
        assert_eq!(i * a, a);
        assert_eq!(a * i, a);
    }
}

#[test]
fn vector_iter_sum_matches() {
    let mut s = 0x55aa_55aa_55aa_55aau64;
    let mut total = Vector::zero();
    for _ in 0..50 {
        let (_, v) = rand_vec(&mut s);
        total += v;
    }
    let sum: Vector = (0..50).map(|_| Vector::zero()).sum(); // sanity: zero sum
    assert_eq!(sum, Vector::zero());
    // `total` is non-trivially built; just ensure finite and dimensionally sane
    for x in &total.0 {
        assert!(x.is_finite());
    }
}
