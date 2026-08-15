use crate::types::Vector;

/// Stateful PRNG handle. The C code relied on libc's process-global
/// `random()`/`srandom()`. To drop that global (and the `libc` RNG dependency)
/// while staying byte-for-byte identical to the reference binary, this is a
/// pure-Rust reimplementation of glibc's `TYPE_3` additive-feedback generator
/// (the default used by `srandom`/`random`): a 31-degree trinomial with
/// separation 3, seeded by the same LCG and warmed up with the same 310
/// discarded draws. A single `unit()` therefore yields the exact same stream
/// as the C `random()` on any platform.
#[derive(Clone, Copy)]
pub struct RngState {
    state: [i32; 32],
    fptr: usize,
    rptr: usize,
}

const RAND_DEG: usize = 31;
const RAND_SEP: usize = 3;
const END_PTR: usize = 32;

impl RngState {
    pub fn new(seed: u32) -> Self {
        let mut rng = RngState {
            state: [0i32; 32],
            fptr: 1 + RAND_SEP,
            rptr: 1,
        };

        let mut seed = seed as i64;
        if seed == 0 {
            seed = 1;
        }
        rng.state[1] = seed as i32;

        let mut word = seed;
        for i in 1..RAND_DEG {
            let hi = word / 127773;
            let lo = word % 127773;
            word = 16807 * lo - 2836 * hi;
            if word < 0 {
                word += 2147483647;
            }
            rng.state[1 + i] = word as i32;
        }

        rng.fptr = 1 + RAND_SEP;
        rng.rptr = 1;

        let mut discard = 0i32;
        for _ in 0..(RAND_DEG * 10) {
            rng.step(&mut discard);
        }
        rng
    }

    fn step(&mut self, result: &mut i32) {
        let val = (self.state[self.fptr] as u32).wrapping_add(self.state[self.rptr] as u32);
        self.state[self.fptr] = val as i32;
        *result = (val >> 1) as i32;
        self.fptr += 1;
        if self.fptr >= END_PTR {
            self.fptr = 1;
            self.rptr += 1;
        } else {
            self.rptr += 1;
            if self.rptr >= END_PTR {
                self.rptr = 1;
            }
        }
    }

    fn unit(&mut self) -> f64 {
        let mut r = 0i32;
        self.step(&mut r);
        (r as f64) / 2147483647.0
    }
}

pub fn xrandom(state: &mut RngState, xl: f64, xh: f64) -> f64 {
    xl + (xh - xl) * state.unit()
}

pub fn grandom(state: &mut RngState, mean: f64, sdev: f64) -> f64 {
    loop {
        let v1 = xrandom(state, -1.0, 1.0);
        let v2 = xrandom(state, -1.0, 1.0);
        let s = v1 * v1 + v2 * v2;
        if s < 1.0 {
            return mean + sdev * v1 * (-2.0 * s.ln() / s).sqrt();
        }
    }
}

pub fn pickshell(state: &mut RngState, vec: &mut Vector, ndim: usize, rad: f32) {
    loop {
        let mut rsq: f32 = 0.0;
        for v in vec.iter_mut().take(ndim) {
            *v = xrandom(state, -1.0, 1.0) as f32;
            rsq += *v * *v;
        }
        if rsq <= 1.0 {
            let rscale = rad / rsq.sqrt();
            for v in vec.iter_mut().take(ndim) {
                *v *= rscale;
            }
            return;
        }
    }
}

pub fn pickball(state: &mut RngState, vec: &mut Vector, ndim: usize, rad: f32) {
    loop {
        let mut rsq: f32 = 0.0;
        for v in vec.iter_mut().take(ndim) {
            *v = xrandom(state, -1.0, 1.0) as f32;
            rsq += *v * *v;
        }
        if rsq <= 1.0 {
            for v in vec.iter_mut().take(ndim) {
                *v *= rad;
            }
            return;
        }
    }
}

pub fn pickbox(state: &mut RngState, vec: &mut Vector, ndim: usize, size: f32) {
    for v in vec.iter_mut().take(ndim) {
        *v = xrandom(state, -(size as f64), size as f64) as f32;
    }
}
