#![allow(clippy::needless_range_loop)]

use std::os::raw::{c_char, c_double, c_int, c_short, c_void};

pub type Real = f32;
pub type Vector = [Real; 3];
pub type Matrix = [[Real; 3]; 3];

pub const BODY: i16 = 0o1;
pub const CELL: i16 = 0o2;
pub const NDIM: usize = 3;
pub const NSUB: usize = 1 << NDIM;

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Node {
    pub node_type: i16,
    pub update: c_short,
    pub mass: Real,
    pub pos: Vector,
    pub next: *mut Node,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Body {
    pub bodynode: Node,
    pub vel: Vector,
    pub acc: Vector,
    pub phi: Real,
}

#[repr(C)]
pub union Sorq {
    pub subp: [*mut Node; NSUB],
    pub quad: Matrix,
}

impl Copy for Sorq {}
impl Clone for Sorq {
    fn clone(&self) -> Self {
        *self
    }
}

impl std::fmt::Debug for Sorq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Sorq(...)")
    }
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Cell {
    pub cellnode: Node,
    pub rcrit2: Real,
    pub more: *mut Node,
    pub sorq: Sorq,
}

#[no_mangle]
pub static mut root: *mut Cell = std::ptr::null_mut();
#[no_mangle]
pub static mut rsize: Real = 0.0;
#[no_mangle]
pub static mut ncell: c_int = 0;
#[no_mangle]
pub static mut tdepth: c_int = 0;
#[no_mangle]
pub static mut cputree: Real = 0.0;
#[no_mangle]
pub static mut theta: Real = 0.0;
#[no_mangle]
pub static mut options: *mut c_char = std::ptr::null_mut();
#[no_mangle]
pub static mut usequad: u8 = 0;
#[no_mangle]
pub static mut eps: Real = 0.0;
#[no_mangle]
pub static mut actmax: c_int = 0;
#[no_mangle]
pub static mut nbbcalc: c_int = 0;
#[no_mangle]
pub static mut nbccalc: c_int = 0;
#[no_mangle]
pub static mut cpuforce: Real = 0.0;
#[no_mangle]
pub static mut infile: *mut c_char = std::ptr::null_mut();
#[no_mangle]
pub static mut outfile: *mut c_char = std::ptr::null_mut();
#[no_mangle]
pub static mut savefile: *mut c_char = std::ptr::null_mut();
#[no_mangle]
pub static mut dtime: Real = 0.0;
#[no_mangle]
pub static mut dtout: Real = 0.0;
#[no_mangle]
pub static mut tstop: Real = 0.0;
#[no_mangle]
pub static mut headline: *mut c_char = std::ptr::null_mut();
#[no_mangle]
pub static mut tnow: Real = 0.0;
#[no_mangle]
pub static mut tout: Real = 0.0;
#[no_mangle]
pub static mut nstep: c_int = 0;
#[no_mangle]
pub static mut nbody: c_int = 0;
#[no_mangle]
pub static mut bodytab: *mut Body = std::ptr::null_mut();

pub static mut MTOT: Real = 0.0;
pub static mut ETOT: [Real; 3] = [0.0; 3];
pub static mut KETEN: Matrix = [[0.0; NDIM]; NDIM];
pub static mut PETEN: Matrix = [[0.0; NDIM]; NDIM];
pub static mut CMPOS: Vector = [0.0; NDIM];
pub static mut CMVEL: Vector = [0.0; NDIM];
pub static mut AMVEC: Vector = [0.0; NDIM];

extern "C" {
    pub fn maketree(btab: *mut Body, nbody: c_int);
    pub fn gravcalc();
    pub fn inputdata();
    pub fn startoutput();
    pub fn forcereport();
    pub fn output();
    pub fn outputdata();
    pub fn savestate(pattern: *const c_char);
    pub fn restorestate(file: *const c_char);
    pub fn initparam(argv: *mut *mut c_char, defv: *mut *mut c_char);
    pub fn getparam(name: *const c_char) -> *mut c_char;
    pub fn getiparam(name: *const c_char) -> c_int;
    pub fn getdparam(name: *const c_char) -> c_double;
    pub fn getbparam(name: *const c_char) -> c_short;
    pub fn getparamstat(name: *const c_char) -> c_int;

    pub fn allocate(nb: c_int) -> *mut c_void;
    pub fn cputime() -> c_double;
    pub fn error(fmt: *const c_char, ...);
    pub fn eprintf(fmt: *const c_char, ...);
    pub fn scanopt(opt: *const c_char, key: *const c_char) -> c_short;
    pub fn stropen(name: *const c_char, mode: *const c_char) -> *mut std::ffi::c_void;

    pub fn fsqr(x: Real) -> Real;
    pub fn fqbe(x: Real) -> Real;
    pub fn flog2(x: Real) -> Real;
    pub fn fexp2(x: Real) -> Real;
    pub fn fdex(x: Real) -> Real;
    pub fn fcbrt(x: f32) -> f32;
    pub fn xrandom(xl: c_double, xh: c_double) -> c_double;
    pub fn grandom(mean: c_double, sdev: c_double) -> c_double;
    pub fn fpickshell(vec: *mut Real, ndim: c_int, rad: Real);
    pub fn fpickball(vec: *mut Real, ndim: c_int, rad: Real);
    pub fn fpickbox(vec: *mut Real, ndim: c_int, rad: Real);

    pub fn treecode_c_main(argc: c_int, argv: *mut *mut c_char) -> c_int;
}

pub fn vector_zero(v: &mut Vector) {
    for i in 0..NDIM {
        v[i] = 0.0;
    }
}

pub fn vector_length(v: &Vector) -> Real {
    let mut sum = 0.0;
    for i in 0..NDIM {
        sum += v[i] * v[i];
    }
    sum.sqrt()
}

pub fn matrix_zero(m: &mut Matrix) {
    for i in 0..NDIM {
        for j in 0..NDIM {
            m[i][j] = 0.0;
        }
    }
}

pub fn matrix_identity(m: &mut Matrix) {
    for i in 0..NDIM {
        for j in 0..NDIM {
            m[i][j] = if i == j { 1.0 } else { 0.0 };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorq_debug_formats() {
        let q: Sorq = Sorq {
            quad: [[1.0; 3]; 3],
        };
        assert_eq!(format!("{:?}", q), "Sorq(...)");
    }
}
