#![allow(non_upper_case_globals)]

use std::os::raw::{c_char, c_int, c_short};

pub fn allocate(nb: usize) -> Result<*mut u8> {
    unsafe {
        let ptr = libc::calloc(nb as libc::size_t, 1) as *mut u8;
        if ptr.is_null() {
            return Err(TreeError::OutOfMemory(nb));
        }
        Ok(ptr)
    }
}

pub fn cputime() -> Result<f64> {
    unsafe {
        let mut buffer: libc::tms = std::mem::zeroed();
        if libc::times(&mut buffer) == -1 {
            return Err(TreeError::CpuTimeFailed);
        }
        let hz = libc::sysconf(libc::_SC_CLK_TCK) as f64;
        Ok((buffer.tms_utime + buffer.tms_stime) as f64 / (60.0 * hz))
    }
}

pub fn scanopt(opt: &str, key: &str) -> bool {
    for word in opt.split(',') {
        if word == key {
            return true;
        }
    }
    false
}

pub use crate::error::eprintf;
pub use crate::error::{Result, TreeError};
pub use crate::vecmath::{matrix_identity, matrix_zero, vector_length, vector_zero};
pub use crate::vecmath::{Matrix, Real, Vector, NDIM};

pub const BODY: i16 = 0o1;
pub const CELL: i16 = 0o2;
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

pub static mut root: *mut Cell = std::ptr::null_mut();
pub static mut rsize: Real = 0.0;
pub static mut ncell: c_int = 0;
pub static mut tdepth: c_int = 0;
pub static mut cputree: Real = 0.0;
pub static mut theta: Real = 0.0;
pub static mut options: *mut c_char = std::ptr::null_mut();
pub static mut usequad: u8 = 0;
pub static mut eps: Real = 0.0;
pub static mut actmax: c_int = 0;
pub static mut nbbcalc: c_int = 0;
pub static mut nbccalc: c_int = 0;
pub static mut cpuforce: Real = 0.0;
pub static mut infile: *mut c_char = std::ptr::null_mut();
pub static mut outfile: *mut c_char = std::ptr::null_mut();
pub static mut savefile: *mut c_char = std::ptr::null_mut();
pub static mut dtime: Real = 0.0;
pub static mut dtout: Real = 0.0;
pub static mut tstop: Real = 0.0;
pub static mut headline: *mut c_char = std::ptr::null_mut();
pub static mut tnow: Real = 0.0;
pub static mut tout: Real = 0.0;
pub static mut nstep: c_int = 0;
pub static mut nbody: c_int = 0;
pub static mut bodytab: *mut Body = std::ptr::null_mut();

pub static mut MTOT: Real = 0.0;
pub static mut ETOT: [Real; 3] = [0.0; 3];
pub static mut KETEN: Matrix = Matrix::zero();
pub static mut PETEN: Matrix = Matrix::zero();
pub static mut CMPOS: Vector = Vector::zero();
pub static mut CMVEL: Vector = Vector::zero();
pub static mut AMVEC: Vector = Vector::zero();

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorq_debug_formats() {
        let q: Sorq = Sorq {
            quad: Matrix([[1.0; 3]; 3]),
        };
        assert_eq!(format!("{:?}", q), "Sorq(...)");
    }
}
