#![allow(clippy::new_without_default)]

use std::os::raw::c_short;

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
pub use crate::error::{error, Result, TreeError};
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

impl Node {
    pub fn new() -> Self {
        Node {
            node_type: 0,
            update: 0,
            mass: 0.0,
            pos: Vector::zero(),
            next: std::ptr::null_mut(),
        }
    }
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Body {
    pub bodynode: Node,
    pub vel: Vector,
    pub acc: Vector,
    pub phi: Real,
}

impl Body {
    pub fn new() -> Self {
        Body {
            bodynode: Node::new(),
            vel: Vector::zero(),
            acc: Vector::zero(),
            phi: 0.0,
        }
    }
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
