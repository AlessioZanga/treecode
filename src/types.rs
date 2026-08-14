#![allow(clippy::new_without_default)]

use std::os::raw::c_short;

pub fn cputime() -> Result<f64> {
    // Safe, dependency-free replacement for the C `times()` call. The process
    // CPU-time counters (utime + stime, in clock ticks) are read from
    // `/proc/self/stat`; on Linux the tick rate is the fixed USER_HZ == 100,
    // the same value `sysconf(_SC_CLK_TCK)` returns, so the result is identical
    // to the original `libc::times` implementation and the diagnostics stay
    // byte-exact. This keeps the crate 100% safe (`#![forbid(unsafe_code)]`)
    // while preserving C-compatible output.
    let stat = std::fs::read_to_string("/proc/self/stat").map_err(|_| TreeError::CpuTimeFailed)?;
    // Field 2 (`comm`) may contain spaces/parens; split after the final ')'.
    let after_comm = stat
        .rsplit_once(')')
        .map(|(_, rest)| rest)
        .ok_or(TreeError::CpuTimeFailed)?;
    let mut utime: Option<f64> = None;
    let mut stime: Option<f64> = None;
    for (i, tok) in after_comm.split_whitespace().enumerate() {
        match i {
            11 => utime = tok.parse().ok(),
            12 => stime = tok.parse().ok(),
            _ => {}
        }
        if i >= 12 {
            break;
        }
    }
    let utime = utime.ok_or(TreeError::CpuTimeFailed)?;
    let stime = stime.ok_or(TreeError::CpuTimeFailed)?;
    const CLK_TCK: f64 = 100.0; // Linux USER_HZ == sysconf(_SC_CLK_TCK)
    Ok((utime + stime) / (60.0 * CLK_TCK))
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

/// Index of a body within [`Tree::bodytab`].
pub type BodyId = usize;
/// Index of a cell within the [`Tree`] cell arena.
pub type CellId = usize;

/// A reference to either a body or a cell, replacing the C `*mut Node`
/// discriminant-and-pointer pair. Both variants carry a plain index, so the
/// tree can be stored in growable vectors instead of a heap of raw pointers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeRef {
    Body(BodyId),
    Cell(CellId),
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Node {
    pub node_type: i16,
    pub update: c_short,
    pub mass: Real,
    pub pos: Vector,
    pub next: Option<NodeRef>,
}

impl Node {
    pub fn new() -> Self {
        Node {
            node_type: 0,
            update: 0,
            mass: 0.0,
            pos: Vector::zero(),
            next: None,
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

/// Replacement for the C `union Sorq`. Cells either hold the `NSUB` subcell
/// pointers (when quadrupole moments are not in use) or the quadrupole moment
/// matrix. The default is the subcell form so freshly allocated cells behave
/// like the original zeroed union.
#[derive(Debug, Clone, Copy)]
pub enum Sorq {
    Subp([Option<NodeRef>; NSUB]),
    Quad(Matrix),
}

impl Sorq {
    pub fn subp(&self) -> &[Option<NodeRef>; NSUB] {
        match self {
            Sorq::Subp(s) => s,
            Sorq::Quad(_) => panic!("Sorq::subp called on Quad variant"),
        }
    }

    pub fn subp_mut(&mut self) -> &mut [Option<NodeRef>; NSUB] {
        match self {
            Sorq::Subp(s) => s,
            Sorq::Quad(_) => panic!("Sorq::subp_mut called on Quad variant"),
        }
    }

    pub fn quad(&self) -> Matrix {
        match self {
            Sorq::Quad(q) => *q,
            Sorq::Subp(_) => Matrix::zero(),
        }
    }

    pub fn quad_mut(&mut self) -> &mut Matrix {
        match self {
            Sorq::Quad(q) => q,
            Sorq::Subp(_) => {
                *self = Sorq::Quad(Matrix::zero());
                match self {
                    Sorq::Quad(q) => q,
                    _ => unreachable!(),
                }
            }
        }
    }
}

impl Default for Sorq {
    fn default() -> Self {
        Sorq::Subp([None; NSUB])
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Cell {
    pub cellnode: Node,
    pub rcrit2: Real,
    pub more: Option<NodeRef>,
    pub sorq: Sorq,
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            cellnode: Node::new(),
            rcrit2: 0.0,
            more: None,
            sorq: Sorq::default(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Interact {
    pub mass: Real,
    pub pos: Vector,
    pub quad: Option<Matrix>,
}

impl Default for Interact {
    fn default() -> Self {
        Interact {
            mass: 0.0,
            pos: Vector::zero(),
            quad: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorq_debug_formats() {
        let q: Sorq = Sorq::Quad(Matrix([[1.0; 3]; 3]));
        assert!(format!("{:?}", q)
            .contains("Matrix([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])"));
    }
}
