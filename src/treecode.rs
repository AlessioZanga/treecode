#![allow(
    clippy::needless_range_loop,
    clippy::unnecessary_cast,
    clippy::let_and_return,
    unused_assignments
)]

use std::os::raw::c_int;

use crate::error::Result;
use crate::error::TreeError;
use crate::getparam;
use crate::mathfns;
use crate::types::{Body, Cell, CellId, Matrix, NodeRef, Real, Vector, BODY, NDIM};

pub const MAXLEVEL: usize = 32;

/// All formerly-global simulation state, owned by the caller and threaded
/// through the algorithm functions. Field names intentionally preserve the
/// original C global names (including the uppercase diagnostics) so the
/// mapping stays 1:1.
#[derive(Debug)]
#[allow(non_snake_case)]
pub struct Tree {
    // tree structure (arena-backed; cells live in `cells`, bodies in `bodytab`)
    pub root: Option<CellId>,
    pub bodytab: Vec<Body>,
    pub cells: Vec<Cell>,

    // scalar state
    pub rsize: Real,
    pub ncell: c_int,
    pub tdepth: c_int,
    pub cputree: Real,
    pub theta: Real,
    pub usequad: u8,
    pub eps: Real,
    pub actmax: c_int,
    pub nbbcalc: c_int,
    pub nbccalc: c_int,
    pub cpuforce: Real,
    pub dtime: Real,
    pub dtout: Real,
    pub tstop: Real,
    pub tnow: Real,
    pub tout: Real,
    pub nstep: c_int,
    pub nbody: c_int,

    // string state (was `*mut c_char`)
    pub options: String,
    pub infile: String,
    pub outfile: String,
    pub savefile: String,
    pub headline: String,

    // diagnostics
    pub MTOT: Real,
    pub ETOT: [Real; 3],
    pub KETEN: Matrix,
    pub PETEN: Matrix,
    pub CMPOS: Vector,
    pub CMVEL: Vector,
    pub AMVEC: Vector,

    // treeload module state (was `static mut`)
    pub freecell: Vec<CellId>,
    pub firstcall: bool,
    pub bh86: bool,
    pub sw94: bool,
    pub cellhist: [i32; MAXLEVEL],
    pub subnhist: [i32; MAXLEVEL],

    // treegrav module state (was `static mut`)
    pub actlen: c_int,
    pub active: Vec<NodeRef>,
    pub interact: Vec<Cell>,
}

impl Default for Tree {
    fn default() -> Self {
        Tree::new()
    }
}

impl Tree {
    pub fn new() -> Self {
        Tree {
            root: None,
            bodytab: Vec::new(),
            cells: Vec::new(),
            rsize: 0.0,
            ncell: 0,
            tdepth: 0,
            cputree: 0.0,
            theta: 0.0,
            usequad: 0,
            eps: 0.0,
            actmax: 0,
            nbbcalc: 0,
            nbccalc: 0,
            cpuforce: 0.0,
            dtime: 0.0,
            dtout: 0.0,
            tstop: 0.0,
            tnow: 0.0,
            tout: 0.0,
            nstep: 0,
            nbody: 0,
            options: String::new(),
            infile: String::new(),
            outfile: String::new(),
            savefile: String::new(),
            headline: String::new(),
            MTOT: 0.0,
            ETOT: [0.0; 3],
            KETEN: Matrix::zero(),
            PETEN: Matrix::zero(),
            CMPOS: Vector::zero(),
            CMVEL: Vector::zero(),
            AMVEC: Vector::zero(),
            freecell: Vec::new(),
            firstcall: true,
            bh86: false,
            sw94: false,
            cellhist: [0; MAXLEVEL],
            subnhist: [0; MAXLEVEL],
            actlen: 0,
            active: Vec::new(),
            interact: Vec::new(),
        }
    }

    unsafe fn treeforce(&mut self) -> Result<()> {
        let nb = self.nbody as usize;
        for i in 0..nb {
            self.bodytab[i].bodynode.update = 1;
        }
        self.maketree(self.nbody)?;
        self.gravcalc()?;
        self.forcereport();
        Ok(())
    }

    unsafe fn stepsystem(&mut self) -> Result<()> {
        let nb = self.nbody as usize;

        for i in 0..nb {
            let p = &mut self.bodytab[i];
            for k in 0..NDIM {
                p.vel[k] += p.acc[k] * 0.5 * self.dtime;
                p.bodynode.pos[k] += p.vel[k] * self.dtime;
            }
        }

        self.treeforce()?;

        for i in 0..nb {
            let p = &mut self.bodytab[i];
            for k in 0..NDIM {
                p.vel[k] += p.acc[k] * 0.5 * self.dtime;
            }
        }

        self.nstep += 1;
        self.tnow += self.dtime;
        Ok(())
    }

    unsafe fn startrun(&mut self) -> Result<()> {
        let in_str = getparam::getparam("in")?;
        let out_str = getparam::getparam("out")?;
        let save_str = getparam::getparam("save")?;

        self.infile = in_str.clone();
        self.outfile = out_str.clone();
        self.savefile = save_str.clone();

        let restore = getparam::getparam("restore")?;

        if restore.is_empty() {
            self.eps = getparam::getdparam("eps")? as Real;

            let dtime_str = getparam::getparam("dtime")?;
            self.dtime = if let Some((n, d)) = dtime_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as Real
            } else {
                getparam::getdparam("dtime")? as Real
            };

            self.theta = getparam::getdparam("theta")? as Real;
            self.usequad = getparam::getbparam("usequad")? as u8;
            self.tstop = getparam::getdparam("tstop")? as Real;

            let dtout_str = getparam::getparam("dtout")?;
            self.dtout = if let Some((n, d)) = dtout_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as Real
            } else {
                getparam::getdparam("dtout")? as Real
            };

            self.options = getparam::getparam("options")?;

            if !in_str.is_empty() {
                self.inputdata()?;
            } else {
                self.nbody = getparam::getiparam("nbody")?;
                if self.nbody < 1 {
                    return Err(TreeError::AbsurdNbody(self.nbody));
                }
                let seed = getparam::getiparam("seed")?;
                extern "C" {
                    fn srandom(seed: u32);
                }
                srandom(seed as u32);
                self.testdata()?;
                self.tnow = 0.0;
            }

            self.rsize = 1.0;
            self.nstep = 0;
            self.tout = self.tnow;
        } else {
            self.restorestate(&restore)?;

            if getparam::getparamstat("eps") & 0o4 != 0 {
                self.eps = getparam::getdparam("eps")? as Real;
            }
            if getparam::getparamstat("theta") & 0o4 != 0 {
                self.theta = getparam::getdparam("theta")? as Real;
            }
            if getparam::getparamstat("usequad") & 0o4 != 0 {
                self.usequad = getparam::getbparam("usequad")? as u8;
            }
            if getparam::getparamstat("options") & 0o4 != 0 {
                self.options = getparam::getparam("options")?;
            }
            if getparam::getparamstat("tstop") & 0o4 != 0 {
                self.tstop = getparam::getdparam("tstop")? as Real;
            }

            let dtout_str = getparam::getparam("dtout")?;
            self.dtout = if let Some((n, d)) = dtout_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as Real
            } else {
                getparam::getdparam("dtout")? as Real
            };

            let opts = getparam::getparam("options")?;
            if crate::types::scanopt(&opts, "new-tout") {
                self.tout = self.tnow + self.dtout;
            }
        }
        Ok(())
    }

    unsafe fn testdata(&mut self) -> Result<()> {
        let nb = self.nbody as usize;

        self.bodytab = (0..nb).map(|_| Body::new()).collect();

        let rsc = 3.0 * std::f32::consts::PI / 16.0;
        let vsc = (1.0 / rsc).sqrt();

        let mut rcm: Vector = Vector::zero();
        let mut vcm: Vector = Vector::zero();

        for i in 0..nb {
            let p = &mut self.bodytab[i];
            p.bodynode.node_type = BODY;
            p.bodynode.mass = (1.0 / nb as f64) as Real;

            let x_f = mathfns::xrandom(0.0, 0.999) as f32;
            let r = (1.0 / (x_f.powf(-2.0 / 3.0) - 1.0).sqrt() as f64) as f32;

            mathfns::pickshell(&mut p.bodynode.pos, NDIM, rsc * r);

            let mut x: f32 = 0.0;
            let mut y: f32 = 0.0;
            loop {
                x = mathfns::xrandom(0.0, 1.0) as f32;
                y = mathfns::xrandom(0.0, 0.1) as f32;
                let term = x * x * (1.0 - x * x).powf(3.5);
                if y <= term {
                    break;
                }
            }

            let a = (1.0 + r * r).sqrt();
            let b = (2.0 / a as f64) as f32;
            let v = x * b.sqrt();
            mathfns::pickshell(&mut p.vel, NDIM, vsc * v);

            for k in 0..NDIM {
                rcm[k] = (rcm[k] as f64 + p.bodynode.pos[k] as f64 * (1.0 / nb as f64)) as f32;
                vcm[k] = (vcm[k] as f64 + p.vel[k] as f64 * (1.0 / nb as f64)) as f32;
            }
        }

        for i in 0..nb {
            let p = &mut self.bodytab[i];
            for k in 0..NDIM {
                p.bodynode.pos[k] -= rcm[k];
                p.vel[k] -= vcm[k];
            }
        }
        Ok(())
    }
}

/// Full simulation entry point. Builds a fresh [`Tree`], runs the whole
/// loop, and returns the resulting state so callers can introspect it
/// (e.g. diagnostics, cell counts).
pub fn run(argv: &[&str]) -> Result<Tree> {
    let defv = [
        ";Hierarchical N-body code (theta scan)",
        "in=",
        "out=",
        "dtime=1/32",
        "eps=0.025",
        "theta=1.0",
        "usequad=false",
        "options=",
        "tstop=2.0",
        "dtout=1/4",
        "nbody=4096",
        "seed=123",
        "save=",
        "restore=",
        "VERSION=1.4",
    ];

    getparam::initparam(argv, &defv)?;

    let mut tree = Tree::new();
    tree.headline = "Hierarchical N-body code (theta scan)".to_string();
    unsafe {
        tree.startrun()?;
        tree.startoutput()?;

        if tree.nstep == 0 {
            tree.treeforce()?;
            tree.output()?;
        }
        while (tree.tstop as f64 - tree.tnow as f64) > 0.01 * tree.dtime as f64 {
            tree.stepsystem()?;
            tree.output()?;
        }
        while (tree.tstop as f64 - tree.tnow as f64) > 0.01 * tree.dtime as f64 {
            tree.stepsystem()?;
            tree.output()?;
        }
    }
    Ok(tree)
}
