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
use crate::rng;
use crate::types::{Body, Cell, CellId, Matrix, Real, Vector, BODY, NDIM};

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

    // parsed command-line parameters (was getparam `static mut` table)
    pub config: getparam::Config,

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
            config: getparam::Config::new(),
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
        }
    }

    fn treeforce(&mut self) -> Result<()> {
        let nb = self.nbody as usize;
        for i in 0..nb {
            self.bodytab[i].bodynode.update = 1;
        }
        self.maketree(self.nbody)?;
        self.gravcalc()?;
        self.forcereport();
        Ok(())
    }

    fn stepsystem(&mut self) -> Result<()> {
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

    fn startrun(&mut self) -> Result<()> {
        let in_str = self.config.getparam("in")?;
        let out_str = self.config.getparam("out")?;
        let save_str = self.config.getparam("save")?;

        self.infile = in_str.clone();
        self.outfile = out_str.clone();
        self.savefile = save_str.clone();

        let restore = self.config.getparam("restore")?;

        if restore.is_empty() {
            self.eps = self.config.getdparam("eps")? as Real;

            let dtime_str = self.config.getparam("dtime")?;
            self.dtime = if let Some((n, d)) = dtime_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as Real
            } else {
                self.config.getdparam("dtime")? as Real
            };

            self.theta = self.config.getdparam("theta")? as Real;
            self.usequad = self.config.getbparam("usequad")? as u8;
            self.tstop = self.config.getdparam("tstop")? as Real;

            let dtout_str = self.config.getparam("dtout")?;
            self.dtout = if let Some((n, d)) = dtout_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as Real
            } else {
                self.config.getdparam("dtout")? as Real
            };

            self.options = self.config.getparam("options")?;

            if !in_str.is_empty() {
                self.inputdata()?;
            } else {
                self.nbody = self.config.getiparam("nbody")?;
                if self.nbody < 1 {
                    return Err(TreeError::AbsurdNbody(self.nbody));
                }
                let seed = self.config.getiparam("seed")?;
                let mut rng = rng::RngState::new(seed as u32);
                self.testdata(&mut rng)?;
                self.tnow = 0.0;
            }

            self.rsize = 1.0;
            self.nstep = 0;
            self.tout = self.tnow;
        } else {
            self.restorestate(&restore)?;

            if self.config.getparamstat("eps") & 0o4 != 0 {
                self.eps = self.config.getdparam("eps")? as Real;
            }
            if self.config.getparamstat("theta") & 0o4 != 0 {
                self.theta = self.config.getdparam("theta")? as Real;
            }
            if self.config.getparamstat("usequad") & 0o4 != 0 {
                self.usequad = self.config.getbparam("usequad")? as u8;
            }
            if self.config.getparamstat("options") & 0o4 != 0 {
                self.options = self.config.getparam("options")?;
            }
            if self.config.getparamstat("tstop") & 0o4 != 0 {
                self.tstop = self.config.getdparam("tstop")? as Real;
            }

            let dtout_str = self.config.getparam("dtout")?;
            self.dtout = if let Some((n, d)) = dtout_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as Real
            } else {
                self.config.getdparam("dtout")? as Real
            };

            let opts = self.config.getparam("options")?;
            if crate::types::scanopt(&opts, "new-tout") {
                self.tout = self.tnow + self.dtout;
            }
        }
        Ok(())
    }

    fn testdata(&mut self, rng: &mut rng::RngState) -> Result<()> {
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

            let x_f = rng::xrandom(rng, 0.0, 0.999) as f32;
            let r = (1.0 / (x_f.powf(-2.0 / 3.0) - 1.0).sqrt() as f64) as f32;

            rng::pickshell(rng, &mut p.bodynode.pos, NDIM, rsc * r);

            let mut x: f32 = 0.0;
            let mut y: f32 = 0.0;
            loop {
                x = rng::xrandom(rng, 0.0, 1.0) as f32;
                y = rng::xrandom(rng, 0.0, 0.1) as f32;
                let term = x * x * (1.0 - x * x).powf(3.5);
                if y <= term {
                    break;
                }
            }

            let a = (1.0 + r * r).sqrt();
            let b = (2.0 / a as f64) as f32;
            let v = x * b.sqrt();
            rng::pickshell(rng, &mut p.vel, NDIM, vsc * v);

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

/// Full simulation entry point. Builds a fresh [`Tree`], runs the whole loop,
/// and returns the resulting state so callers can introspect it (e.g.
/// diagnostics, cell counts). Equivalent to constructing a [`Simulation`] and
/// calling [`Simulation::run`].
///
/// [`Simulation::run`]: Simulation::run
pub fn run(argv: &[&str]) -> Result<Tree> {
    let mut sim = Simulation::from_argv(argv)?;
    sim.run()?;
    Ok(sim.into_inner())
}

impl Tree {
    /// Build a fresh [`Tree`] from CLI args and run `startrun` (parameter
    /// parsing, input/test-data loading) without executing the simulation
    /// loop. This is the setup half of [`run`]; see [`Simulation`] for a
    /// granular, ergonomic wrapper.
    ///
    /// [`run`]: crate::treecode::run
    pub fn new_simulation(argv: &[&str]) -> Result<Self> {
        let config = getparam::Config::initparam(argv, &default_defv())?;
        let mut tree = Tree::new();
        tree.config = config;
        tree.headline = "Hierarchical N-body code (theta scan)".to_string();
        tree.startrun()?;
        Ok(tree)
    }
}

fn default_defv() -> [&'static str; 15] {
    [
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
    ]
}

/// Ergonomic, high-level handle to a treecode run, built on the faithful 1:1
/// internals ([`Tree`] and its `maketree`/`gravcalc`/`stepsystem`/… methods).
///
/// [`Simulation::new`] parses parameters and sets up the body data; [`run`]
/// executes the full integration loop (with diagnostics/output). For callers
/// that want finer control, [`step`] advances a single timestep and [`output`]
/// emits diagnostics once.
///
/// [`run`]: Simulation::run
/// [`step`]: Simulation::step
/// [`output`]: Simulation::output
pub struct Simulation(Tree);

impl Simulation {
    /// Parse `params` (as if supplied on the command line) and set up the body
    /// data. No integration is performed yet — call [`run`] to execute the full
    /// loop, or drive it manually with [`step`]/[`output`].
    ///
    /// [`run`]: Simulation::run
    /// [`step`]: Simulation::step
    /// [`output`]: Simulation::output
    ///
    /// # Examples
    /// ```
    /// let mut sim = treecode::Simulation::new([
    ///     "nbody=16",
    ///     "tstop=0.01",
    ///     "dtout=0.005",
    /// ])?;
    /// sim.run()?;
    /// assert!(sim.state().nstep > 0);
    /// # Ok::<(), treecode::error::TreeError>(())
    /// ```
    pub fn new(params: impl IntoIterator<Item = impl Into<String>>) -> Result<Self> {
        let argv: Vec<String> = std::iter::once("treecode".to_string())
            .chain(params.into_iter().map(Into::into))
            .collect();
        let arg_strs: Vec<&str> = argv.iter().map(|s| s.as_str()).collect();
        Simulation::from_argv(&arg_strs)
    }

    fn from_argv(argv: &[&str]) -> Result<Self> {
        Ok(Simulation(Tree::new_simulation(argv)?))
    }

    /// Run the full integration loop: emit the initial diagnostics/output, then
    /// repeatedly advance the system and output until the stop time is reached.
    pub fn run(&mut self) -> Result<()> {
        let tree = &mut self.0;
        tree.startoutput()?;
        if tree.nstep == 0 {
            tree.treeforce()?;
            tree.output()?;
        }
        while (tree.tstop as f64 - tree.tnow as f64) > 0.01 * tree.dtime as f64 {
            tree.stepsystem()?;
            tree.output()?;
        }
        Ok(())
    }

    /// Advance the system by a single timestep (`stepsystem`).
    pub fn step(&mut self) -> Result<()> {
        self.0.stepsystem()
    }

    /// Emit diagnostics and write the current output (`output`).
    pub fn output(&mut self) -> Result<()> {
        self.0.output()
    }

    /// Borrow the underlying [`Tree`] for introspection (diagnostics, counts).
    pub fn state(&self) -> &Tree {
        &self.0
    }

    /// Consume the handle and return the underlying [`Tree`].
    pub fn into_inner(self) -> Tree {
        self.0
    }
}
