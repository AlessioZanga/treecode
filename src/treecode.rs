use crate::{
    error::{Result, TreeError},
    getparam, rng,
    types::{BODY, Body, Cell, CellId, Matrix, NDIM, Vector, scanopt},
};

pub const MAXLEVEL: usize = 32;

/// All formerly-global simulation state, owned by the caller and threaded
/// through the algorithm functions. Field names intentionally preserve the
/// original C global names (including the uppercase diagnostics) so the
/// mapping stays 1:1.
#[derive(Debug)]
pub struct Tree {
    // tree structure (arena-backed; cells live in `cells`, bodies in `bodytab`)
    pub root: Option<CellId>,
    pub bodytab: Vec<Body>,
    pub cells: Vec<Cell>,

    // scalar state
    pub rsize: f32,
    pub ncell: usize,
    pub tdepth: usize,
    pub cputree: f32,
    pub theta: f32,
    pub usequad: u8,
    pub eps: f32,
    pub actmax: usize,
    pub nbbcalc: usize,
    pub nbccalc: usize,
    pub cpuforce: f32,
    pub dtime: f32,
    pub dtout: f32,
    pub tstop: f32,
    pub tnow: f32,
    pub tout: f32,
    pub nstep: usize,
    pub nbody: usize,

    // derived run-constant values, set by `refresh_derived`
    pub inv_nbody: f64,
    pub eps2: f32,
    pub half_dt: f32,
    pub theta_pow_m2_5: f32,
    pub theta2: f32,

    // string state (was `*mut c_char`)
    pub options: String,
    pub infile: String,
    pub outfile: String,
    pub savefile: String,
    pub headline: String,

    // parsed command-line parameters (was getparam `static mut` table)
    pub config: getparam::Config,

    // diagnostics
    pub mtot: f32,
    pub etot: [f32; 3],
    pub keten: Matrix,
    pub peten: Matrix,
    pub cmpos: Vector,
    pub cmvel: Vector,
    pub amvec: Vector,

    // treeload module state (was `static mut`)
    pub freecell: Vec<CellId>,
    pub firstcall: bool,
    pub bh86: bool,
    pub sw94: bool,
    pub cellhist: [usize; MAXLEVEL],
    pub subnhist: [usize; MAXLEVEL],

    // treegrav module state (was `static mut`)
    pub actlen: usize,
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
            inv_nbody: 0.0,
            eps2: 0.0,
            half_dt: 0.0,
            theta_pow_m2_5: 0.0,
            theta2: 0.0,
            options: String::new(),
            infile: String::new(),
            outfile: String::new(),
            savefile: String::new(),
            headline: String::new(),
            config: getparam::Config::new(),
            mtot: 0.0,
            etot: [0.0; 3],
            keten: Matrix::zero(),
            peten: Matrix::zero(),
            cmpos: Vector::zero(),
            cmvel: Vector::zero(),
            amvec: Vector::zero(),
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
        for b in &mut self.bodytab {
            b.bodynode.update = 1;
        }
        self.maketree(self.nbody)?;
        self.gravcalc()?;
        self.forcereport();
        Ok(())
    }

    fn stepsystem(&mut self) -> Result<()> {
        // `0.5 * dtime` is exact (multiplying by a power of two never rounds), so
        // `acc * (0.5 * dtime)` == `(acc * 0.5) * dtime` bit-for-bit. It is
        // precomputed once as `self.half_dt` in `refresh_derived`, turning the
        // per-component kick from two `vmulss` into one, with identical f32
        // results.
        let half_dt = self.half_dt;

        for p in &mut self.bodytab {
            for k in 0..NDIM {
                p.vel[k] += p.acc[k] * half_dt;
                p.bodynode.pos[k] += p.vel[k] * self.dtime;
            }
        }

        self.treeforce()?;

        for p in &mut self.bodytab {
            for k in 0..NDIM {
                p.vel[k] += p.acc[k] * half_dt;
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
            self.eps = self.config.getdparam("eps")? as f32;

            let dtime_str = self.config.getparam("dtime")?;
            self.dtime = if let Some((n, d)) = dtime_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as f32
            } else {
                self.config.getdparam("dtime")? as f32
            };

            self.theta = self.config.getdparam("theta")? as f32;
            self.usequad = self.config.getbparam("usequad")? as u8;
            self.tstop = self.config.getdparam("tstop")? as f32;

            let dtout_str = self.config.getparam("dtout")?;
            self.dtout = if let Some((n, d)) = dtout_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as f32
            } else {
                self.config.getdparam("dtout")? as f32
            };

            self.options = self.config.getparam("options")?;

            if !in_str.is_empty() {
                self.inputdata()?;
                self.refresh_derived();
            } else {
                self.nbody = self.config.getiparam("nbody")? as usize;
                if self.nbody < 1 {
                    return Err(TreeError::AbsurdNbody(self.nbody as i32));
                }
                self.refresh_derived();
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
                self.eps = self.config.getdparam("eps")? as f32;
            }
            if self.config.getparamstat("theta") & 0o4 != 0 {
                self.theta = self.config.getdparam("theta")? as f32;
            }
            if self.config.getparamstat("usequad") & 0o4 != 0 {
                self.usequad = self.config.getbparam("usequad")? as u8;
            }
            if self.config.getparamstat("options") & 0o4 != 0 {
                self.options = self.config.getparam("options")?;
            }
            if self.config.getparamstat("tstop") & 0o4 != 0 {
                self.tstop = self.config.getdparam("tstop")? as f32;
            }

            let dtout_str = self.config.getparam("dtout")?;
            self.dtout = if let Some((n, d)) = dtout_str.split_once('/') {
                let n: f64 = n.parse().unwrap_or(0.0);
                let d: f64 = d.parse().unwrap_or(1.0);
                (n / d) as f32
            } else {
                self.config.getdparam("dtout")? as f32
            };

            let opts = self.config.getparam("options")?;
            if scanopt(&opts, "new-tout") {
                self.tout = self.tnow + self.dtout;
            }
        }
        self.refresh_derived();
        Ok(())
    }

    /// Recompute run-constant derived values once `nbody`/`eps`/`dtime`/`theta`
    /// are known. Each is a pure function of values that do not change during a
    /// run, so computing them here (instead of per body / per force call) is
    /// bit-identical to the original in-loop expressions.
    fn refresh_derived(&mut self) {
        // `1.0 / nbody` must stay a `double` division (C promotes the `1.0`
        // literal to `double` and truncates back to `float` on each store), so
        // `inv_nbody` is kept as `f64`.
        self.inv_nbody = 1.0 / self.nbody as f64;
        self.eps2 = self.eps * self.eps;
        self.half_dt = 0.5 * self.dtime;
        self.theta_pow_m2_5 = self.theta.powf(-2.5);
        self.theta2 = self.theta * self.theta;
    }

    fn testdata(&mut self, rng: &mut rng::RngState) -> Result<()> {
        let nb = self.nbody;

        self.bodytab = (0..nb).map(|_| Body::new()).collect();

        let rsc = 3.0 * std::f32::consts::PI / 16.0;
        let vsc = (1.0 / rsc).sqrt();

        let mut rcm: Vector = Vector::zero();
        let mut vcm: Vector = Vector::zero();

        // `self.inv_nbody` is the loop-invariant `1.0 / nbody` divisor, kept as
        // `f64` (C promotes the `1.0` literal to `double` and truncates back to
        // `float` on each store) — see `refresh_derived`.
        for p in &mut self.bodytab {
            p.bodynode.node_type = BODY;
            p.bodynode.mass = self.inv_nbody as f32;

            let x_f = rng::xrandom(rng, 0.0, 0.999) as f32;
            let r = (1.0 / (x_f.powf(-2.0 / 3.0) - 1.0).sqrt() as f64) as f32;

            rng::pickshell(rng, &mut p.bodynode.pos, NDIM, rsc * r);

            let mut x: f32;
            let mut y: f32;
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
                rcm[k] = (rcm[k] as f64 + p.bodynode.pos[k] as f64 * self.inv_nbody) as f32;
                vcm[k] = (vcm[k] as f64 + p.vel[k] as f64 * self.inv_nbody) as f32;
            }
        }

        for p in &mut self.bodytab {
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
