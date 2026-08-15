use std::{
    fmt::Write as _,
    io::{BufReader, BufWriter, Read, Write},
};

use crate::{
    error::{Result, TreeError},
    getparam,
    treecode::Tree,
    types::{BODY, Body, NDIM, Vector, cputime, scanopt},
    vecmath::{matrix_zero, vector_zero},
};

fn getargv0(config: &getparam::Config) -> Result<String> {
    config.getparam("argv0")
}

fn getversion(config: &getparam::Config) -> Result<String> {
    config.getparam("VERSION")
}

const E_WIDTH: usize = 14;

/// Format a float the way the C reference does with `"%14.7E"`: a 14-wide
/// right-justified field whose exponent is always two digits with an explicit
/// sign (`E+00`/`E-05`). Rust's `{:>14.7E}` omits the sign and zero-pad for a
/// zero exponent (`E0`), so the exponent must be re-emitted explicitly to stay
/// byte-identical to the C output.
fn fmt_e14(v: f32) -> String {
    let s = format!("{:.7E}", v);
    let (mant, exp) = s.split_once('E').unwrap_or((&s, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let sign = if exp < 0 { '-' } else { '+' };
    let body = format!("{}E{}{:02}", mant, sign, exp.abs());
    let mut out = String::with_capacity(E_WIDTH);
    for _ in 0..E_WIDTH.saturating_sub(body.len()) {
        out.push(' ');
    }
    out.push_str(&body);
    out
}

fn parse_count(toks: &[String], i: &mut usize) -> Result<usize> {
    let s = toks.get(*i).ok_or(TreeError::InputIntConversion)?;
    *i += 1;
    s.parse().map_err(|_| TreeError::InputIntConversion)
}

fn parse_f64(toks: &[String], i: &mut usize) -> Result<f64> {
    let s = toks.get(*i).ok_or(TreeError::InputFloatConversion)?;
    *i += 1;
    s.parse().map_err(|_| TreeError::InputFloatConversion)
}

fn out_int(f: &mut impl Write, v: i32) -> Result<()> {
    writeln!(f, " {}", v).map_err(|_| TreeError::WriteFailed)
}

fn out_real(f: &mut impl Write, v: f32) -> Result<()> {
    writeln!(f, " {}", fmt_e14(v)).map_err(|_| TreeError::WriteFailed)
}

fn out_vector(f: &mut impl Write, v: Vector) -> Result<()> {
    writeln!(f, " {} {} {}", fmt_e14(v[0]), fmt_e14(v[1]), fmt_e14(v[2]))
        .map_err(|_| TreeError::WriteFailed)
}

fn write_bytes(f: &mut impl Write, buf: &[u8]) -> Result<()> {
    f.write_all(buf).map_err(|_| TreeError::WriteFailed)
}

fn write_int(f: &mut impl Write, v: i32) -> Result<()> {
    write_bytes(f, &v.to_ne_bytes())
}

fn write_real(f: &mut impl Write, v: f32) -> Result<()> {
    write_bytes(f, &v.to_ne_bytes())
}

fn write_string(f: &mut impl Write, s: &str) -> Result<()> {
    write_int(f, (s.len() + 1) as i32)?;
    write_bytes(f, s.as_bytes())?;
    write_bytes(f, &[0u8])
}

fn read_bytes(f: &mut impl Read, buf: &mut [u8]) -> Result<()> {
    f.read_exact(buf).map_err(|_| TreeError::ReadFailed)
}

fn read_int(f: &mut impl Read) -> Result<i32> {
    let mut b = [0u8; 4];
    read_bytes(f, &mut b)?;
    Ok(i32::from_ne_bytes(b))
}

fn read_real(f: &mut impl Read) -> Result<f32> {
    let mut b = [0u8; 4];
    read_bytes(f, &mut b)?;
    Ok(f32::from_ne_bytes(b))
}

fn read_string(f: &mut impl Read) -> Result<String> {
    let nchars = read_int(f)? as usize;
    let mut buf = vec![0u8; nchars];
    read_bytes(f, &mut buf)?;
    while buf.last() == Some(&0) {
        buf.pop();
    }
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

impl Tree {
    pub fn inputdata(&mut self) -> Result<()> {
        let toks = self.read_input_tokens()?;
        let mut i = 0usize;

        let nb = parse_count(&toks, &mut i)?;
        if nb < 1 {
            return Err(TreeError::AbsurdNbody(nb as i32));
        }
        let ndim = parse_count(&toks, &mut i)?;
        if ndim != NDIM {
            return Err(TreeError::BadNdim {
                got: ndim as i32,
                expected: NDIM,
            });
        }
        let t = parse_f64(&toks, &mut i)? as f32;

        self.tnow = t;
        self.nbody = nb;
        self.bodytab = (0..nb).map(|_| Body::new()).collect();

        for body in &mut self.bodytab {
            body.bodynode.mass = parse_f64(&toks, &mut i)? as f32;
        }
        for body in &mut self.bodytab {
            for p in &mut body.bodynode.pos.0 {
                *p = parse_f64(&toks, &mut i)? as f32;
            }
        }
        for body in &mut self.bodytab {
            for p in &mut body.vel.0 {
                *p = parse_f64(&toks, &mut i)? as f32;
            }
        }
        for body in &mut self.bodytab {
            body.bodynode.node_type = BODY;
        }

        if scanopt(&self.options, "reset-time") {
            self.tnow = 0.0;
        }
        Ok(())
    }

    fn read_input_tokens(&self) -> Result<Vec<String>> {
        let in_str = self.infile.clone();
        let contents =
            std::fs::read_to_string(&in_str).map_err(|_| TreeError::FileOpen(in_str.clone()))?;
        Ok(contents.split_whitespace().map(|s| s.to_string()).collect())
    }

    pub fn startoutput(&self) -> Result<()> {
        let use_quad_str = if self.usequad != 0 { "true" } else { "false" };
        let opts = &self.options;
        let mut s = String::new();
        writeln!(s, "\n{}", self.headline).map_err(|_| TreeError::WriteFailed)?;
        writeln!(
            s,
            "\n{:>8}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}",
            "nbody", "dtime", "eps", "theta", "usequad", "dtout", "tstop"
        )
        .map_err(|_| TreeError::WriteFailed)?;
        writeln!(
            s,
            "{:8}{:10.5}{:10.4}{:10.2}{:>10}{:10.5}{:10.4}",
            self.nbody, self.dtime, self.eps, self.theta, use_quad_str, self.dtout, self.tstop
        )
        .map_err(|_| TreeError::WriteFailed)?;
        if !opts.is_empty() {
            writeln!(s, "\n\toptions: {}", opts).map_err(|_| TreeError::WriteFailed)?;
        }
        print!("{s}");
        let save_str = self.savefile.clone();
        if !save_str.is_empty() {
            self.savestate(&save_str)?;
        }
        Ok(())
    }

    pub fn forcereport(&self) {
        let ftree = (self.nbody + self.ncell - 1) as f32 / self.ncell as f32;
        let mut s = String::new();
        let _ = writeln!(
            s,
            "\n\t{:>8}{:>8}{:>8}{:>8}{:>10}{:>10}{:>8}",
            "rsize", "tdepth", "ftree", "actmax", "nbbtot", "nbctot", "CPUfc"
        );
        let _ = writeln!(
            s,
            "\t{:8.1}{:8}{:8.3}{:8}{:10}{:10}{:8.3}",
            self.rsize, self.tdepth, ftree, self.actmax, self.nbbcalc, self.nbccalc, self.cpuforce
        );
        print!("{s}");
    }

    pub fn output(&mut self) -> Result<()> {
        self.diagnostics();

        let mut cmabs: f32 = 0.0;
        for k in 0..NDIM {
            cmabs += self.cmvel[k] * self.cmvel[k];
        }
        cmabs = cmabs.sqrt();
        let mut amabs: f32 = 0.0;
        for k in 0..NDIM {
            amabs += self.amvec[k] * self.amvec[k];
        }
        amabs = amabs.sqrt();

        let mut s = String::new();
        writeln!(
            s,
            "\n    {:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            "time", "|T+U|", "T", "-U", "-T/U", "|Vcom|", "|Jtot|", "CPUtot"
        )
        .map_err(|_| TreeError::WriteFailed)?;
        writeln!(
            s,
            "    {:8.3}{:8.5}{:8.5}{:8.5}{:8.5}{:8.5}{:8.5}{:8.3}",
            self.tnow,
            self.etot[0].abs(),
            self.etot[1],
            -self.etot[2],
            -self.etot[1] / self.etot[2],
            cmabs,
            amabs,
            cputime()?
        )
        .map_err(|_| TreeError::WriteFailed)?;
        print!("{s}");

        let teff = self.tnow + self.dtime / 8.0;
        let out_str = self.outfile.clone();
        if !out_str.is_empty() && teff >= self.tout {
            self.outputdata()?;
        }
        let save_str = self.savefile.clone();
        if !save_str.is_empty() {
            self.savestate(&save_str)?;
        }
        Ok(())
    }

    pub fn outputdata(&mut self) -> Result<()> {
        let out_str = self.outfile.clone();
        let name = out_str.replace("%d", &self.nstep.to_string());

        let exists = std::path::Path::new(&name).exists();
        let file = if exists {
            std::fs::OpenOptions::new().append(true).open(&name)
        } else {
            std::fs::File::create(&name)
        };
        let mut f = BufWriter::new(file.map_err(|_| TreeError::OutputFileOpen)?);
        self.outputdata_to(&mut f)?;
        println!("\n\tdata output to file {} at time {:.6}", name, self.tnow);
        self.tout += self.dtout;
        Ok(())
    }

    /// Write the body snapshot to an arbitrary [`Write`] sink. This is the
    /// injectable core of [`Tree::outputdata`]; the binary opens the `out=`
    /// file and passes a `BufWriter`, while tests can capture into a `Vec<u8>`.
    pub fn outputdata_to(&self, f: &mut impl Write) -> Result<()> {
        out_int(f, self.nbody as i32)?;
        out_int(f, NDIM as i32)?;
        out_real(f, self.tnow)?;
        for b in &self.bodytab {
            out_real(f, b.bodynode.mass)?;
        }
        for b in &self.bodytab {
            out_vector(f, b.bodynode.pos)?;
        }
        for b in &self.bodytab {
            out_vector(f, b.vel)?;
        }
        let opts = &self.options;
        if scanopt(opts, "out-phi") {
            for b in &self.bodytab {
                out_real(f, b.phi)?;
            }
        }
        if scanopt(opts, "out-acc") {
            for b in &self.bodytab {
                out_vector(f, b.acc)?;
            }
        }
        Ok(())
    }

    fn diagnostics(&mut self) {
        self.mtot = 0.0;
        self.etot[1] = 0.0;
        self.etot[2] = 0.0;
        matrix_zero(&mut self.keten);
        matrix_zero(&mut self.peten);
        vector_zero(&mut self.amvec);
        vector_zero(&mut self.cmpos);
        vector_zero(&mut self.cmvel);

        for p in &self.bodytab {
            let m = p.bodynode.mass;
            self.mtot += m;

            let mut velsq: f32 = 0.0;
            for k in 0..NDIM {
                velsq += p.vel[k] * p.vel[k];
            }
            self.etot[1] += 0.5 * m * velsq;
            self.etot[2] += 0.5 * m * p.phi;

            for i in 0..NDIM {
                for k in 0..NDIM {
                    self.keten[i][k] += (0.5 * m * p.vel[i]) * p.vel[k];
                    self.peten[i][k] += (m * p.bodynode.pos[i]) * p.acc[k];
                }
            }

            for i in 0..NDIM {
                let ii = (i + 1) % NDIM;
                let jj = (i + 2) % NDIM;
                self.amvec[i] +=
                    m * (p.vel[ii] * p.bodynode.pos[jj] - p.vel[jj] * p.bodynode.pos[ii]);
            }

            for k in 0..NDIM {
                self.cmpos[k] += m * p.bodynode.pos[k];
                self.cmvel[k] += m * p.vel[k];
            }
        }

        self.etot[0] = self.etot[1] + self.etot[2];
        for k in 0..NDIM {
            self.cmpos[k] /= self.mtot;
            self.cmvel[k] /= self.mtot;
        }
    }

    pub fn savestate(&self, pattern: &str) -> Result<()> {
        let name = if pattern.contains("%d") {
            pattern.replace("%d", &((self.nstep) & 1).to_string())
        } else {
            pattern.to_string()
        };

        let f = std::fs::File::create(&name).map_err(|_| TreeError::FileCreate(name))?;
        let mut f = BufWriter::new(f);
        self.savestate_to(&mut f)
    }

    /// Serialize the full simulation snapshot to an arbitrary [`Write`] sink.
    /// This is the injectable core of [`Tree::savestate`]; the binary opens the
    /// `save=` file and passes a `BufWriter`, while tests can capture into a
    /// `Vec<u8>`.
    pub fn savestate_to(&self, f: &mut impl Write) -> Result<()> {
        write_string(f, &getargv0(&self.config)?)?;
        write_string(f, &getversion(&self.config)?)?;
        write_real(f, self.dtime)?;
        write_real(f, self.theta)?;
        write_bytes(f, &[self.usequad])?;
        write_real(f, self.eps)?;
        write_string(f, &self.options)?;
        write_real(f, self.tstop)?;
        write_real(f, self.dtout)?;
        write_real(f, self.tnow)?;
        write_real(f, self.tout)?;
        write_int(f, self.nstep as i32)?;
        write_real(f, self.rsize)?;
        write_int(f, self.nbody as i32)?;
        self.write_bodytab(f)?;
        Ok(())
    }

    fn write_bodytab(&self, f: &mut impl Write) -> Result<()> {
        for b in &self.bodytab {
            write_int(f, b.bodynode.node_type as i32)?;
            write_int(f, b.bodynode.update as i32)?;
            write_real(f, b.bodynode.mass)?;
            for &p in &b.bodynode.pos.0 {
                write_real(f, p)?;
            }
            for &p in &b.vel.0 {
                write_real(f, p)?;
            }
            for &p in &b.acc.0 {
                write_real(f, p)?;
            }
            write_real(f, b.phi)?;
        }
        Ok(())
    }

    pub fn restorestate(&mut self, file: &str) -> Result<()> {
        let f = std::fs::File::open(file).map_err(|_| TreeError::FileOpen(file.to_string()))?;
        let mut f = BufReader::new(f);
        self.restorestate_from(&mut f)
    }

    /// Deserialize a simulation snapshot from an arbitrary [`Read`] source. This
    /// is the injectable core of [`Tree::restorestate`]; the binary opens the
    /// `restore=` file and passes a `BufReader`, while tests can replay from a
    /// `Cursor<Vec<u8>>`.
    pub fn restorestate_from(&mut self, f: &mut impl Read) -> Result<()> {
        let program = read_string(f)?;
        let version = read_string(f)?;
        // The saved program/version are compared only as a best-effort warning;
        // a freshly built tree (e.g. `Tree::new()`) may not carry those params.
        if let (Ok(argv0), Ok(ver)) = (getargv0(&self.config), getversion(&self.config)) {
            if program != argv0 || version != ver {
                println!("warning: state file may be outdated\n\n");
            }
        }

        self.dtime = read_real(f)?;
        self.theta = read_real(f)?;
        let mut uq = [0u8; 1];
        read_bytes(f, &mut uq)?;
        self.usequad = uq[0];
        self.eps = read_real(f)?;
        self.options = read_string(f)?;
        self.tstop = read_real(f)?;
        self.dtout = read_real(f)?;
        self.tnow = read_real(f)?;
        self.tout = read_real(f)?;
        self.nstep = read_int(f)? as usize;
        self.rsize = read_real(f)?;
        self.nbody = read_int(f)? as usize;
        self.bodytab = (0..self.nbody).map(|_| Body::new()).collect();
        for b in &mut self.bodytab {
            b.bodynode.node_type = read_int(f)? as i16;
            b.bodynode.update = read_int(f)? as i16;
            b.bodynode.mass = read_real(f)?;
            for p in &mut b.bodynode.pos.0 {
                *p = read_real(f)?;
            }
            for p in &mut b.vel.0 {
                *p = read_real(f)?;
            }
            for p in &mut b.acc.0 {
                *p = read_real(f)?;
            }
            b.phi = read_real(f)?;
        }
        Ok(())
    }
}
