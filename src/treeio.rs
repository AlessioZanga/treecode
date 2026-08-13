#![allow(clippy::needless_range_loop, clippy::manual_memcpy, static_mut_refs)]

use std::ffi::CStr;
use std::io::{BufReader, BufWriter, Read, Write};
use std::os::raw::c_char;

use crate::error::Result;
use crate::error::TreeError;
use crate::getparam;
use crate::types::{
    actmax, allocate, bodytab, cpuforce, cputime, dtime, dtout, eps, headline, infile, matrix_zero,
    nbbcalc, nbccalc, nbody, ncell, nstep, options, outfile, rsize, savefile, scanopt, tdepth,
    theta, tnow, tout, tstop, usequad, vector_zero, Body, Real, Vector, AMVEC, BODY, CMPOS, CMVEL,
    ETOT, KETEN, MTOT, NDIM, PETEN,
};

const E_WIDTH: usize = 14;

fn c_str_to_string(ptr: *const c_char) -> String {
    if ptr.is_null() {
        return String::new();
    }
    unsafe { CStr::from_ptr(ptr).to_string_lossy().into_owned() }
}

fn getargv0() -> String {
    getparam::getparam("argv0").unwrap()
}

fn getversion() -> String {
    getparam::getparam("VERSION").unwrap()
}

fn alloc_c_string(s: &str) -> Result<*mut c_char> {
    let p = allocate(s.len() + 1)?;
    unsafe {
        std::ptr::copy_nonoverlapping(s.as_ptr(), p, s.len());
        *p.add(s.len()) = 0;
    }
    Ok(p as *mut c_char)
}

fn fmt_e14(v: Real) -> String {
    let s = format!("{:.7E}", v);
    let (mant, exp) = s.split_once('E').unwrap_or((&s, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let sign = if exp < 0 { '-' } else { '+' };
    let body = format!("{}E{}{:02}", mant, sign, exp.abs());
    let mut out = String::new();
    for _ in 0..E_WIDTH.saturating_sub(body.len()) {
        out.push(' ');
    }
    out.push_str(&body);
    out
}

pub fn inputdata() -> Result<()> {
    let toks = read_input_tokens()?;
    let mut i = 0usize;

    let nb = parse_i32(&toks, &mut i)?;
    if nb < 1 {
        return Err(TreeError::AbsurdNbody(nb));
    }
    let ndim = parse_i32(&toks, &mut i)?;
    if ndim != NDIM as i32 {
        return Err(TreeError::BadNdim {
            got: ndim,
            expected: NDIM,
        });
    }
    let t = parse_f64(&toks, &mut i)? as Real;

    unsafe {
        tnow = t;
        nbody = nb;
        bodytab = allocate(nb as usize * std::mem::size_of::<Body>())? as *mut Body;

        for j in 0..nb as usize {
            let p = &mut *bodytab.add(j);
            p.bodynode.mass = parse_f64(&toks, &mut i)? as Real;
        }
        for j in 0..nb as usize {
            let p = &mut *bodytab.add(j);
            for k in 0..NDIM {
                p.bodynode.pos[k] = parse_f64(&toks, &mut i)? as Real;
            }
        }
        for j in 0..nb as usize {
            let p = &mut *bodytab.add(j);
            for k in 0..NDIM {
                p.vel[k] = parse_f64(&toks, &mut i)? as Real;
            }
        }
        for j in 0..nb as usize {
            let p = &mut *bodytab.add(j);
            p.bodynode.node_type = BODY;
        }

        let opts = c_str_to_string(options);
        if scanopt(&opts, "reset-time") {
            tnow = 0.0;
        }
    }
    Ok(())
}

fn read_input_tokens() -> Result<Vec<String>> {
    let in_str = c_str_to_string(unsafe { infile });
    let contents =
        std::fs::read_to_string(&in_str).map_err(|_| TreeError::FileOpen(in_str.clone()))?;
    Ok(contents.split_whitespace().map(|s| s.to_string()).collect())
}

fn parse_i32(toks: &[String], i: &mut usize) -> Result<i32> {
    let s = toks.get(*i).ok_or(TreeError::InputIntConversion)?;
    *i += 1;
    s.parse().map_err(|_| TreeError::InputIntConversion)
}

fn parse_f64(toks: &[String], i: &mut usize) -> Result<f64> {
    let s = toks.get(*i).ok_or(TreeError::InputFloatConversion)?;
    *i += 1;
    s.parse().map_err(|_| TreeError::InputFloatConversion)
}

pub fn startoutput() -> Result<()> {
    let headline_str = c_str_to_string(unsafe { headline });
    println!("\n{}", headline_str);
    let use_quad_str = if unsafe { usequad } != 0 {
        "true"
    } else {
        "false"
    };
    println!(
        "\n{:>8}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}",
        "nbody", "dtime", "eps", "theta", "usequad", "dtout", "tstop"
    );
    unsafe {
        println!(
            "{:8}{:10.5}{:10.4}{:10.2}{:>10}{:10.5}{:10.4}",
            nbody, dtime, eps, theta, use_quad_str, dtout, tstop
        );
    }
    let opts = c_str_to_string(unsafe { options });
    if !opts.is_empty() {
        println!("\n\toptions: {}", opts);
    }
    let save_str = c_str_to_string(unsafe { savefile });
    if !save_str.is_empty() {
        savestate(&save_str)?;
    }
    Ok(())
}

pub fn forcereport() {
    unsafe {
        let ftree = (nbody + ncell - 1) as f32 / ncell as f32;
        println!(
            "\n\t{:>8}{:>8}{:>8}{:>8}{:>10}{:>10}{:>8}",
            "rsize", "tdepth", "ftree", "actmax", "nbbtot", "nbctot", "CPUfc"
        );
        println!(
            "\t{:8.1}{:8}{:8.3}{:8}{:10}{:10}{:8.3}",
            rsize, tdepth, ftree, actmax, nbbcalc, nbccalc, cpuforce
        );
    }
}

pub fn output() -> Result<()> {
    unsafe {
        diagnostics();

        let mut cmabs: Real = 0.0;
        for k in 0..NDIM {
            cmabs += CMVEL[k] * CMVEL[k];
        }
        cmabs = cmabs.sqrt();
        let mut amabs: Real = 0.0;
        for k in 0..NDIM {
            amabs += AMVEC[k] * AMVEC[k];
        }
        amabs = amabs.sqrt();

        println!(
            "\n    {:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            "time", "|T+U|", "T", "-U", "-T/U", "|Vcom|", "|Jtot|", "CPUtot"
        );
        println!(
            "    {:8.3}{:8.5}{:8.5}{:8.5}{:8.5}{:8.5}{:8.5}{:8.3}",
            tnow,
            ETOT[0].abs(),
            ETOT[1],
            -ETOT[2],
            -ETOT[1] / ETOT[2],
            cmabs,
            amabs,
            cputime()?
        );

        let teff = tnow + dtime / 8.0;
        let out_str = c_str_to_string(outfile);
        if !out_str.is_empty() && teff >= tout {
            outputdata()?;
        }
        let save_str = c_str_to_string(savefile);
        if !save_str.is_empty() {
            savestate(&save_str)?;
        }
    }
    Ok(())
}

pub fn outputdata() -> Result<()> {
    let out_str = c_str_to_string(unsafe { outfile });
    let name = out_str.replace("%d", &unsafe { nstep }.to_string());

    let exists = std::path::Path::new(&name).exists();
    let file = if exists {
        std::fs::OpenOptions::new().append(true).open(&name)
    } else {
        std::fs::File::create(&name)
    };
    let mut f = BufWriter::new(file.map_err(|_| TreeError::OutputFileOpen)?);

    unsafe {
        out_int(&mut f, nbody)?;
        out_int(&mut f, NDIM as i32)?;
        out_real(&mut f, tnow)?;
        for j in 0..nbody as usize {
            out_real(&mut f, (*bodytab.add(j)).bodynode.mass)?;
        }
        for j in 0..nbody as usize {
            out_vector(&mut f, (*bodytab.add(j)).bodynode.pos)?;
        }
        for j in 0..nbody as usize {
            out_vector(&mut f, (*bodytab.add(j)).vel)?;
        }
        let opts = c_str_to_string(options);
        if scanopt(&opts, "out-phi") {
            for j in 0..nbody as usize {
                out_real(&mut f, (*bodytab.add(j)).phi)?;
            }
        }
        if scanopt(&opts, "out-acc") {
            for j in 0..nbody as usize {
                out_vector(&mut f, (*bodytab.add(j)).acc)?;
            }
        }
    }

    println!("\n\tdata output to file {} at time {:.6}", name, unsafe {
        tnow
    });
    unsafe {
        tout += dtout;
    }
    Ok(())
}

fn out_int(f: &mut impl Write, v: i32) -> Result<()> {
    let line = format!(" {}\n", v);
    f.write_all(line.as_bytes())
        .map_err(|_| TreeError::WriteFailed)
}

fn out_real(f: &mut impl Write, v: Real) -> Result<()> {
    let line = format!(" {}\n", fmt_e14(v));
    f.write_all(line.as_bytes())
        .map_err(|_| TreeError::WriteFailed)
}

fn out_vector(f: &mut impl Write, v: Vector) -> Result<()> {
    let line = format!(" {} {} {}\n", fmt_e14(v[0]), fmt_e14(v[1]), fmt_e14(v[2]));
    f.write_all(line.as_bytes())
        .map_err(|_| TreeError::WriteFailed)
}

unsafe fn diagnostics() {
    MTOT = 0.0;
    ETOT[1] = 0.0;
    ETOT[2] = 0.0;
    matrix_zero(&mut KETEN);
    matrix_zero(&mut PETEN);
    vector_zero(&mut AMVEC);
    vector_zero(&mut CMPOS);
    vector_zero(&mut CMVEL);

    for j in 0..nbody as usize {
        let p = &*bodytab.add(j);
        let m = p.bodynode.mass;
        MTOT += m;

        let mut velsq: Real = 0.0;
        for k in 0..NDIM {
            velsq += p.vel[k] * p.vel[k];
        }
        ETOT[1] += 0.5 * m * velsq;
        ETOT[2] += 0.5 * m * p.phi;

        for i in 0..NDIM {
            for k in 0..NDIM {
                KETEN[i][k] += (0.5 * m * p.vel[i]) * p.vel[k];
                PETEN[i][k] += (m * p.bodynode.pos[i]) * p.acc[k];
            }
        }

        for i in 0..NDIM {
            let ii = (i + 1) % NDIM;
            let jj = (i + 2) % NDIM;
            AMVEC[i] += m * (p.vel[ii] * p.bodynode.pos[jj] - p.vel[jj] * p.bodynode.pos[ii]);
        }

        for k in 0..NDIM {
            CMPOS[k] += m * p.bodynode.pos[k];
            CMVEL[k] += m * p.vel[k];
        }
    }

    ETOT[0] = ETOT[1] + ETOT[2];
    for k in 0..NDIM {
        CMPOS[k] /= MTOT;
        CMVEL[k] /= MTOT;
    }
}

pub fn savestate(pattern: &str) -> Result<()> {
    let name = if pattern.contains("%d") {
        pattern.replace("%d", &unsafe { nstep & 1 }.to_string())
    } else {
        pattern.to_string()
    };

    let f = std::fs::File::create(&name).map_err(|_| TreeError::FileCreate(name))?;
    let mut f = BufWriter::new(f);

    write_string(&mut f, &getargv0())?;
    write_string(&mut f, &getversion())?;
    unsafe {
        write_real(&mut f, dtime)?;
        write_real(&mut f, theta)?;
        write_bytes(&mut f, &[usequad])?;
        write_real(&mut f, eps)?;
        write_string(&mut f, &c_str_to_string(options))?;
        write_real(&mut f, tstop)?;
        write_real(&mut f, dtout)?;
        write_real(&mut f, tnow)?;
        write_real(&mut f, tout)?;
        write_int(&mut f, nstep)?;
        write_real(&mut f, rsize)?;
        write_int(&mut f, nbody)?;
        write_bodytab(&mut f)?;
    }
    Ok(())
}

fn write_bytes(f: &mut impl Write, buf: &[u8]) -> Result<()> {
    f.write_all(buf).map_err(|_| TreeError::WriteFailed)
}

fn write_int(f: &mut impl Write, v: i32) -> Result<()> {
    write_bytes(f, &v.to_ne_bytes())
}

fn write_real(f: &mut impl Write, v: Real) -> Result<()> {
    write_bytes(f, &v.to_ne_bytes())
}

fn write_string(f: &mut impl Write, s: &str) -> Result<()> {
    write_int(f, (s.len() + 1) as i32)?;
    write_bytes(f, s.as_bytes())?;
    write_bytes(f, &[0u8])
}

fn write_bodytab(f: &mut impl Write) -> Result<()> {
    let nb = unsafe { nbody } as usize;
    let slice = unsafe {
        std::slice::from_raw_parts(bodytab as *const u8, nb * std::mem::size_of::<Body>())
    };
    write_bytes(f, slice)
}

pub fn restorestate(file: &str) -> Result<()> {
    let f = std::fs::File::open(file).map_err(|_| TreeError::FileOpen(file.to_string()))?;
    let mut f = BufReader::new(f);

    let program = read_string(&mut f)?;
    let version = read_string(&mut f)?;
    if program != getargv0() || version != getversion() {
        println!("warning: state file may be outdated\n\n");
    }

    unsafe {
        dtime = read_real(&mut f)?;
        theta = read_real(&mut f)?;
        let mut uq = [0u8; 1];
        read_bytes(&mut f, &mut uq)?;
        usequad = uq[0];
        eps = read_real(&mut f)?;
        let opts = read_string(&mut f)?;
        options = alloc_c_string(&opts)?;
        tstop = read_real(&mut f)?;
        dtout = read_real(&mut f)?;
        tnow = read_real(&mut f)?;
        tout = read_real(&mut f)?;
        nstep = read_int(&mut f)?;
        rsize = read_real(&mut f)?;
        nbody = read_int(&mut f)?;
        let nb = nbody as usize;
        bodytab = allocate(nb * std::mem::size_of::<Body>())? as *mut Body;
        let slice =
            std::slice::from_raw_parts_mut(bodytab as *mut u8, nb * std::mem::size_of::<Body>());
        f.read_exact(slice).map_err(|_| TreeError::ReadFailed)?;
    }
    Ok(())
}

fn read_bytes(f: &mut impl Read, buf: &mut [u8]) -> Result<()> {
    f.read_exact(buf).map_err(|_| TreeError::ReadFailed)
}

fn read_int(f: &mut impl Read) -> Result<i32> {
    let mut b = [0u8; 4];
    read_bytes(f, &mut b)?;
    Ok(i32::from_ne_bytes(b))
}

fn read_real(f: &mut impl Read) -> Result<Real> {
    let mut b = [0u8; 4];
    read_bytes(f, &mut b)?;
    Ok(Real::from_ne_bytes(b))
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
