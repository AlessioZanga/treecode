#![allow(
    clippy::needless_range_loop,
    clippy::unnecessary_cast,
    clippy::let_and_return,
    unused_assignments,
    static_mut_refs
)]

use crate::error::Result;
use crate::error::TreeError;
use crate::getparam;
use crate::mathfns;
use crate::treegrav;
use crate::treeio;
use crate::treeload;
use crate::types::{
    allocate, bodytab, dtime, dtout, eps, headline, infile, nbody, nstep, options, outfile, rsize,
    savefile, scanopt, theta, tnow, tout, tstop, usequad, Body, Real, Vector, BODY, NDIM,
};

const MFRAC: f64 = 0.999;

pub fn run(argv: &[&str]) -> Result<()> {
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

    unsafe {
        let headline_str = "Hierarchical N-body code (theta scan)";
        headline = allocate(headline_str.len() + 1)? as *mut std::os::raw::c_char;
        std::ptr::copy_nonoverlapping(
            headline_str.as_ptr(),
            headline as *mut u8,
            headline_str.len(),
        );
        *headline.add(headline_str.len()) = 0;

        startrun()?;
        treeio::startoutput()?;

        if nstep == 0 {
            treeforce()?;
            treeio::output()?;
        }
        while tstop as f64 - tnow as f64 > 0.01 * dtime as f64 {
            stepsystem()?;
            treeio::output()?;
        }
        while tstop as f64 - tnow as f64 > 0.01 * dtime as f64 {
            stepsystem()?;
            treeio::output()?;
        }
    }
    Ok(())
}

unsafe fn treeforce() -> Result<()> {
    let nb = nbody as usize;
    for i in 0..nb {
        (*bodytab.add(i)).bodynode.update = 1;
    }
    treeload::maketree(std::slice::from_raw_parts_mut(bodytab, nb), nbody)?;
    treegrav::gravcalc()?;
    treeio::forcereport();
    Ok(())
}

unsafe fn stepsystem() -> Result<()> {
    let nb = nbody as usize;

    for i in 0..nb {
        let p = &mut *bodytab.add(i);
        for k in 0..NDIM {
            p.vel[k] += p.acc[k] * 0.5 * dtime;
            p.bodynode.pos[k] += p.vel[k] * dtime;
        }
    }

    treeforce()?;

    for i in 0..nb {
        let p = &mut *bodytab.add(i);
        for k in 0..NDIM {
            p.vel[k] += p.acc[k] * 0.5 * dtime;
        }
    }

    nstep += 1;
    tnow += dtime;
    Ok(())
}

unsafe fn startrun() -> Result<()> {
    let in_str = getparam::getparam("in")?;
    let out_str = getparam::getparam("out")?;
    let save_str = getparam::getparam("save")?;

    infile = allocate(in_str.len() + 1)? as *mut std::os::raw::c_char;
    std::ptr::copy_nonoverlapping(in_str.as_ptr(), infile as *mut u8, in_str.len());
    *infile.add(in_str.len()) = 0;

    outfile = allocate(out_str.len() + 1)? as *mut std::os::raw::c_char;
    std::ptr::copy_nonoverlapping(out_str.as_ptr(), outfile as *mut u8, out_str.len());
    *outfile.add(out_str.len()) = 0;

    savefile = allocate(save_str.len() + 1)? as *mut std::os::raw::c_char;
    std::ptr::copy_nonoverlapping(save_str.as_ptr(), savefile as *mut u8, save_str.len());
    *savefile.add(save_str.len()) = 0;

    let restore = getparam::getparam("restore")?;

    if restore.is_empty() {
        eps = getparam::getdparam("eps")? as Real;

        let dtime_str = getparam::getparam("dtime")?;
        dtime = if let Some((n, d)) = dtime_str.split_once('/') {
            let n: f64 = n.parse().unwrap_or(0.0);
            let d: f64 = d.parse().unwrap_or(1.0);
            (n / d) as Real
        } else {
            getparam::getdparam("dtime")? as Real
        };

        theta = getparam::getdparam("theta")? as Real;
        usequad = getparam::getbparam("usequad")? as u8;
        tstop = getparam::getdparam("tstop")? as Real;

        let dtout_str = getparam::getparam("dtout")?;
        dtout = if let Some((n, d)) = dtout_str.split_once('/') {
            let n: f64 = n.parse().unwrap_or(0.0);
            let d: f64 = d.parse().unwrap_or(1.0);
            (n / d) as Real
        } else {
            getparam::getdparam("dtout")? as Real
        };

        let opts = getparam::getparam("options")?;
        let opts_c = std::ffi::CString::new(opts).unwrap();
        options = allocate(opts_c.as_bytes().len() + 1)? as *mut std::os::raw::c_char;
        std::ptr::copy_nonoverlapping(
            opts_c.as_ptr() as *mut u8,
            options as *mut u8,
            opts_c.as_bytes().len() + 1,
        );

        if !in_str.is_empty() {
            treeio::inputdata()?;
        } else {
            nbody = getparam::getiparam("nbody")?;
            if nbody < 1 {
                return Err(TreeError::AbsurdNbody(nbody));
            }
            let seed = getparam::getiparam("seed")?;
            extern "C" {
                fn srandom(seed: u32);
            }
            srandom(seed as u32);
            testdata()?;
            tnow = 0.0;
        }

        rsize = 1.0;
        nstep = 0;
        tout = tnow;
    } else {
        treeio::restorestate(&restore)?;

        if getparam::getparamstat("eps") & 0o4 != 0 {
            eps = getparam::getdparam("eps")? as Real;
        }
        if getparam::getparamstat("theta") & 0o4 != 0 {
            theta = getparam::getdparam("theta")? as Real;
        }
        if getparam::getparamstat("usequad") & 0o4 != 0 {
            usequad = getparam::getbparam("usequad")? as u8;
        }
        if getparam::getparamstat("options") & 0o4 != 0 {
            let opts = getparam::getparam("options")?;
            let opts_c = std::ffi::CString::new(opts).unwrap();
            options = allocate(opts_c.as_bytes().len() + 1)? as *mut std::os::raw::c_char;
            std::ptr::copy_nonoverlapping(
                opts_c.as_ptr() as *mut u8,
                options as *mut u8,
                opts_c.as_bytes().len() + 1,
            );
        }
        if getparam::getparamstat("tstop") & 0o4 != 0 {
            tstop = getparam::getdparam("tstop")? as Real;
        }

        let dtout_str = getparam::getparam("dtout")?;
        dtout = if let Some((n, d)) = dtout_str.split_once('/') {
            let n: f64 = n.parse().unwrap_or(0.0);
            let d: f64 = d.parse().unwrap_or(1.0);
            (n / d) as Real
        } else {
            getparam::getdparam("dtout")? as Real
        };

        let opts = getparam::getparam("options")?;
        if scanopt(&opts, "new-tout") {
            tout = tnow + dtout;
        }
    }
    Ok(())
}

unsafe fn testdata() -> Result<()> {
    let nb = nbody as usize;

    bodytab = allocate(nb * std::mem::size_of::<Body>())? as *mut Body;

    let rsc = 3.0 * std::f32::consts::PI / 16.0;
    let vsc = (1.0 / rsc).sqrt();

    let mut rcm: Vector = Vector::zero();
    let mut vcm: Vector = Vector::zero();

    for i in 0..nb {
        let p = &mut *bodytab.add(i);
        p.bodynode.node_type = BODY;
        p.bodynode.mass = (1.0 / nb as f64) as Real;

        let x_f = mathfns::xrandom(0.0, MFRAC) as f32;
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
        let p = &mut *bodytab.add(i);
        for k in 0..NDIM {
            p.bodynode.pos[k] -= rcm[k];
            p.vel[k] -= vcm[k];
        }
    }
    Ok(())
}
