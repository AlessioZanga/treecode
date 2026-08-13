#![allow(
    clippy::needless_range_loop,
    clippy::unnecessary_cast,
    clippy::let_and_return,
    unused_assignments,
    static_mut_refs
)]

use crate::mathfns;
use crate::types;

const MFRAC: f64 = 0.999;

pub fn run(argv: &[&str]) {
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

    crate::getparam::initparam(argv, &defv);

    unsafe {
        let headline_str = "Hierarchical N-body code (theta scan)";
        types::headline =
            crate::clib::allocate(headline_str.len() + 1) as *mut std::os::raw::c_char;
        std::ptr::copy_nonoverlapping(
            headline_str.as_ptr(),
            types::headline as *mut u8,
            headline_str.len(),
        );
        *types::headline.add(headline_str.len()) = 0;

        startrun();
        crate::wrapper::startoutput();

        if types::nstep == 0 {
            treeforce();
            crate::wrapper::output();
        }
        while types::tstop as f64 - types::tnow as f64 > 0.01 * types::dtime as f64 {
            stepsystem();
            crate::wrapper::output();
        }
        while types::tstop as f64 - types::tnow as f64 > 0.01 * types::dtime as f64 {
            stepsystem();
            crate::wrapper::output();
        }
    }
}

unsafe fn treeforce() {
    let nbody = types::nbody as usize;
    let bodytab = types::bodytab;
    for i in 0..nbody {
        (*bodytab.add(i)).bodynode.update = 1;
    }
    crate::wrapper::maketree(
        std::slice::from_raw_parts_mut(types::bodytab, nbody),
        types::nbody,
    );
    crate::wrapper::gravcalc();
    crate::wrapper::forcereport();
}

unsafe fn stepsystem() {
    let dtime = types::dtime;
    let nbody = types::nbody as usize;
    let bodytab = types::bodytab;

    for i in 0..nbody {
        let p = &mut *bodytab.add(i);
        for k in 0..types::NDIM {
            p.vel[k] += p.acc[k] * 0.5 * dtime;
            p.bodynode.pos[k] += p.vel[k] * dtime;
        }
    }

    treeforce();

    for i in 0..nbody {
        let p = &mut *bodytab.add(i);
        for k in 0..types::NDIM {
            p.vel[k] += p.acc[k] * 0.5 * dtime;
        }
    }

    types::nstep += 1;
    types::tnow += dtime;
}

unsafe fn startrun() {
    let infile = crate::getparam::getparam("in");
    let outfile = crate::getparam::getparam("out");
    let savefile = crate::getparam::getparam("save");

    types::infile = crate::clib::allocate(infile.len() + 1) as *mut std::os::raw::c_char;
    std::ptr::copy_nonoverlapping(infile.as_ptr(), types::infile as *mut u8, infile.len());
    *types::infile.add(infile.len()) = 0;

    types::outfile = crate::clib::allocate(outfile.len() + 1) as *mut std::os::raw::c_char;
    std::ptr::copy_nonoverlapping(outfile.as_ptr(), types::outfile as *mut u8, outfile.len());
    *types::outfile.add(outfile.len()) = 0;

    types::savefile = crate::clib::allocate(savefile.len() + 1) as *mut std::os::raw::c_char;
    std::ptr::copy_nonoverlapping(
        savefile.as_ptr(),
        types::savefile as *mut u8,
        savefile.len(),
    );
    *types::savefile.add(savefile.len()) = 0;

    let restore = crate::getparam::getparam("restore");

    if restore.is_empty() {
        types::eps = crate::getparam::getdparam("eps") as types::Real;

        let dtime_str = crate::getparam::getparam("dtime");
        types::dtime = if let Some((n, d)) = dtime_str.split_once('/') {
            let n: f64 = n.parse().unwrap_or(0.0);
            let d: f64 = d.parse().unwrap_or(1.0);
            (n / d) as types::Real
        } else {
            crate::getparam::getdparam("dtime") as types::Real
        };

        types::theta = crate::getparam::getdparam("theta") as types::Real;
        types::usequad = crate::getparam::getbparam("usequad") as u8;
        types::tstop = crate::getparam::getdparam("tstop") as types::Real;

        let dtout_str = crate::getparam::getparam("dtout");
        types::dtout = if let Some((n, d)) = dtout_str.split_once('/') {
            let n: f64 = n.parse().unwrap_or(0.0);
            let d: f64 = d.parse().unwrap_or(1.0);
            (n / d) as types::Real
        } else {
            crate::getparam::getdparam("dtout") as types::Real
        };

        let opts = crate::getparam::getparam("options");
        let opts_c = std::ffi::CString::new(opts).unwrap();
        types::options =
            crate::clib::allocate(opts_c.as_bytes().len() + 1) as *mut std::os::raw::c_char;
        std::ptr::copy_nonoverlapping(
            opts_c.as_ptr() as *mut u8,
            types::options as *mut u8,
            opts_c.as_bytes().len() + 1,
        );

        if !infile.is_empty() {
            crate::wrapper::inputdata();
        } else {
            types::nbody = crate::getparam::getiparam("nbody");
            if types::nbody < 1 {
                crate::clib::error("startrun: absurd value for nbody\n");
            }
            let seed = crate::getparam::getiparam("seed");
            extern "C" {
                fn srandom(seed: u32);
            }
            srandom(seed as u32);
            testdata();
            types::tnow = 0.0;
        }

        types::rsize = 1.0;
        types::nstep = 0;
        types::tout = types::tnow;
    } else {
        crate::wrapper::restorestate(&restore);

        if crate::getparam::getparamstat("eps") & 0o4 != 0 {
            types::eps = crate::getparam::getdparam("eps") as types::Real;
        }
        if crate::getparam::getparamstat("theta") & 0o4 != 0 {
            types::theta = crate::getparam::getdparam("theta") as types::Real;
        }
        if crate::getparam::getparamstat("usequad") & 0o4 != 0 {
            types::usequad = crate::getparam::getbparam("usequad") as u8;
        }
        if crate::getparam::getparamstat("options") & 0o4 != 0 {
            let opts = crate::getparam::getparam("options");
            let opts_c = std::ffi::CString::new(opts).unwrap();
            types::options =
                crate::clib::allocate(opts_c.as_bytes().len() + 1) as *mut std::os::raw::c_char;
            std::ptr::copy_nonoverlapping(
                opts_c.as_ptr() as *mut u8,
                types::options as *mut u8,
                opts_c.as_bytes().len() + 1,
            );
        }
        if crate::getparam::getparamstat("tstop") & 0o4 != 0 {
            types::tstop = crate::getparam::getdparam("tstop") as types::Real;
        }

        let dtout_str = crate::getparam::getparam("dtout");
        types::dtout = if let Some((n, d)) = dtout_str.split_once('/') {
            let n: f64 = n.parse().unwrap_or(0.0);
            let d: f64 = d.parse().unwrap_or(1.0);
            (n / d) as types::Real
        } else {
            crate::getparam::getdparam("dtout") as types::Real
        };

        let opts = crate::getparam::getparam("options");
        if crate::clib::scanopt(&opts, "new-tout") {
            types::tout = types::tnow + types::dtout;
        }
    }
}

unsafe fn testdata() {
    let nbody = types::nbody as usize;
    let nb = types::nbody;

    types::bodytab =
        crate::clib::allocate(nbody * std::mem::size_of::<types::Body>()) as *mut types::Body;

    let rsc = ((3.0 * std::f64::consts::PI) / 16.0) as f32;
    let vsc = ((1.0 / rsc as f64) as f32).sqrt();

    let mut rcm: types::Vector = [0.0; types::NDIM];
    let mut vcm: types::Vector = [0.0; types::NDIM];

    for i in 0..nbody {
        let p = &mut *types::bodytab.add(i);
        p.bodynode.node_type = types::BODY;
        p.bodynode.mass = (1.0 / nb as f64) as types::Real;

        let x_f = mathfns::xrandom(0.0, MFRAC) as f32;
        let r = (1.0 / (x_f.powf(-2.0 / 3.0) - 1.0).sqrt() as f64) as f32;

        mathfns::pickshell(&mut p.bodynode.pos, types::NDIM, rsc * r);

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
        mathfns::pickshell(&mut p.vel, types::NDIM, vsc * v);

        for k in 0..types::NDIM {
            rcm[k] = (rcm[k] as f64 + p.bodynode.pos[k] as f64 * (1.0 / nb as f64)) as f32;
            vcm[k] = (vcm[k] as f64 + p.vel[k] as f64 * (1.0 / nb as f64)) as f32;
        }
    }

    for i in 0..nbody {
        let p = &mut *types::bodytab.add(i);
        for k in 0..types::NDIM {
            p.bodynode.pos[k] -= rcm[k];
            p.vel[k] -= vcm[k];
        }
    }
}
