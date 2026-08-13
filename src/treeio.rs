#![allow(clippy::needless_range_loop, clippy::manual_memcpy, static_mut_refs)]

use std::ffi::CStr;
use std::io::{BufReader, BufWriter, Read, Write};
use std::os::raw::c_char;

use crate::types;

const E_WIDTH: usize = 14;

fn c_str_to_string(ptr: *const c_char) -> String {
    if ptr.is_null() {
        return String::new();
    }
    unsafe { CStr::from_ptr(ptr).to_string_lossy().into_owned() }
}

fn getargv0() -> String {
    crate::getparam::getparam("argv0")
}

fn getversion() -> String {
    crate::getparam::getparam("VERSION")
}

fn alloc_c_string(s: &str) -> *mut c_char {
    let p = crate::types::allocate(s.len() + 1);
    unsafe {
        std::ptr::copy_nonoverlapping(s.as_ptr(), p, s.len());
        *p.add(s.len()) = 0;
    }
    p as *mut c_char
}

fn fmt_e14(v: types::Real) -> String {
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

pub fn inputdata() {
    let toks = read_input_tokens();
    let mut i = 0usize;

    let nbody = parse_i32(&toks, &mut i);
    if nbody < 1 {
        crate::types::error(&format!("inputdata: nbody = {} is absurd\n", nbody));
    }
    let ndim = parse_i32(&toks, &mut i);
    if ndim != types::NDIM as i32 {
        crate::types::error(&format!(
            "inputdata: ndim = {}; expected {}\n",
            ndim,
            types::NDIM
        ));
    }
    let tnow = parse_f64(&toks, &mut i) as types::Real;

    unsafe {
        types::tnow = tnow;
        types::nbody = nbody;
        types::bodytab = crate::types::allocate(nbody as usize * std::mem::size_of::<types::Body>())
            as *mut types::Body;

        for j in 0..nbody as usize {
            let p = &mut *types::bodytab.add(j);
            p.bodynode.mass = parse_f64(&toks, &mut i) as types::Real;
        }
        for j in 0..nbody as usize {
            let p = &mut *types::bodytab.add(j);
            for k in 0..types::NDIM {
                p.bodynode.pos[k] = parse_f64(&toks, &mut i) as types::Real;
            }
        }
        for j in 0..nbody as usize {
            let p = &mut *types::bodytab.add(j);
            for k in 0..types::NDIM {
                p.vel[k] = parse_f64(&toks, &mut i) as types::Real;
            }
        }
        for j in 0..nbody as usize {
            let p = &mut *types::bodytab.add(j);
            p.bodynode.node_type = types::BODY;
        }

        let opts = c_str_to_string(types::options);
        if crate::types::scanopt(&opts, "reset-time") {
            types::tnow = 0.0;
        }
    }
}

fn read_input_tokens() -> Vec<String> {
    let infile = c_str_to_string(unsafe { types::infile });
    let contents = std::fs::read_to_string(&infile).unwrap_or_else(|_| {
        crate::types::error(&format!("inputdata: cannot open file \"{}\"\n", infile));
        unreachable!()
    });
    contents.split_whitespace().map(|s| s.to_string()).collect()
}

fn parse_i32(toks: &[String], i: &mut usize) -> i32 {
    let s = toks.get(*i).unwrap_or_else(|| {
        crate::types::error("in_int: input conversion error\n");
        unreachable!()
    });
    *i += 1;
    s.parse().unwrap_or_else(|_| {
        crate::types::error("in_int: input conversion error\n");
        unreachable!()
    })
}

fn parse_f64(toks: &[String], i: &mut usize) -> f64 {
    let s = toks.get(*i).unwrap_or_else(|| {
        crate::types::error("in_real: input conversion error\n");
        unreachable!()
    });
    *i += 1;
    s.parse().unwrap_or_else(|_| {
        crate::types::error("in_real: input conversion error\n");
        unreachable!()
    })
}

pub fn startoutput() {
    let headline = c_str_to_string(unsafe { types::headline });
    println!("\n{}", headline);
    let usequad = if unsafe { types::usequad } != 0 {
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
            types::nbody,
            types::dtime,
            types::eps,
            types::theta,
            usequad,
            types::dtout,
            types::tstop
        );
    }
    let opts = c_str_to_string(unsafe { types::options });
    if !opts.is_empty() {
        println!("\n\toptions: {}", opts);
    }
    let savefile = c_str_to_string(unsafe { types::savefile });
    if !savefile.is_empty() {
        savestate(&savefile);
    }
}

pub fn forcereport() {
    unsafe {
        let ftree = (types::nbody + types::ncell - 1) as f32 / types::ncell as f32;
        println!(
            "\n\t{:>8}{:>8}{:>8}{:>8}{:>10}{:>10}{:>8}",
            "rsize", "tdepth", "ftree", "actmax", "nbbtot", "nbctot", "CPUfc"
        );
        println!(
            "\t{:8.1}{:8}{:8.3}{:8}{:10}{:10}{:8.3}",
            types::rsize,
            types::tdepth,
            ftree,
            types::actmax,
            types::nbbcalc,
            types::nbccalc,
            types::cpuforce
        );
    }
}

pub fn output() {
    unsafe {
        diagnostics();

        let mut cmabs: types::Real = 0.0;
        for k in 0..types::NDIM {
            cmabs += types::CMVEL[k] * types::CMVEL[k];
        }
        cmabs = cmabs.sqrt();
        let mut amabs: types::Real = 0.0;
        for k in 0..types::NDIM {
            amabs += types::AMVEC[k] * types::AMVEC[k];
        }
        amabs = amabs.sqrt();

        println!(
            "\n    {:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            "time", "|T+U|", "T", "-U", "-T/U", "|Vcom|", "|Jtot|", "CPUtot"
        );
        println!(
            "    {:8.3}{:8.5}{:8.5}{:8.5}{:8.5}{:8.5}{:8.5}{:8.3}",
            types::tnow,
            types::ETOT[0].abs(),
            types::ETOT[1],
            -types::ETOT[2],
            -types::ETOT[1] / types::ETOT[2],
            cmabs,
            amabs,
            crate::types::cputime()
        );

        let teff = types::tnow + types::dtime / 8.0;
        let outfile = c_str_to_string(types::outfile);
        if !outfile.is_empty() && teff >= types::tout {
            outputdata();
        }
        let savefile = c_str_to_string(types::savefile);
        if !savefile.is_empty() {
            savestate(&savefile);
        }
    }
}

pub fn outputdata() {
    let outfile = c_str_to_string(unsafe { types::outfile });
    let name = outfile.replace("%d", &unsafe { types::nstep }.to_string());

    let exists = std::path::Path::new(&name).exists();
    let file = if exists {
        std::fs::OpenOptions::new().append(true).open(&name)
    } else {
        std::fs::File::create(&name)
    };
    let mut f = BufWriter::new(match file {
        Ok(f) => f,
        Err(_) => {
            crate::types::error("outputdata: cannot open output file\n");
            unreachable!()
        }
    });

    unsafe {
        out_int(&mut f, types::nbody);
        out_int(&mut f, types::NDIM as i32);
        out_real(&mut f, types::tnow);
        for j in 0..types::nbody as usize {
            out_real(&mut f, (*types::bodytab.add(j)).bodynode.mass);
        }
        for j in 0..types::nbody as usize {
            out_vector(&mut f, (*types::bodytab.add(j)).bodynode.pos);
        }
        for j in 0..types::nbody as usize {
            out_vector(&mut f, (*types::bodytab.add(j)).vel);
        }
        let opts = c_str_to_string(types::options);
        if crate::types::scanopt(&opts, "out-phi") {
            for j in 0..types::nbody as usize {
                out_real(&mut f, (*types::bodytab.add(j)).phi);
            }
        }
        if crate::types::scanopt(&opts, "out-acc") {
            for j in 0..types::nbody as usize {
                out_vector(&mut f, (*types::bodytab.add(j)).acc);
            }
        }
    }

    println!("\n\tdata output to file {} at time {:.6}", name, unsafe {
        types::tnow
    });
    unsafe {
        types::tout += types::dtout;
    }
}

fn out_int(f: &mut impl Write, v: i32) {
    let line = format!(" {}\n", v);
    if f.write_all(line.as_bytes()).is_err() {
        crate::types::error("out_int: fprintf failed\n");
    }
}

fn out_real(f: &mut impl Write, v: types::Real) {
    let line = format!(" {}\n", fmt_e14(v));
    if f.write_all(line.as_bytes()).is_err() {
        crate::types::error("out_real: fprintf failed\n");
    }
}

fn out_vector(f: &mut impl Write, v: types::Vector) {
    let line = format!(" {} {} {}\n", fmt_e14(v[0]), fmt_e14(v[1]), fmt_e14(v[2]));
    if f.write_all(line.as_bytes()).is_err() {
        crate::types::error("out_vector: fprintf failed\n");
    }
}

unsafe fn diagnostics() {
    types::MTOT = 0.0;
    types::ETOT[1] = 0.0;
    types::ETOT[2] = 0.0;
    types::matrix_zero(&mut types::KETEN);
    types::matrix_zero(&mut types::PETEN);
    types::vector_zero(&mut types::AMVEC);
    types::vector_zero(&mut types::CMPOS);
    types::vector_zero(&mut types::CMVEL);

    for j in 0..types::nbody as usize {
        let p = &*types::bodytab.add(j);
        let m = p.bodynode.mass;
        types::MTOT += m;

        let mut velsq: types::Real = 0.0;
        for k in 0..types::NDIM {
            velsq += p.vel[k] * p.vel[k];
        }
        types::ETOT[1] += 0.5 * m * velsq;
        types::ETOT[2] += 0.5 * m * p.phi;

        for i in 0..types::NDIM {
            for k in 0..types::NDIM {
                types::KETEN[i][k] += (0.5 * m * p.vel[i]) * p.vel[k];
                types::PETEN[i][k] += (m * p.bodynode.pos[i]) * p.acc[k];
            }
        }

        for i in 0..types::NDIM {
            let ii = (i + 1) % types::NDIM;
            let jj = (i + 2) % types::NDIM;
            types::AMVEC[i] +=
                m * (p.vel[ii] * p.bodynode.pos[jj] - p.vel[jj] * p.bodynode.pos[ii]);
        }

        for k in 0..types::NDIM {
            types::CMPOS[k] += m * p.bodynode.pos[k];
            types::CMVEL[k] += m * p.vel[k];
        }
    }

    types::ETOT[0] = types::ETOT[1] + types::ETOT[2];
    for k in 0..types::NDIM {
        types::CMPOS[k] /= types::MTOT;
        types::CMVEL[k] /= types::MTOT;
    }
}

pub fn savestate(pattern: &str) {
    let name = if pattern.contains("%d") {
        pattern.replace("%d", &unsafe { types::nstep & 1 }.to_string())
    } else {
        pattern.to_string()
    };

    let file = std::fs::File::create(&name);
    let mut f = BufWriter::new(match file {
        Ok(f) => f,
        Err(_) => {
            crate::types::error("savestate: cannot create file\n");
            unreachable!()
        }
    });

    write_string(&mut f, &getargv0());
    write_string(&mut f, &getversion());
    unsafe {
        write_real(&mut f, types::dtime);
        write_real(&mut f, types::theta);
        write_bytes(&mut f, &[types::usequad]);
        write_real(&mut f, types::eps);
        write_string(&mut f, &c_str_to_string(types::options));
        write_real(&mut f, types::tstop);
        write_real(&mut f, types::dtout);
        write_real(&mut f, types::tnow);
        write_real(&mut f, types::tout);
        write_int(&mut f, types::nstep);
        write_real(&mut f, types::rsize);
        write_int(&mut f, types::nbody);
        write_bodytab(&mut f);
    }
}

fn write_bytes(f: &mut impl Write, buf: &[u8]) {
    if f.write_all(buf).is_err() {
        crate::types::error("savestate: fwrite failed\n");
    }
}

fn write_int(f: &mut impl Write, v: i32) {
    write_bytes(f, &v.to_ne_bytes());
}

fn write_real(f: &mut impl Write, v: types::Real) {
    write_bytes(f, &v.to_ne_bytes());
}

fn write_string(f: &mut impl Write, s: &str) {
    write_int(f, (s.len() + 1) as i32);
    write_bytes(f, s.as_bytes());
    write_bytes(f, &[0u8]);
}

fn write_bodytab(f: &mut impl Write) {
    let nb = unsafe { types::nbody } as usize;
    let slice = unsafe {
        std::slice::from_raw_parts(
            types::bodytab as *const u8,
            nb * std::mem::size_of::<types::Body>(),
        )
    };
    write_bytes(f, slice);
}

pub fn restorestate(file: &str) {
    let f = std::fs::File::open(file);
    let mut f = BufReader::new(match f {
        Ok(f) => f,
        Err(_) => {
            crate::types::error("restorestate: cannot open file\n");
            unreachable!()
        }
    });

    let program = read_string(&mut f);
    let version = read_string(&mut f);
    if program != getargv0() || version != getversion() {
        println!("warning: state file may be outdated\n\n");
    }

    unsafe {
        types::dtime = read_real(&mut f);
        types::theta = read_real(&mut f);
        let mut uq = [0u8; 1];
        read_bytes(&mut f, &mut uq);
        types::usequad = uq[0];
        types::eps = read_real(&mut f);
        let opts = read_string(&mut f);
        types::options = alloc_c_string(&opts);
        types::tstop = read_real(&mut f);
        types::dtout = read_real(&mut f);
        types::tnow = read_real(&mut f);
        types::tout = read_real(&mut f);
        types::nstep = read_int(&mut f);
        types::rsize = read_real(&mut f);
        types::nbody = read_int(&mut f);
        let nb = types::nbody as usize;
        types::bodytab =
            crate::types::allocate(nb * std::mem::size_of::<types::Body>()) as *mut types::Body;
        let slice = std::slice::from_raw_parts_mut(
            types::bodytab as *mut u8,
            nb * std::mem::size_of::<types::Body>(),
        );
        if f.read_exact(slice).is_err() {
            crate::types::error("restorestate: fread failed\n");
        }
    }
}

fn read_bytes(f: &mut impl Read, buf: &mut [u8]) {
    if f.read_exact(buf).is_err() {
        crate::types::error("restorestate: fread failed\n");
    }
}

fn read_int(f: &mut impl Read) -> i32 {
    let mut b = [0u8; 4];
    read_bytes(f, &mut b);
    i32::from_ne_bytes(b)
}

fn read_real(f: &mut impl Read) -> types::Real {
    let mut b = [0u8; 4];
    read_bytes(f, &mut b);
    types::Real::from_ne_bytes(b)
}

fn read_string(f: &mut impl Read) -> String {
    let nchars = read_int(f) as usize;
    let mut buf = vec![0u8; nchars];
    read_bytes(f, &mut buf);
    while buf.last() == Some(&0) {
        buf.pop();
    }
    String::from_utf8_lossy(&buf).into_owned()
}
