#![allow(clippy::needless_range_loop, static_mut_refs)]

use crate::types;

const FACTIVE: types::Real = 0.75;

static mut ACTLEN: i32 = 0;
static mut ACTIVE: *mut *mut types::Node = std::ptr::null_mut();
static mut INTERACT: *mut types::Cell = std::ptr::null_mut();

pub fn gravcalc() {
    unsafe {
        let rmid: types::Vector = [0.0; types::NDIM];

        ACTLEN = estimate_active_length();
        ACTIVE = crate::types::allocate(ACTLEN as usize * std::mem::size_of::<*mut types::Node>())
            as *mut *mut types::Node;
        INTERACT = crate::types::allocate(ACTLEN as usize * std::mem::size_of::<types::Cell>())
            as *mut types::Cell;
        let cpustart = crate::types::cputime();
        types::actmax = 0;
        types::nbbcalc = 0;
        types::nbccalc = 0;
        *ACTIVE = types::root as *mut types::Node;
        walktree(
            ACTIVE,
            ACTIVE.add(1),
            INTERACT,
            INTERACT.add(ACTLEN as usize),
            types::root as *mut types::Node,
            types::rsize,
            rmid,
        );
        types::cpuforce = (crate::types::cputime() - cpustart) as types::Real;
        libc::free(ACTIVE as *mut libc::c_void);
        libc::free(INTERACT as *mut libc::c_void);
    }
}

unsafe fn estimate_active_length() -> i32 {
    let base = (FACTIVE * 216.0 * types::tdepth as types::Real) as i32;
    (base as types::Real * types::theta.powf(-2.5)) as i32
}

unsafe fn walktree(
    aptr: *mut *mut types::Node,
    nptr: *mut *mut types::Node,
    cptr: *mut types::Cell,
    bptr: *mut types::Cell,
    p: *mut types::Node,
    psize: types::Real,
    pmid: types::Vector,
) {
    if (*p).update != 0 {
        let mut np = nptr;
        let actsafe = ACTLEN - types::NSUB as i32;
        let mut cptr = cptr;
        let mut bptr = bptr;
        let mut ap = aptr;
        while ap < nptr {
            let apnode = *ap;
            if (*apnode).node_type == types::CELL {
                if accept(apnode, psize, pmid) {
                    let acell = apnode as *mut types::Cell;
                    (*cptr).cellnode.mass = (*apnode).mass;
                    (*cptr).cellnode.pos = (*apnode).pos;
                    (*cptr).sorq.quad = (*acell).sorq.quad;
                    cptr = cptr.add(1);
                } else {
                    if np.offset_from(ACTIVE) >= actsafe as isize {
                        crate::types::error("walktree: active list overflow\n");
                    }
                    let mut q = (*(apnode as *mut types::Cell)).more;
                    while q != (*apnode).next {
                        *np = q;
                        np = np.add(1);
                        q = (*q).next;
                    }
                }
            } else if apnode != p {
                bptr = bptr.sub(1);
                (*bptr).cellnode.mass = (*apnode).mass;
                (*bptr).cellnode.pos = (*apnode).pos;
            }
            ap = ap.add(1);
        }
        let nact = np.offset_from(ACTIVE) as i32;
        if nact > types::actmax {
            types::actmax = nact;
        }
        if np != nptr {
            walksub(nptr, np, cptr, bptr, p, psize, pmid);
        } else {
            if (*p).node_type != types::BODY {
                crate::types::error("walktree: recursion terminated with cell\n");
            }
            gravsum(p as *mut types::Body, cptr, bptr);
        }
    }
}

unsafe fn accept(c: *mut types::Node, psize: types::Real, pmid: types::Vector) -> bool {
    let mut dmax = psize;
    let mut dsq: types::Real = 0.0;
    for k in 0..types::NDIM {
        let mut dk = (*c).pos[k] - pmid[k];
        if dk < 0.0 {
            dk = -dk;
        }
        if dk > dmax {
            dmax = dk;
        }
        dk -= 0.5 * psize;
        if dk > 0.0 {
            dsq += dk * dk;
        }
    }
    let rcrit2 = (*(c as *mut types::Cell)).rcrit2;
    dsq > rcrit2 && dmax > 1.5 * psize
}

unsafe fn walksub(
    nptr: *mut *mut types::Node,
    np: *mut *mut types::Node,
    cptr: *mut types::Cell,
    bptr: *mut types::Cell,
    p: *mut types::Node,
    psize: types::Real,
    pmid: types::Vector,
) {
    let poff = psize / 4.0;
    if (*p).node_type == types::CELL {
        let mut q = (*(p as *mut types::Cell)).more;
        while q != (*p).next {
            let nmid = next_midpoint(pmid, (*q).pos, poff);
            walktree(nptr, np, cptr, bptr, q, psize / 2.0, nmid);
            q = (*q).next;
        }
    } else {
        let nmid = next_midpoint(pmid, (*p).pos, poff);
        walktree(nptr, np, cptr, bptr, p, psize / 2.0, nmid);
    }
}

fn next_midpoint(pmid: types::Vector, pos: types::Vector, poff: types::Real) -> types::Vector {
    let mut nmid = [0.0; types::NDIM];
    for k in 0..types::NDIM {
        nmid[k] = pmid[k] + if pos[k] < pmid[k] { -poff } else { poff };
    }
    nmid
}

unsafe fn gravsum(p0: *mut types::Body, cptr: *mut types::Cell, bptr: *mut types::Cell) {
    let pos0 = (*p0).bodynode.pos;
    let mut phi0: types::Real = 0.0;
    let mut acc0: types::Vector = [0.0; types::NDIM];
    if types::usequad != 0 {
        sumcell(INTERACT, cptr, pos0, &mut phi0, &mut acc0);
    } else {
        sumnode(INTERACT, cptr, pos0, &mut phi0, &mut acc0);
    }
    sumnode(
        bptr,
        INTERACT.add(ACTLEN as usize),
        pos0,
        &mut phi0,
        &mut acc0,
    );
    (*p0).phi = phi0;
    (*p0).acc = acc0;
    types::nbbcalc += INTERACT.add(ACTLEN as usize).offset_from(bptr) as i32;
    types::nbccalc += cptr.offset_from(INTERACT) as i32;
}

unsafe fn sumnode(
    start: *mut types::Cell,
    finish: *mut types::Cell,
    pos0: types::Vector,
    phi0: &mut types::Real,
    acc0: &mut types::Vector,
) {
    let eps2 = types::eps * types::eps;
    let mut p = start;
    while p < finish {
        let (dr, mut dr2) = separation(p, pos0);
        dr2 += eps2;
        let drab = dr2.sqrt();
        let phi_p = (*p).cellnode.mass / drab;
        *phi0 -= phi_p;
        let mr3i = phi_p / dr2;
        add_mul_acc(acc0, &dr, mr3i);
        p = p.add(1);
    }
}

unsafe fn sumcell(
    start: *mut types::Cell,
    finish: *mut types::Cell,
    pos0: types::Vector,
    phi0: &mut types::Real,
    acc0: &mut types::Vector,
) {
    let eps2 = types::eps * types::eps;
    let mut p = start;
    while p < finish {
        let (dr, mut dr2) = separation(p, pos0);
        dr2 += eps2;
        let drab = dr2.sqrt();
        let phi_p = (*p).cellnode.mass / drab;
        let mut mr3i = phi_p / dr2;
        let (qdr, drqdr) = quad_dot(p, &dr);
        let dr5i = 1.0 / (dr2 * dr2 * drab);
        let phi_q = 0.5 * dr5i * drqdr;
        *phi0 -= phi_p + phi_q;
        mr3i += 5.0 * phi_q / dr2;
        add_mul_acc2(acc0, &dr, mr3i, &qdr, -dr5i);
        p = p.add(1);
    }
}

unsafe fn separation(p: *mut types::Cell, pos0: types::Vector) -> (types::Vector, types::Real) {
    let mut dr = [0.0; types::NDIM];
    let mut dr2: types::Real = 0.0;
    for k in 0..types::NDIM {
        dr[k] = (*p).cellnode.pos[k] - pos0[k];
        dr2 += dr[k] * dr[k];
    }
    (dr, dr2)
}

unsafe fn quad_dot(p: *mut types::Cell, dr: &types::Vector) -> (types::Vector, types::Real) {
    let mut qdr = [0.0; types::NDIM];
    let mut drqdr: types::Real = 0.0;
    for i in 0..types::NDIM {
        for j in 0..types::NDIM {
            qdr[i] += (*p).sorq.quad[i][j] * dr[j];
        }
        drqdr += qdr[i] * dr[i];
    }
    (qdr, drqdr)
}

fn add_mul_acc(acc0: &mut types::Vector, dr: &types::Vector, s: types::Real) {
    for k in 0..types::NDIM {
        acc0[k] += dr[k] * s;
    }
}

fn add_mul_acc2(
    acc0: &mut types::Vector,
    dr: &types::Vector,
    s: types::Real,
    w: &types::Vector,
    r: types::Real,
) {
    for k in 0..types::NDIM {
        acc0[k] += dr[k] * s + w[k] * r;
    }
}

pub fn force_max_active() -> i32 {
    unsafe { types::actmax }
}

pub fn force_bb_calc() -> i32 {
    unsafe { types::nbbcalc }
}

pub fn force_bc_calc() -> i32 {
    unsafe { types::nbccalc }
}

pub fn force_cpu_time() -> f64 {
    unsafe { types::cpuforce as f64 }
}
