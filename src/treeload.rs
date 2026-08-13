#![allow(clippy::needless_range_loop, clippy::manual_memcpy, static_mut_refs)]

use crate::mathfns;
use crate::types;

const MAXLEVEL: usize = 32;

static mut FREECELL: *mut types::Node = std::ptr::null_mut();
static mut FIRSTCALL: bool = true;
static mut BH86: bool = false;
static mut SW94: bool = false;
static mut CELLHIST: [i32; MAXLEVEL] = [0; MAXLEVEL];
static mut SUBNHIST: [i32; MAXLEVEL] = [0; MAXLEVEL];

pub fn maketree(btab: &mut [types::Body], nbody: i32) {
    unsafe {
        let cpustart = crate::clib::cputime();
        newtree();
        types::root = makecell();
        types::vector_zero(&mut (*types::root).cellnode.pos);
        expandbox(btab.as_mut_ptr(), nbody);
        load_all_bodies(btab.as_mut_ptr(), nbody);
        parse_options();
        reset_stats();
        hackcofm(types::root, types::rsize, 0);
        threadtree(types::root as *mut types::Node, std::ptr::null_mut());
        if types::usequad != 0 {
            hackquad(types::root);
        }
        types::cputree = (crate::clib::cputime() - cpustart) as types::Real;
    }
}

unsafe fn load_all_bodies(btab: *mut types::Body, nbody: i32) {
    for i in 0..nbody as usize {
        loadbody(btab.add(i));
    }
}

unsafe fn parse_options() {
    let opts = crate::getparam::getparam("options");
    BH86 = crate::clib::scanopt(&opts, "bh86");
    SW94 = crate::clib::scanopt(&opts, "sw94");
    if BH86 && SW94 {
        crate::clib::error("maketree: incompatible options bh86 and sw94\n");
    }
}

unsafe fn reset_stats() {
    types::tdepth = 0;
    for i in 0..MAXLEVEL {
        CELLHIST[i] = 0;
        SUBNHIST[i] = 0;
    }
}

unsafe fn newtree() {
    if !FIRSTCALL {
        let mut p = types::root as *mut types::Node;
        while !p.is_null() {
            if (*p).node_type == types::CELL {
                (*p).next = FREECELL;
                FREECELL = p;
                p = (*(p as *mut types::Cell)).more;
            } else {
                p = (*p).next;
            }
        }
    } else {
        FIRSTCALL = false;
    }
    types::root = std::ptr::null_mut();
    types::ncell = 0;
}

unsafe fn makecell() -> *mut types::Cell {
    let c: *mut types::Cell;
    if FREECELL.is_null() {
        c = crate::clib::allocate(std::mem::size_of::<types::Cell>()) as *mut types::Cell;
    } else {
        c = FREECELL as *mut types::Cell;
        FREECELL = (*FREECELL).next;
    }
    (*c).cellnode.node_type = types::CELL;
    (*c).cellnode.update = 0;
    for i in 0..types::NSUB {
        (*c).sorq.subp[i] = std::ptr::null_mut();
    }
    types::ncell += 1;
    c
}

unsafe fn expandbox(btab: *mut types::Body, nbody: i32) {
    let mut dmax: types::Real = 0.0;
    for i in 0..nbody as usize {
        let p = &*btab.add(i);
        for k in 0..types::NDIM {
            let d = (p.bodynode.pos[k] - (*types::root).cellnode.pos[k]).abs();
            if d > dmax {
                dmax = d;
            }
        }
    }
    while types::rsize < 2.0 * dmax {
        types::rsize *= 2.0;
    }
}

unsafe fn loadbody(p: *mut types::Body) {
    let mut q = types::root;
    let mut qind = subindex(p, q);
    let mut qsize = types::rsize;

    while !(*q).sorq.subp[qind].is_null() {
        if (*(*q).sorq.subp[qind]).node_type == types::BODY {
            let other = (*q).sorq.subp[qind] as *mut types::Body;
            require_distinct_positions(p, other);
            let c = makecell();
            set_cell_midpoint(c, q, p, qsize);
            let sub = subindex(other, c);
            (*c).sorq.subp[sub] = other as *mut types::Node;
            (*q).sorq.subp[qind] = c as *mut types::Node;
        }
        q = (*q).sorq.subp[qind] as *mut types::Cell;
        qind = subindex(p, q);
        qsize /= 2.0;
    }
    (*q).sorq.subp[qind] = p as *mut types::Node;
}

unsafe fn require_distinct_positions(p: *mut types::Body, other: *mut types::Body) {
    let mut dist2: types::Real = 0.0;
    for k in 0..types::NDIM {
        let d = (*p).bodynode.pos[k] - (*other).bodynode.pos[k];
        dist2 += d * d;
    }
    if dist2 == 0.0 {
        crate::clib::error("loadbody: two bodies have same position\n");
    }
}

unsafe fn set_cell_midpoint(
    c: *mut types::Cell,
    q: *mut types::Cell,
    p: *mut types::Body,
    qsize: types::Real,
) {
    for k in 0..types::NDIM {
        let offset = if (*p).bodynode.pos[k] < (*q).cellnode.pos[k] {
            -qsize
        } else {
            qsize
        } / 4.0;
        (*c).cellnode.pos[k] = (*q).cellnode.pos[k] + offset;
    }
}

unsafe fn subindex(p: *mut types::Body, q: *mut types::Cell) -> usize {
    let mut ind: usize = 0;
    for k in 0..types::NDIM {
        if (*q).cellnode.pos[k] <= (*p).bodynode.pos[k] {
            ind += types::NSUB >> (k + 1);
        }
    }
    ind
}

unsafe fn hackcofm(p: *mut types::Cell, psize: types::Real, lev: i32) {
    let mut cmpos: types::Vector = [0.0; types::NDIM];

    update_depth(lev);
    CELLHIST[lev as usize] += 1;
    (*p).cellnode.mass = 0.0;

    accumulate_subnodes(p, psize, lev, &mut cmpos);
    set_center_of_mass(p, &mut cmpos);
    verify_center(p, &cmpos, psize);
    setrcrit(p, &cmpos, psize);
    for k in 0..types::NDIM {
        (*p).cellnode.pos[k] = cmpos[k];
    }
}

unsafe fn update_depth(lev: i32) {
    if lev > types::tdepth {
        types::tdepth = lev;
    }
}

unsafe fn accumulate_subnodes(
    p: *mut types::Cell,
    psize: types::Real,
    lev: i32,
    cmpos: &mut types::Vector,
) {
    let mut tmpv: types::Vector = [0.0; types::NDIM];
    for i in 0..types::NSUB {
        let q = (*p).sorq.subp[i];
        if !q.is_null() {
            SUBNHIST[lev as usize] += 1;
            if (*q).node_type == types::CELL {
                hackcofm(q as *mut types::Cell, psize / 2.0, lev + 1);
            }
            (*p).cellnode.update |= (*q).update;
            (*p).cellnode.mass += (*q).mass;
            for k in 0..types::NDIM {
                tmpv[k] = (*q).pos[k] * (*q).mass;
                cmpos[k] += tmpv[k];
            }
        }
    }
}

unsafe fn set_center_of_mass(p: *mut types::Cell, cmpos: &mut types::Vector) {
    if (*p).cellnode.mass > 0.0 {
        for k in 0..types::NDIM {
            cmpos[k] /= (*p).cellnode.mass;
        }
    } else {
        for k in 0..types::NDIM {
            cmpos[k] = (*p).cellnode.pos[k];
        }
    }
}

unsafe fn verify_center(p: *mut types::Cell, cmpos: &types::Vector, psize: types::Real) {
    for k in 0..types::NDIM {
        if cmpos[k] < (*p).cellnode.pos[k] - psize / 2.0
            || (*p).cellnode.pos[k] + psize / 2.0 <= cmpos[k]
        {
            crate::clib::error("hackcofm: tree structure error\n");
        }
    }
}

unsafe fn setrcrit(p: *mut types::Cell, cmpos: &types::Vector, psize: types::Real) {
    if types::theta == 0.0 {
        set_rcrit_exact(p);
    } else if SW94 {
        set_rcrit_sw94(p, cmpos, psize);
    } else if BH86 {
        set_rcrit_bh86(p, psize);
    } else {
        set_rcrit_default(p, cmpos, psize);
    }
}

unsafe fn set_rcrit_exact(p: *mut types::Cell) {
    (*p).rcrit2 = mathfns::rsqr(2.0 * types::rsize);
}

unsafe fn set_rcrit_sw94(p: *mut types::Cell, cmpos: &types::Vector, psize: types::Real) {
    let mut bmax2: types::Real = 0.0;
    for k in 0..types::NDIM {
        let d = cmpos[k] - (*p).cellnode.pos[k] + psize / 2.0;
        bmax2 += mathfns::rsqr(d.max(psize - d));
    }
    (*p).rcrit2 = bmax2 / mathfns::rsqr(types::theta);
}

unsafe fn set_rcrit_bh86(p: *mut types::Cell, psize: types::Real) {
    (*p).rcrit2 = mathfns::rsqr(psize / types::theta);
}

unsafe fn set_rcrit_default(p: *mut types::Cell, cmpos: &types::Vector, psize: types::Real) {
    let mut d: types::Real = 0.0;
    for k in 0..types::NDIM {
        let dk = cmpos[k] - (*p).cellnode.pos[k];
        d += dk * dk;
    }
    (*p).rcrit2 = mathfns::rsqr(psize / types::theta + d.sqrt());
}

unsafe fn threadtree(p: *mut types::Node, n: *mut types::Node) {
    (*p).next = n;
    if (*p).node_type == types::CELL {
        let c = p as *mut types::Cell;
        let mut desc: [*mut types::Node; types::NSUB + 1] = [std::ptr::null_mut(); types::NSUB + 1];
        let ndesc = collect_descendants(c, &mut desc);
        (*c).more = desc[0];
        desc[ndesc] = n;
        for i in 0..ndesc {
            threadtree(desc[i], desc[i + 1]);
        }
    }
}

unsafe fn hackquad(p: *mut types::Cell) {
    let mut desc: [*mut types::Node; types::NSUB] = [std::ptr::null_mut(); types::NSUB];
    let ndesc = collect_descendants(p, &mut desc);

    types::matrix_zero(&mut (*p).sorq.quad);

    for i in 0..ndesc {
        let q = desc[i];
        if (*q).node_type == types::CELL {
            hackquad(q as *mut types::Cell);
        }
        accumulate_moment(p, q);
    }
}

unsafe fn collect_descendants(c: *mut types::Cell, desc: &mut [*mut types::Node]) -> usize {
    let mut ndesc: usize = 0;
    for i in 0..types::NSUB {
        if !(*c).sorq.subp[i].is_null() {
            desc[ndesc] = (*c).sorq.subp[i];
            ndesc += 1;
        }
    }
    ndesc
}

unsafe fn accumulate_moment(p: *mut types::Cell, q: *mut types::Node) {
    let dr = displacement(q, p);
    let mut tmpm = quadrupole_tensor(q, &dr);
    if (*q).node_type == types::CELL {
        let qm = &*(q as *mut types::Cell);
        for j in 0..types::NDIM {
            for k in 0..types::NDIM {
                tmpm[j][k] += qm.sorq.quad[j][k];
            }
        }
    }
    for j in 0..types::NDIM {
        for k in 0..types::NDIM {
            (*p).sorq.quad[j][k] += tmpm[j][k];
        }
    }
}

unsafe fn displacement(q: *mut types::Node, p: *mut types::Cell) -> types::Vector {
    let mut dr: types::Vector = [0.0; types::NDIM];
    for k in 0..types::NDIM {
        dr[k] = (*q).pos[k] - (*p).cellnode.pos[k];
    }
    dr
}

unsafe fn quadrupole_tensor(q: *mut types::Node, dr: &types::Vector) -> types::Matrix {
    let drsq = dot_product(dr);
    let mut tmpm: types::Matrix = [[0.0; types::NDIM]; types::NDIM];
    for j in 0..types::NDIM {
        for k in 0..types::NDIM {
            let id_rsq = if j == k { drsq } else { 0.0 };
            tmpm[j][k] = (3.0 * dr[j] * dr[k] - id_rsq) * (*q).mass;
        }
    }
    tmpm
}

fn dot_product(dr: &types::Vector) -> types::Real {
    let mut drsq: types::Real = 0.0;
    for k in 0..types::NDIM {
        drsq += dr[k] * dr[k];
    }
    drsq
}

pub fn tree_depth() -> i32 {
    unsafe { types::tdepth }
}

pub fn cell_count() -> i32 {
    unsafe { types::ncell }
}

pub fn tree_build_time() -> f64 {
    unsafe { types::cputree as f64 }
}
