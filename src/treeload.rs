#![allow(
    clippy::needless_range_loop,
    clippy::manual_memcpy,
    static_mut_refs,
    // Raw-pointer tree access is temporary (removed by the Phase 2 arena);
    // indexing a `Deref` field through a raw-pointer deref is intentional here.
    dangerous_implicit_autorefs
)]

use crate::error::Result;
use crate::error::TreeError;
use crate::getparam;
use crate::mathfns;
use crate::types::{
    allocate, cputime, cputree, matrix_zero, ncell, root, rsize, scanopt, tdepth, theta, usequad,
    vector_zero, Body, Cell, Matrix, Node, Real, Vector, BODY, CELL, NDIM, NSUB,
};

const MAXLEVEL: usize = 32;

static mut FREECELL: *mut Node = std::ptr::null_mut();
static mut FIRSTCALL: bool = true;
static mut BH86: bool = false;
static mut SW94: bool = false;
static mut CELLHIST: [i32; MAXLEVEL] = [0; MAXLEVEL];
static mut SUBNHIST: [i32; MAXLEVEL] = [0; MAXLEVEL];

pub fn maketree(btab: &mut [Body], nbody: i32) -> Result<()> {
    unsafe {
        let cpustart = cputime()?;
        newtree();
        root = makecell()?;
        vector_zero(&mut (*root).cellnode.pos);
        expandbox(btab.as_mut_ptr(), nbody);
        load_all_bodies(btab.as_mut_ptr(), nbody)?;
        parse_options()?;
        reset_stats();
        hackcofm(root, rsize, 0)?;
        threadtree(root as *mut Node, std::ptr::null_mut());
        if usequad != 0 {
            hackquad(root);
        }
        cputree = (cputime()? - cpustart) as Real;
    }
    Ok(())
}

unsafe fn load_all_bodies(btab: *mut Body, nbody: i32) -> Result<()> {
    for i in 0..nbody as usize {
        loadbody(btab.add(i))?;
    }
    Ok(())
}

unsafe fn parse_options() -> Result<()> {
    let opts = getparam::getparam("options")?;
    BH86 = scanopt(&opts, "bh86");
    SW94 = scanopt(&opts, "sw94");
    if BH86 && SW94 {
        return Err(TreeError::IncompatibleOptions);
    }
    Ok(())
}

unsafe fn reset_stats() {
    tdepth = 0;
    for i in 0..MAXLEVEL {
        CELLHIST[i] = 0;
        SUBNHIST[i] = 0;
    }
}

unsafe fn newtree() {
    if !FIRSTCALL {
        let mut p = root as *mut Node;
        while !p.is_null() {
            if (*p).node_type == CELL {
                (*p).next = FREECELL;
                FREECELL = p;
                p = (*(p as *mut Cell)).more;
            } else {
                p = (*p).next;
            }
        }
    } else {
        FIRSTCALL = false;
    }
    root = std::ptr::null_mut();
    ncell = 0;
}

unsafe fn makecell() -> Result<*mut Cell> {
    let c: *mut Cell;
    if FREECELL.is_null() {
        c = allocate(std::mem::size_of::<Cell>())? as *mut Cell;
    } else {
        c = FREECELL as *mut Cell;
        FREECELL = (*FREECELL).next;
    }
    (*c).cellnode.node_type = CELL;
    (*c).cellnode.update = 0;
    for i in 0..NSUB {
        (*c).sorq.subp[i] = std::ptr::null_mut();
    }
    ncell += 1;
    Ok(c)
}

unsafe fn expandbox(btab: *mut Body, nbody: i32) {
    let mut dmax: Real = 0.0;
    for i in 0..nbody as usize {
        let p = &*btab.add(i);
        for k in 0..NDIM {
            let d = (p.bodynode.pos[k] - (*root).cellnode.pos[k]).abs();
            if d > dmax {
                dmax = d;
            }
        }
    }
    while rsize < 2.0 * dmax {
        rsize *= 2.0;
    }
}

unsafe fn loadbody(p: *mut Body) -> Result<()> {
    let mut q = root;
    let mut qind = subindex(p, q);
    let mut qsize = rsize;

    while !(*q).sorq.subp[qind].is_null() {
        if (*(*q).sorq.subp[qind]).node_type == BODY {
            let other = (*q).sorq.subp[qind] as *mut Body;
            require_distinct_positions(p, other)?;
            let c = makecell()?;
            set_cell_midpoint(c, q, p, qsize);
            let sub = subindex(other, c);
            (*c).sorq.subp[sub] = other as *mut Node;
            (*q).sorq.subp[qind] = c as *mut Node;
        }
        q = (*q).sorq.subp[qind] as *mut Cell;
        qind = subindex(p, q);
        qsize /= 2.0;
    }
    (*q).sorq.subp[qind] = p as *mut Node;
    Ok(())
}

unsafe fn require_distinct_positions(p: *mut Body, other: *mut Body) -> Result<()> {
    let mut dist2: Real = 0.0;
    for k in 0..NDIM {
        let d = (*p).bodynode.pos[k] - (*other).bodynode.pos[k];
        dist2 += d * d;
    }
    if dist2 == 0.0 {
        return Err(TreeError::CoincidentBodies);
    }
    Ok(())
}

unsafe fn set_cell_midpoint(c: *mut Cell, q: *mut Cell, p: *mut Body, qsize: Real) {
    for k in 0..NDIM {
        let offset = if (*p).bodynode.pos[k] < (*q).cellnode.pos[k] {
            -qsize
        } else {
            qsize
        } / 4.0;
        (*c).cellnode.pos[k] = (*q).cellnode.pos[k] + offset;
    }
}

unsafe fn subindex(p: *mut Body, q: *mut Cell) -> usize {
    let mut ind: usize = 0;
    for k in 0..NDIM {
        if (*q).cellnode.pos[k] <= (*p).bodynode.pos[k] {
            ind += NSUB >> (k + 1);
        }
    }
    ind
}

unsafe fn hackcofm(p: *mut Cell, psize: Real, lev: i32) -> Result<()> {
    let mut cmpos: Vector = Vector::zero();

    update_depth(lev);
    CELLHIST[lev as usize] += 1;
    (*p).cellnode.mass = 0.0;

    accumulate_subnodes(p, psize, lev, &mut cmpos)?;
    set_center_of_mass(p, &mut cmpos);
    verify_center(p, &cmpos, psize)?;
    setrcrit(p, &cmpos, psize);
    for k in 0..NDIM {
        (*p).cellnode.pos[k] = cmpos[k];
    }
    Ok(())
}

unsafe fn update_depth(lev: i32) {
    if lev > tdepth {
        tdepth = lev;
    }
}

unsafe fn accumulate_subnodes(
    p: *mut Cell,
    psize: Real,
    lev: i32,
    cmpos: &mut Vector,
) -> Result<()> {
    let mut tmpv: Vector = Vector::zero();
    for i in 0..NSUB {
        let q = (*p).sorq.subp[i];
        if !q.is_null() {
            SUBNHIST[lev as usize] += 1;
            if (*q).node_type == CELL {
                hackcofm(q as *mut Cell, psize / 2.0, lev + 1)?;
            }
            (*p).cellnode.update |= (*q).update;
            (*p).cellnode.mass += (*q).mass;
            for k in 0..NDIM {
                tmpv[k] = (*q).pos[k] * (*q).mass;
                cmpos[k] += tmpv[k];
            }
        }
    }
    Ok(())
}

unsafe fn set_center_of_mass(p: *mut Cell, cmpos: &mut Vector) {
    if (*p).cellnode.mass > 0.0 {
        for k in 0..NDIM {
            cmpos[k] /= (*p).cellnode.mass;
        }
    } else {
        for k in 0..NDIM {
            cmpos[k] = (*p).cellnode.pos[k];
        }
    }
}

unsafe fn verify_center(p: *mut Cell, cmpos: &Vector, psize: Real) -> Result<()> {
    for k in 0..NDIM {
        if cmpos[k] < (*p).cellnode.pos[k] - psize / 2.0
            || (*p).cellnode.pos[k] + psize / 2.0 <= cmpos[k]
        {
            return Err(TreeError::TreeStructure);
        }
    }
    Ok(())
}

unsafe fn setrcrit(p: *mut Cell, cmpos: &Vector, psize: Real) {
    if theta == 0.0 {
        set_rcrit_exact(p);
    } else if SW94 {
        set_rcrit_sw94(p, cmpos, psize);
    } else if BH86 {
        set_rcrit_bh86(p, psize);
    } else {
        set_rcrit_default(p, cmpos, psize);
    }
}

unsafe fn set_rcrit_exact(p: *mut Cell) {
    (*p).rcrit2 = mathfns::rsqr(2.0 * rsize);
}

unsafe fn set_rcrit_sw94(p: *mut Cell, cmpos: &Vector, psize: Real) {
    let mut bmax2: Real = 0.0;
    for k in 0..NDIM {
        let d = cmpos[k] - (*p).cellnode.pos[k] + psize / 2.0;
        bmax2 += mathfns::rsqr(d.max(psize - d));
    }
    (*p).rcrit2 = bmax2 / mathfns::rsqr(theta);
}

unsafe fn set_rcrit_bh86(p: *mut Cell, psize: Real) {
    (*p).rcrit2 = mathfns::rsqr(psize / theta);
}

unsafe fn set_rcrit_default(p: *mut Cell, cmpos: &Vector, psize: Real) {
    let mut d: Real = 0.0;
    for k in 0..NDIM {
        let dk = cmpos[k] - (*p).cellnode.pos[k];
        d += dk * dk;
    }
    (*p).rcrit2 = mathfns::rsqr(psize / theta + d.sqrt());
}

unsafe fn threadtree(p: *mut Node, n: *mut Node) {
    (*p).next = n;
    if (*p).node_type == CELL {
        let c = p as *mut Cell;
        let mut desc: [*mut Node; NSUB + 1] = [std::ptr::null_mut(); NSUB + 1];
        let ndesc = collect_descendants(c, &mut desc);
        (*c).more = desc[0];
        desc[ndesc] = n;
        for i in 0..ndesc {
            threadtree(desc[i], desc[i + 1]);
        }
    }
}

unsafe fn hackquad(p: *mut Cell) {
    let mut desc: [*mut Node; NSUB] = [std::ptr::null_mut(); NSUB];
    let ndesc = collect_descendants(p, &mut desc);

    matrix_zero(&mut (*p).sorq.quad);

    for i in 0..ndesc {
        let q = desc[i];
        if (*q).node_type == CELL {
            hackquad(q as *mut Cell);
        }
        accumulate_moment(p, q);
    }
}

unsafe fn collect_descendants(c: *mut Cell, desc: &mut [*mut Node]) -> usize {
    let mut ndesc: usize = 0;
    for i in 0..NSUB {
        if !(*c).sorq.subp[i].is_null() {
            desc[ndesc] = (*c).sorq.subp[i];
            ndesc += 1;
        }
    }
    ndesc
}

unsafe fn accumulate_moment(p: *mut Cell, q: *mut Node) {
    let dr = displacement(q, p);
    let mut tmpm = quadrupole_tensor(q, &dr);
    if (*q).node_type == CELL {
        let qm = &*(q as *mut Cell);
        for j in 0..NDIM {
            for k in 0..NDIM {
                tmpm[j][k] += qm.sorq.quad[j][k];
            }
        }
    }
    for j in 0..NDIM {
        for k in 0..NDIM {
            (*p).sorq.quad[j][k] += tmpm[j][k];
        }
    }
}

unsafe fn displacement(q: *mut Node, p: *mut Cell) -> Vector {
    let mut dr: Vector = Vector::zero();
    for k in 0..NDIM {
        dr[k] = (*q).pos[k] - (*p).cellnode.pos[k];
    }
    dr
}

unsafe fn quadrupole_tensor(q: *mut Node, dr: &Vector) -> Matrix {
    let drsq = dot_product(dr);
    let mut tmpm: Matrix = Matrix::zero();
    for j in 0..NDIM {
        for k in 0..NDIM {
            let id_rsq = if j == k { drsq } else { 0.0 };
            tmpm[j][k] = (3.0 * dr[j] * dr[k] - id_rsq) * (*q).mass;
        }
    }
    tmpm
}

fn dot_product(dr: &Vector) -> Real {
    let mut drsq: Real = 0.0;
    for k in 0..NDIM {
        drsq += dr[k] * dr[k];
    }
    drsq
}

pub fn tree_depth() -> i32 {
    unsafe { tdepth }
}

pub fn cell_count() -> i32 {
    unsafe { ncell }
}

pub fn tree_build_time() -> f64 {
    unsafe { cputree as f64 }
}
