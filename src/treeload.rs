#![allow(
    clippy::needless_range_loop,
    clippy::manual_memcpy,
    // Raw-pointer tree access is temporary (removed by the Phase 4 arena);
    // indexing a `Deref` field through a raw-pointer deref is intentional here.
    dangerous_implicit_autorefs
)]

use crate::error::Result;
use crate::error::TreeError;
use crate::mathfns;
use crate::treecode::{Tree, MAXLEVEL};
use crate::types::{allocate, cputime, Body, Cell, Node, Real, Vector, BODY, CELL, NDIM, NSUB};
use crate::vecmath::{matrix_zero, vector_zero, Matrix};

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

unsafe fn set_rcrit_exact(p: *mut Cell, rsize: Real) {
    (*p).rcrit2 = mathfns::rsqr(2.0 * rsize);
}

unsafe fn set_rcrit_sw94(p: *mut Cell, cmpos: &Vector, psize: Real, theta: Real) {
    let mut bmax2: Real = 0.0;
    for k in 0..NDIM {
        let d = cmpos[k] - (*p).cellnode.pos[k] + psize / 2.0;
        bmax2 += mathfns::rsqr(d.max(psize - d));
    }
    (*p).rcrit2 = bmax2 / mathfns::rsqr(theta);
}

unsafe fn set_rcrit_bh86(p: *mut Cell, psize: Real, theta: Real) {
    (*p).rcrit2 = mathfns::rsqr(psize / theta);
}

unsafe fn set_rcrit_default(p: *mut Cell, cmpos: &Vector, psize: Real, theta: Real) {
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

impl Tree {
    pub fn maketree(&mut self, nbody: i32) -> Result<()> {
        unsafe {
            let cpustart = cputime()?;
            let btab = self.bodytab.as_mut_ptr();
            self.newtree();
            self.root = self.makecell()?;
            vector_zero(&mut (*self.root).cellnode.pos);
            self.expandbox(btab, nbody);
            self.load_all_bodies(btab, nbody)?;
            self.parse_options()?;
            self.reset_stats();
            self.hackcofm(self.root, self.rsize, 0)?;
            threadtree(self.root as *mut Node, std::ptr::null_mut());
            if self.usequad != 0 {
                hackquad(self.root);
            }
            self.cputree = (cputime()? - cpustart) as Real;
        }
        Ok(())
    }

    unsafe fn load_all_bodies(&mut self, btab: *mut Body, nbody: i32) -> Result<()> {
        for i in 0..nbody as usize {
            self.loadbody(btab.add(i))?;
        }
        Ok(())
    }

    unsafe fn parse_options(&mut self) -> Result<()> {
        self.bh86 = crate::types::scanopt(&self.options, "bh86");
        self.sw94 = crate::types::scanopt(&self.options, "sw94");
        if self.bh86 && self.sw94 {
            return Err(TreeError::IncompatibleOptions);
        }
        Ok(())
    }

    unsafe fn reset_stats(&mut self) {
        self.tdepth = 0;
        for i in 0..MAXLEVEL {
            self.cellhist[i] = 0;
            self.subnhist[i] = 0;
        }
    }

    unsafe fn newtree(&mut self) {
        if !self.firstcall {
            let mut p = self.root as *mut Node;
            while !p.is_null() {
                if (*p).node_type == CELL {
                    (*p).next = self.freecell;
                    self.freecell = p;
                    p = (*(p as *mut Cell)).more;
                } else {
                    p = (*p).next;
                }
            }
        } else {
            self.firstcall = false;
        }
        self.root = std::ptr::null_mut();
        self.ncell = 0;
    }

    unsafe fn makecell(&mut self) -> Result<*mut Cell> {
        let c: *mut Cell;
        if self.freecell.is_null() {
            c = allocate(std::mem::size_of::<Cell>())? as *mut Cell;
        } else {
            c = self.freecell as *mut Cell;
            self.freecell = (*self.freecell).next;
        }
        (*c).cellnode.node_type = CELL;
        (*c).cellnode.update = 0;
        for i in 0..NSUB {
            (*c).sorq.subp[i] = std::ptr::null_mut();
        }
        self.ncell += 1;
        Ok(c)
    }

    unsafe fn expandbox(&mut self, btab: *mut Body, nbody: i32) {
        let mut dmax: Real = 0.0;
        for i in 0..nbody as usize {
            let p = &*btab.add(i);
            for k in 0..NDIM {
                let d = (p.bodynode.pos[k] - (*self.root).cellnode.pos[k]).abs();
                if d > dmax {
                    dmax = d;
                }
            }
        }
        while self.rsize < 2.0 * dmax {
            self.rsize *= 2.0;
        }
    }

    unsafe fn loadbody(&mut self, p: *mut Body) -> Result<()> {
        let mut q = self.root;
        let mut qind = subindex(p, q);
        let mut qsize = self.rsize;

        while !(*q).sorq.subp[qind].is_null() {
            if (*(*q).sorq.subp[qind]).node_type == BODY {
                let other = (*q).sorq.subp[qind] as *mut Body;
                require_distinct_positions(p, other)?;
                let c = self.makecell()?;
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

    unsafe fn hackcofm(&mut self, p: *mut Cell, psize: Real, lev: i32) -> Result<()> {
        let mut cmpos: Vector = Vector::zero();

        self.update_depth(lev);
        self.cellhist[lev as usize] += 1;
        (*p).cellnode.mass = 0.0;

        self.accumulate_subnodes(p, psize, lev, &mut cmpos)?;
        set_center_of_mass(p, &mut cmpos);
        verify_center(p, &cmpos, psize)?;
        self.setrcrit(p, &cmpos, psize);
        for k in 0..NDIM {
            (*p).cellnode.pos[k] = cmpos[k];
        }
        Ok(())
    }

    unsafe fn update_depth(&mut self, lev: i32) {
        if lev > self.tdepth {
            self.tdepth = lev;
        }
    }

    unsafe fn accumulate_subnodes(
        &mut self,
        p: *mut Cell,
        psize: Real,
        lev: i32,
        cmpos: &mut Vector,
    ) -> Result<()> {
        let mut tmpv: Vector = Vector::zero();
        for i in 0..NSUB {
            let q = (*p).sorq.subp[i];
            if !q.is_null() {
                self.subnhist[lev as usize] += 1;
                if (*q).node_type == CELL {
                    self.hackcofm(q as *mut Cell, psize / 2.0, lev + 1)?;
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

    unsafe fn setrcrit(&self, p: *mut Cell, cmpos: &Vector, psize: Real) {
        if self.theta == 0.0 {
            set_rcrit_exact(p, self.rsize);
        } else if self.sw94 {
            set_rcrit_sw94(p, cmpos, psize, self.theta);
        } else if self.bh86 {
            set_rcrit_bh86(p, psize, self.theta);
        } else {
            set_rcrit_default(p, cmpos, psize, self.theta);
        }
    }

    pub fn tree_depth(&self) -> i32 {
        self.tdepth
    }

    pub fn cell_count(&self) -> i32 {
        self.ncell
    }

    pub fn tree_build_time(&self) -> f64 {
        self.cputree as f64
    }
}
