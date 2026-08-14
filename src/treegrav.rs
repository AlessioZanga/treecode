#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    // Raw-pointer tree access is temporary (removed by the Phase 4 arena);
    // indexing a `Deref` field through a raw-pointer deref is intentional here.
    dangerous_implicit_autorefs
)]

use crate::error::Result;
use crate::error::TreeError;
use crate::treecode::Tree;
use crate::types::{allocate, cputime, Body, Cell, Node, Real, Vector, BODY, CELL, NDIM, NSUB};

const FACTIVE: Real = 0.75;

fn next_midpoint(pmid: Vector, pos: Vector, poff: Real) -> Vector {
    let mut nmid = Vector::zero();
    for k in 0..NDIM {
        nmid[k] = pmid[k] + if pos[k] < pmid[k] { -poff } else { poff };
    }
    nmid
}

unsafe fn sumnode(
    eps: Real,
    start: *mut Cell,
    finish: *mut Cell,
    pos0: Vector,
    phi0: &mut Real,
    acc0: &mut Vector,
) {
    let eps2 = eps * eps;
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
    eps: Real,
    start: *mut Cell,
    finish: *mut Cell,
    pos0: Vector,
    phi0: &mut Real,
    acc0: &mut Vector,
) {
    let eps2 = eps * eps;
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

unsafe fn separation(p: *mut Cell, pos0: Vector) -> (Vector, Real) {
    let mut dr = Vector::zero();
    let mut dr2: Real = 0.0;
    for k in 0..NDIM {
        dr[k] = (*p).cellnode.pos[k] - pos0[k];
        dr2 += dr[k] * dr[k];
    }
    (dr, dr2)
}

unsafe fn quad_dot(p: *mut Cell, dr: &Vector) -> (Vector, Real) {
    let mut qdr = Vector::zero();
    let mut drqdr: Real = 0.0;
    for i in 0..NDIM {
        for j in 0..NDIM {
            qdr[i] += (*p).sorq.quad[i][j] * dr[j];
        }
        drqdr += qdr[i] * dr[i];
    }
    (qdr, drqdr)
}

fn add_mul_acc(acc0: &mut Vector, dr: &Vector, s: Real) {
    for k in 0..NDIM {
        acc0[k] += dr[k] * s;
    }
}

fn add_mul_acc2(acc0: &mut Vector, dr: &Vector, s: Real, w: &Vector, r: Real) {
    for k in 0..NDIM {
        acc0[k] += dr[k] * s + w[k] * r;
    }
}

impl Tree {
    pub fn gravcalc(&mut self) -> Result<()> {
        unsafe {
            let rmid: Vector = Vector::zero();

            self.actlen = self.estimate_active_length();
            self.active = allocate(self.actlen as usize * std::mem::size_of::<*mut Node>())?
                as *mut *mut Node;
            self.interact =
                allocate(self.actlen as usize * std::mem::size_of::<Cell>())? as *mut Cell;
            let cpustart = cputime()?;
            self.actmax = 0;
            self.nbbcalc = 0;
            self.nbccalc = 0;
            *self.active = self.root as *mut Node;
            self.walktree(
                self.active,
                self.active.add(1),
                self.interact,
                self.interact.add(self.actlen as usize),
                self.root as *mut Node,
                self.rsize,
                rmid,
            )?;
            self.cpuforce = (cputime()? - cpustart) as Real;
            libc::free(self.active as *mut libc::c_void);
            libc::free(self.interact as *mut libc::c_void);
            self.active = std::ptr::null_mut();
            self.interact = std::ptr::null_mut();
        }
        Ok(())
    }

    unsafe fn estimate_active_length(&self) -> i32 {
        let base = (FACTIVE * 216.0 * self.tdepth as Real) as i32;
        (base as Real * self.theta.powf(-2.5)) as i32
    }

    unsafe fn walktree(
        &mut self,
        aptr: *mut *mut Node,
        nptr: *mut *mut Node,
        cptr: *mut Cell,
        bptr: *mut Cell,
        p: *mut Node,
        psize: Real,
        pmid: Vector,
    ) -> Result<()> {
        if (*p).update != 0 {
            let mut np = nptr;
            let actsafe = self.actlen - NSUB as i32;
            let mut cptr = cptr;
            let mut bptr = bptr;
            let mut ap = aptr;
            while ap < nptr {
                let apnode = *ap;
                if (*apnode).node_type == CELL {
                    if Self::accept(apnode, psize, pmid) {
                        let acell = apnode as *mut Cell;
                        (*cptr).cellnode.mass = (*apnode).mass;
                        (*cptr).cellnode.pos = (*apnode).pos;
                        (*cptr).sorq.quad = (*acell).sorq.quad;
                        cptr = cptr.add(1);
                    } else {
                        if np.offset_from(self.active) >= actsafe as isize {
                            return Err(TreeError::ActiveListOverflow);
                        }
                        let mut q = (*(apnode as *mut Cell)).more;
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
            let nact = np.offset_from(self.active) as i32;
            if nact > self.actmax {
                self.actmax = nact;
            }
            if np != nptr {
                self.walksub(nptr, np, cptr, bptr, p, psize, pmid)?;
            } else {
                if (*p).node_type != BODY {
                    return Err(TreeError::RecursionTerminated);
                }
                self.gravsum(p as *mut Body, cptr, bptr);
            }
        }
        Ok(())
    }

    unsafe fn accept(c: *mut Node, psize: Real, pmid: Vector) -> bool {
        let mut dmax = psize;
        let mut dsq: Real = 0.0;
        for k in 0..NDIM {
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
        let rcrit2 = (*(c as *mut Cell)).rcrit2;
        dsq > rcrit2 && dmax > 1.5 * psize
    }

    unsafe fn walksub(
        &mut self,
        nptr: *mut *mut Node,
        np: *mut *mut Node,
        cptr: *mut Cell,
        bptr: *mut Cell,
        p: *mut Node,
        psize: Real,
        pmid: Vector,
    ) -> Result<()> {
        let poff = psize / 4.0;
        if (*p).node_type == CELL {
            let mut q = (*(p as *mut Cell)).more;
            while q != (*p).next {
                let nmid = next_midpoint(pmid, (*q).pos, poff);
                self.walktree(nptr, np, cptr, bptr, q, psize / 2.0, nmid)?;
                q = (*q).next;
            }
        } else {
            let nmid = next_midpoint(pmid, (*p).pos, poff);
            self.walktree(nptr, np, cptr, bptr, p, psize / 2.0, nmid)?;
        }
        Ok(())
    }

    unsafe fn gravsum(&mut self, p0: *mut Body, cptr: *mut Cell, bptr: *mut Cell) {
        let pos0 = (*p0).bodynode.pos;
        let mut phi0: Real = 0.0;
        let mut acc0: Vector = Vector::zero();
        if self.usequad != 0 {
            sumcell(self.eps, self.interact, cptr, pos0, &mut phi0, &mut acc0);
        } else {
            sumnode(self.eps, self.interact, cptr, pos0, &mut phi0, &mut acc0);
        }
        sumnode(
            self.eps,
            bptr,
            self.interact.add(self.actlen as usize),
            pos0,
            &mut phi0,
            &mut acc0,
        );
        (*p0).phi = phi0;
        (*p0).acc = acc0;
        self.nbbcalc += self.interact.add(self.actlen as usize).offset_from(bptr) as i32;
        self.nbccalc += cptr.offset_from(self.interact) as i32;
    }

    pub fn force_max_active(&self) -> i32 {
        self.actmax
    }

    pub fn force_bb_calc(&self) -> i32 {
        self.nbbcalc
    }

    pub fn force_bc_calc(&self) -> i32 {
        self.nbccalc
    }

    pub fn force_cpu_time(&self) -> f64 {
        self.cpuforce as f64
    }
}
