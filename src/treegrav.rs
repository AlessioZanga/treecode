#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::many_single_char_names
)]

use crate::error::Result;
use crate::error::TreeError;
use crate::treecode::Tree;
use crate::types::{cputime, Cell, NodeRef, Real, Sorq, Vector, BODY, CELL, NDIM, NSUB};

const FACTIVE: Real = 0.75;

fn next_midpoint(pmid: Vector, pos: Vector, poff: Real) -> Vector {
    let mut nmid = Vector::zero();
    for k in 0..NDIM {
        nmid[k] = pmid[k] + if pos[k] < pmid[k] { -poff } else { poff };
    }
    nmid
}

fn sumnode(
    eps: Real,
    interact: &[Cell],
    start: usize,
    finish: usize,
    pos0: Vector,
    phi0: &mut Real,
    acc0: &mut Vector,
) {
    let eps2 = eps * eps;
    for idx in start..finish {
        let c = &interact[idx];
        let (dr, mut dr2) = separation(c, pos0);
        dr2 += eps2;
        let drab = dr2.sqrt();
        let phi_p = c.cellnode.mass / drab;
        *phi0 -= phi_p;
        let mr3i = phi_p / dr2;
        add_mul_acc(acc0, &dr, mr3i);
    }
}

fn sumcell(
    eps: Real,
    interact: &[Cell],
    start: usize,
    finish: usize,
    pos0: Vector,
    phi0: &mut Real,
    acc0: &mut Vector,
) {
    let eps2 = eps * eps;
    for idx in start..finish {
        let c = &interact[idx];
        let (dr, mut dr2) = separation(c, pos0);
        dr2 += eps2;
        let drab = dr2.sqrt();
        let phi_p = c.cellnode.mass / drab;
        let mut mr3i = phi_p / dr2;
        let (qdr, drqdr) = quad_dot(c, &dr);
        let dr5i = 1.0 / (dr2 * dr2 * drab);
        let phi_q = 0.5 * dr5i * drqdr;
        *phi0 -= phi_p + phi_q;
        mr3i += 5.0 * phi_q / dr2;
        add_mul_acc2(acc0, &dr, mr3i, &qdr, -dr5i);
    }
}

fn separation(c: &Cell, pos0: Vector) -> (Vector, Real) {
    let mut dr = Vector::zero();
    let mut dr2: Real = 0.0;
    for k in 0..NDIM {
        dr[k] = c.cellnode.pos[k] - pos0[k];
        dr2 += dr[k] * dr[k];
    }
    (dr, dr2)
}

fn quad_dot(c: &Cell, dr: &Vector) -> (Vector, Real) {
    let mut qdr = Vector::zero();
    let mut drqdr: Real = 0.0;
    for i in 0..NDIM {
        for j in 0..NDIM {
            qdr[i] += c.sorq.quad()[i][j] * dr[j];
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
        let rmid: Vector = Vector::zero();

        self.actlen = self.estimate_active_length();
        let n = if self.actlen > 0 {
            self.actlen as usize
        } else {
            1
        };
        let mut active: Vec<NodeRef> = Vec::new();
        active
            .try_reserve(n)
            .map_err(|_| TreeError::OutOfMemory(n))?;
        active.resize(n, NodeRef::Body(0));
        let mut interact: Vec<Cell> = Vec::new();
        interact
            .try_reserve(n)
            .map_err(|_| TreeError::OutOfMemory(n))?;
        interact.resize(n, Cell::default());
        let cpustart = cputime()?;
        self.actmax = 0;
        self.nbbcalc = 0;
        self.nbccalc = 0;
        let root = self.root.ok_or(TreeError::TreeStructure)?;
        active[0] = NodeRef::Cell(root);
        self.active = active;
        self.interact = interact;
        self.walktree(
            0,
            1,
            0,
            self.actlen as usize,
            NodeRef::Cell(root),
            self.rsize,
            rmid,
        )?;
        self.cpuforce = (cputime()? - cpustart) as Real;
        self.active = Vec::new();
        self.interact = Vec::new();
        Ok(())
    }

    fn estimate_active_length(&self) -> i32 {
        let base = (FACTIVE * 216.0 * self.tdepth as Real) as i32;
        (base as Real * self.theta.powf(-2.5)) as i32
    }

    fn walktree(
        &mut self,
        aptr: usize,
        nptr: usize,
        cptr: usize,
        bptr: usize,
        p: NodeRef,
        psize: Real,
        pmid: Vector,
    ) -> Result<()> {
        if self.node(p).update != 0 {
            let mut np = nptr;
            let actsafe = self.actlen - NSUB as i32;
            let mut cptr = cptr;
            let mut bptr = bptr;
            let mut ap = aptr;
            while ap < nptr {
                let apnode = self.active[ap];
                if self.node(apnode).node_type == CELL {
                    if self.accept(apnode, psize, pmid) {
                        let cid = match apnode {
                            NodeRef::Cell(c) => c,
                            _ => unreachable!(),
                        };
                        let c = &self.cells[cid];
                        let mass = c.cellnode.mass;
                        let pos = c.cellnode.pos;
                        let sq = if self.usequad != 0 {
                            Sorq::Quad(c.sorq.quad())
                        } else {
                            Sorq::Subp([None; NSUB])
                        };
                        self.interact[cptr].cellnode.mass = mass;
                        self.interact[cptr].cellnode.pos = pos;
                        self.interact[cptr].sorq = sq;
                        cptr += 1;
                    } else {
                        if np as i32 >= actsafe {
                            return Err(TreeError::ActiveListOverflow);
                        }
                        let cid = match apnode {
                            NodeRef::Cell(c) => c,
                            _ => unreachable!(),
                        };
                        let pnext = self.node(apnode).next;
                        let mut q = self.cells[cid].more;
                        while q != pnext {
                            let qr = q.ok_or(TreeError::TreeStructure)?;
                            self.active[np] = qr;
                            np += 1;
                            q = self.node(qr).next;
                        }
                    }
                } else if apnode != p {
                    bptr -= 1;
                    let (mass, pos) = match apnode {
                        NodeRef::Body(b) => {
                            (self.bodytab[b].bodynode.mass, self.bodytab[b].bodynode.pos)
                        }
                        _ => unreachable!(),
                    };
                    self.interact[bptr].cellnode.mass = mass;
                    self.interact[bptr].cellnode.pos = pos;
                }
                ap += 1;
            }
            let nact = np as i32;
            if nact > self.actmax {
                self.actmax = nact;
            }
            if np != nptr {
                self.walksub(nptr, np, cptr, bptr, p, psize, pmid)?;
            } else {
                if self.node(p).node_type != BODY {
                    return Err(TreeError::RecursionTerminated);
                }
                self.gravsum(p, cptr, bptr);
            }
        }
        Ok(())
    }

    fn accept(&self, c: NodeRef, psize: Real, pmid: Vector) -> bool {
        let cn = self.node(c);
        let mut dmax = psize;
        let mut dsq: Real = 0.0;
        for k in 0..NDIM {
            let mut dk = cn.pos[k] - pmid[k];
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
        let rcrit2 = match c {
            NodeRef::Cell(cid) => self.cells[cid].rcrit2,
            _ => unreachable!(),
        };
        dsq > rcrit2 && dmax > 1.5 * psize
    }

    fn walksub(
        &mut self,
        nptr: usize,
        np: usize,
        cptr: usize,
        bptr: usize,
        p: NodeRef,
        psize: Real,
        pmid: Vector,
    ) -> Result<()> {
        let poff = psize / 4.0;
        if let NodeRef::Cell(pid) = p {
            let pnext = self.node(p).next;
            let mut q = self.cells[pid].more;
            while q != pnext {
                let qr = q.ok_or(TreeError::TreeStructure)?;
                let nmid = next_midpoint(pmid, self.node(qr).pos, poff);
                self.walktree(nptr, np, cptr, bptr, qr, psize / 2.0, nmid)?;
                q = self.node(qr).next;
            }
        } else {
            let nmid = next_midpoint(pmid, self.node(p).pos, poff);
            self.walktree(nptr, np, cptr, bptr, p, psize / 2.0, nmid)?;
        }
        Ok(())
    }

    fn gravsum(&mut self, p0: NodeRef, cptr: usize, bptr: usize) {
        let pos0 = match p0 {
            NodeRef::Body(b) => self.bodytab[b].bodynode.pos,
            _ => unreachable!(),
        };
        let eps = self.eps;
        let usequad = self.usequad != 0;
        let actlen = self.actlen as usize;
        let mut phi0: Real = 0.0;
        let mut acc0: Vector = Vector::zero();
        if usequad {
            sumcell(eps, &self.interact, 0, cptr, pos0, &mut phi0, &mut acc0);
        } else {
            sumnode(eps, &self.interact, 0, cptr, pos0, &mut phi0, &mut acc0);
        }
        sumnode(
            eps,
            &self.interact,
            bptr,
            actlen,
            pos0,
            &mut phi0,
            &mut acc0,
        );
        if let NodeRef::Body(b) = p0 {
            self.bodytab[b].phi = phi0;
            self.bodytab[b].acc = acc0;
        }
        self.nbbcalc += actlen as i32 - bptr as i32;
        self.nbccalc += cptr as i32;
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
