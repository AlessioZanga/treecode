#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::many_single_char_names
)]

use crate::error::Result;
use crate::error::TreeError;
use crate::treecode::Tree;
use crate::types::{cputime, Interact, Matrix, NodeRef, Real, Vector, BODY, CELL, NDIM, NSUB};
use std::sync::Mutex;

const FACTIVE: Real = 0.75;

/// Counters accumulated during the force walk. They are exact integer
/// diagnostics (independent of the schedule), so they must match the C
/// reference bit-for-bit even when the walk is parallelized.
#[derive(Default)]
struct WalkCounters {
    actmax: i32,
    nbbcalc: i32,
    nbccalc: i32,
}

#[inline]
fn next_midpoint(pmid: Vector, pos: Vector, poff: Real) -> Vector {
    let mut nmid = Vector::zero();
    for k in 0..NDIM {
        let s = poff * (2.0 * usize::from(pos[k] >= pmid[k]) as Real - 1.0);
        nmid[k] = pmid[k] + s;
    }
    nmid
}

#[inline]
fn sumnode(
    eps: Real,
    interact: &[Interact],
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
        let phi_p = c.mass / drab;
        *phi0 -= phi_p;
        let mr3i = phi_p / dr2;
        add_mul_acc(acc0, &dr, mr3i);
    }
}

#[inline]
fn sumcell(
    eps: Real,
    interact: &[Interact],
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
        let phi_p = c.mass / drab;
        let mut mr3i = phi_p / dr2;
        let (qdr, drqdr) = quad_dot(c, &dr);
        let dr5i = 1.0 / (dr2 * dr2 * drab);
        let phi_q = 0.5 * dr5i * drqdr;
        *phi0 -= phi_p + phi_q;
        mr3i += 5.0 * phi_q / dr2;
        add_mul_acc2(acc0, &dr, mr3i, &qdr, -dr5i);
    }
}

#[inline]
fn separation(c: &Interact, pos0: Vector) -> (Vector, Real) {
    let mut dr = Vector::zero();
    let mut dr2: Real = 0.0;
    for k in 0..NDIM {
        dr[k] = c.pos[k] - pos0[k];
        dr2 += dr[k] * dr[k];
    }
    (dr, dr2)
}

#[inline]
fn quad_dot(c: &Interact, dr: &Vector) -> (Vector, Real) {
    let quad = c.quad.unwrap_or(Matrix::zero());
    let mut qdr = Vector::zero();
    let mut drqdr: Real = 0.0;
    for i in 0..NDIM {
        for j in 0..NDIM {
            qdr[i] += quad[i][j] * dr[j];
        }
        drqdr += qdr[i] * dr[i];
    }
    (qdr, drqdr)
}

#[inline(always)]
fn add_mul_acc(acc0: &mut Vector, dr: &Vector, s: Real) {
    for k in 0..NDIM {
        acc0[k] += dr[k] * s;
    }
}

#[inline(always)]
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
        let nb = self.nbody as usize;
        let results = Mutex::new(vec![(0.0 as Real, Vector::zero()); nb]);
        let counters = Mutex::new(WalkCounters::default());
        let cpustart = cputime()?;
        let root = self.root.ok_or(TreeError::TreeStructure)?;

        // Private scratch for the top-level walk. Each parallel subtree gets its
        // own copy, so the per-body arithmetic and summation order are identical
        // to the sequential C `gravcalc` -> byte-exact output.
        let mut active: Vec<NodeRef> = Vec::new();
        active
            .try_reserve(n)
            .map_err(|_| TreeError::OutOfMemory(n))?;
        active.resize(n, NodeRef::Body(0));
        let mut interact: Vec<Interact> = Vec::new();
        interact
            .try_reserve(n)
            .map_err(|_| TreeError::OutOfMemory(n))?;
        interact.resize(n, Interact::default());
        active[0] = NodeRef::Cell(root);

        self.walktree(
            &mut active,
            &mut interact,
            0,
            1,
            0,
            n,
            NodeRef::Cell(root),
            self.rsize,
            rmid,
            &results,
            &counters,
            true,
        )?;

        self.cpuforce = (cputime()? - cpustart) as Real;
        let c = counters.into_inner().unwrap_or_default();
        self.actmax = c.actmax;
        self.nbbcalc = c.nbbcalc;
        self.nbccalc = c.nbccalc;

        let res = results.into_inner().unwrap_or_default();
        for (b, (phi, acc)) in res.into_iter().enumerate() {
            self.bodytab[b].phi = phi;
            self.bodytab[b].acc = acc;
        }
        Ok(())
    }

    fn estimate_active_length(&self) -> i32 {
        let base = (FACTIVE * 216.0 * self.tdepth as Real) as i32;
        (base as Real * self.theta.powf(-2.5)) as i32
    }

    #[allow(clippy::too_many_arguments)]
    fn walktree(
        &self,
        active: &mut Vec<NodeRef>,
        interact: &mut Vec<Interact>,
        aptr: usize,
        nptr: usize,
        cptr: usize,
        bptr: usize,
        p: NodeRef,
        psize: Real,
        pmid: Vector,
        results: &Mutex<Vec<(Real, Vector)>>,
        counters: &Mutex<WalkCounters>,
        parallel: bool,
    ) -> Result<()> {
        let pnode = *self.node(p);
        if pnode.update != 0 {
            let mut np = nptr;
            let actsafe = self.actlen - NSUB as i32;
            let mut cptr = cptr;
            let mut bptr = bptr;
            let mut ap = aptr;
            while ap < nptr {
                let apnode = active[ap];
                let anode = *self.node(apnode);
                if anode.node_type == CELL {
                    if self.accept(apnode, psize, pmid) {
                        let cid = match apnode {
                            NodeRef::Cell(c) => c,
                            _ => unreachable!(),
                        };
                        let c = &self.cells[cid];
                        let mass = c.cellnode.mass;
                        let pos = c.cellnode.pos;
                        let quad = if self.usequad != 0 {
                            Some(c.sorq.quad())
                        } else {
                            None
                        };
                        interact[cptr].mass = mass;
                        interact[cptr].pos = pos;
                        interact[cptr].quad = quad;
                        cptr += 1;
                    } else {
                        if np as i32 >= actsafe {
                            return Err(TreeError::ActiveListOverflow);
                        }
                        let cid = match apnode {
                            NodeRef::Cell(c) => c,
                            _ => unreachable!(),
                        };
                        let pnext = anode.next;
                        let mut q = self.cells[cid].more;
                        while q != pnext {
                            let qr = q.ok_or(TreeError::TreeStructure)?;
                            active[np] = qr;
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
                    interact[bptr].mass = mass;
                    interact[bptr].pos = pos;
                }
                ap += 1;
            }
            let nact = np as i32;
            if let Ok(mut c) = counters.lock() {
                if nact > c.actmax {
                    c.actmax = nact;
                }
            }
            if np != nptr {
                self.walksub(
                    active, interact, nptr, np, cptr, bptr, p, psize, pmid, results, counters,
                    parallel,
                )?;
            } else {
                if pnode.node_type != BODY {
                    return Err(TreeError::RecursionTerminated);
                }
                self.gravsum(interact, p, cptr, bptr, results, counters);
            }
        }
        Ok(())
    }

    #[inline]
    fn accept(&self, c: NodeRef, psize: Real, pmid: Vector) -> bool {
        let cn = self.node(c);
        let mut dmax = psize;
        let mut dsq: Real = 0.0;
        for k in 0..NDIM {
            let dk = (cn.pos[k] - pmid[k]).abs();
            dmax = dmax.max(dk);
            let d = dk - 0.5 * psize;
            if d > 0.0 {
                dsq += d * d;
            }
        }
        let rcrit2 = match c {
            NodeRef::Cell(cid) => self.cells[cid].rcrit2,
            _ => unreachable!(),
        };
        dsq > rcrit2 && dmax > 1.5 * psize
    }

    #[allow(clippy::too_many_arguments)]
    fn walksub(
        &self,
        active: &mut Vec<NodeRef>,
        interact: &mut Vec<Interact>,
        nptr: usize,
        np: usize,
        cptr: usize,
        bptr: usize,
        p: NodeRef,
        psize: Real,
        pmid: Vector,
        results: &Mutex<Vec<(Real, Vector)>>,
        counters: &Mutex<WalkCounters>,
        parallel: bool,
    ) -> Result<()> {
        let poff = psize / 4.0;
        let pnode = *self.node(p);
        if let NodeRef::Cell(pid) = p {
            let pnext = pnode.next;
            let mut q = self.cells[pid].more;
            // Collect children first so the active/interact scratch can be moved
            // into each parallel task without borrow conflicts.
            let mut children = Vec::with_capacity(NSUB);
            while q != pnext {
                let qr = q.ok_or(TreeError::TreeStructure)?;
                children.push(qr);
                q = self.node(qr).next;
            }
            if parallel {
                let actlen = self.actlen as usize;
                let err_flag = Mutex::new(None::<TreeError>);
                std::thread::scope(|s| {
                    for qr in children {
                        let nmid = next_midpoint(pmid, self.node(qr).pos, poff);
                        let mut active2 = Vec::new();
                        if active2.try_reserve(actlen).is_err() {
                            if let Ok(mut g) = err_flag.lock() {
                                *g = Some(TreeError::OutOfMemory(actlen));
                            }
                            continue;
                        }
                        active2.resize(actlen, NodeRef::Body(0));
                        active2[0..(np - nptr)].copy_from_slice(&active[nptr..np]);
                        let mut interact2 = Vec::new();
                        if interact2.try_reserve(actlen).is_err() {
                            if let Ok(mut g) = err_flag.lock() {
                                *g = Some(TreeError::OutOfMemory(actlen));
                            }
                            continue;
                        }
                        interact2.resize(actlen, Interact::default());
                        interact2[0..cptr].copy_from_slice(&interact[0..cptr]);
                        interact2[bptr..actlen].copy_from_slice(&interact[bptr..actlen]);
                        let err_flag = &err_flag;
                        s.spawn(move || {
                            let r = self.walktree(
                                &mut active2,
                                &mut interact2,
                                0,
                                np - nptr,
                                cptr,
                                bptr,
                                qr,
                                psize / 2.0,
                                nmid,
                                results,
                                counters,
                                false,
                            );
                            if let Err(e) = r {
                                if let Ok(mut g) = err_flag.lock() {
                                    *g = Some(e);
                                }
                            }
                        });
                    }
                });
                if let Some(e) = err_flag.into_inner().unwrap_or_default() {
                    return Err(e);
                }
                Ok(())
            } else {
                for qr in children {
                    let nmid = next_midpoint(pmid, self.node(qr).pos, poff);
                    self.walktree(
                        active,
                        interact,
                        nptr,
                        np,
                        cptr,
                        bptr,
                        qr,
                        psize / 2.0,
                        nmid,
                        results,
                        counters,
                        false,
                    )?;
                }
                Ok(())
            }
        } else {
            let nmid = next_midpoint(pmid, pnode.pos, poff);
            self.walktree(
                active,
                interact,
                nptr,
                np,
                cptr,
                bptr,
                p,
                psize / 2.0,
                nmid,
                results,
                counters,
                false,
            )?;
            Ok(())
        }
    }

    fn gravsum(
        &self,
        interact: &[Interact],
        p0: NodeRef,
        cptr: usize,
        bptr: usize,
        results: &Mutex<Vec<(Real, Vector)>>,
        counters: &Mutex<WalkCounters>,
    ) {
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
            sumcell(eps, interact, 0, cptr, pos0, &mut phi0, &mut acc0);
        } else {
            sumnode(eps, interact, 0, cptr, pos0, &mut phi0, &mut acc0);
        }
        sumnode(eps, interact, bptr, actlen, pos0, &mut phi0, &mut acc0);
        if let NodeRef::Body(b) = p0 {
            if let Ok(mut g) = results.lock() {
                g[b] = (phi0, acc0);
            }
            if let Ok(mut c) = counters.lock() {
                c.nbbcalc += actlen as i32 - bptr as i32;
                c.nbccalc += cptr as i32;
            }
        }
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
