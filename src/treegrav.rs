use std::sync::Mutex;

#[cfg(feature = "simd")]
use wide::f32x8;

use crate::{
    error::{Result, TreeError},
    treecode::Tree,
    types::{BODY, CELL, Interact, Matrix, NDIM, NSUB, NodeRef, Vector, cputime},
};

const FACTIVE: f32 = 0.75;

// Number of tree levels to parallelize via rayon fan-out. Each fanned-out cell
// spawns its children as rayon tasks; cells deeper than this walk sequentially.
// Bounded so task count stays ~8^PARALLEL_FANOUT_DEPTH instead of exploding.
const PARALLEL_FANOUT_DEPTH: u8 = 3;

/// Counters accumulated during the force walk. They are exact integer
/// diagnostics (independent of the schedule), so they must match the C
/// reference bit-for-bit even when the walk is parallelized.
#[derive(Default, Clone, Copy)]
struct WalkCounters {
    actmax: usize,
    nbbcalc: usize,
    nbccalc: usize,
}

/// Mutable, per-level state for the force walk, bundled into one object so every
/// walk function takes a single context argument (keeps arity under the clippy
/// `too_many_arguments` limit and the parallel fan-out closures short/flat).
struct WalkContext<'a> {
    active: &'a mut Vec<NodeRef>,
    interact: &'a mut Vec<Interact>,
    aptr: usize,
    nptr: usize,
    cptr: usize,
    bptr: usize,
    np: usize,
    p: NodeRef,
    psize: f32,
    pmid: Vector,
    results: &'a Mutex<Vec<(f32, Vector)>>,
    counters: &'a Mutex<WalkCounters>,
    parallel_depth: u8,
}

#[inline]
fn next_midpoint(pmid: Vector, pos: Vector, poff: f32) -> Vector {
    let mut nmid = Vector::zero();
    for k in 0..NDIM {
        let s = poff * (2.0 * usize::from(pos[k] >= pmid[k]) as f32 - 1.0);
        nmid[k] = pmid[k] + s;
    }
    nmid
}

#[inline]
fn set_result(err_flag: &Mutex<Option<TreeError>>, r: Result<()>) {
    if let Err(e) = r
        && let Ok(mut g) = err_flag.lock()
    {
        *g = Some(e);
    }
}

#[inline]
fn sumnode(
    eps2: f32,
    interact: &[Interact],
    start: usize,
    finish: usize,
    pos0: Vector,
    phi0: &mut f32,
    acc0: &mut Vector,
) {
    interact[start..finish].iter().for_each(|c| {
        let (dr, mut dr2) = separation(c, pos0);
        dr2 += eps2;
        let drab = dr2.sqrt();
        let phi_p = c.mass / drab;
        *phi0 -= phi_p;
        let mr3i = phi_p / dr2;
        add_mul_acc(acc0, &dr, mr3i);
    });
}

#[inline]
fn sumcell(
    eps2: f32,
    interact: &[Interact],
    start: usize,
    finish: usize,
    pos0: Vector,
    phi0: &mut f32,
    acc0: &mut Vector,
) {
    interact[start..finish].iter().for_each(|c| {
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
    });
}

#[inline]
fn separation(c: &Interact, pos0: Vector) -> (Vector, f32) {
    let mut dr = Vector::zero();
    let mut dr2: f32 = 0.0;
    for k in 0..NDIM {
        dr[k] = c.pos[k] - pos0[k];
        dr2 += dr[k] * dr[k];
    }
    (dr, dr2)
}

#[inline]
fn quad_dot(c: &Interact, dr: &Vector) -> (Vector, f32) {
    let quad = c.quad;
    let mut qdr = Vector::zero();
    let mut drqdr: f32 = 0.0;
    for i in 0..NDIM {
        for j in 0..NDIM {
            qdr[i] += quad[i][j] * dr[j];
        }
        drqdr += qdr[i] * dr[i];
    }
    (qdr, drqdr)
}

#[inline(always)]
fn add_mul_acc(acc0: &mut Vector, dr: &Vector, s: f32) {
    for k in 0..NDIM {
        acc0[k] += dr[k] * s;
    }
}

#[inline(always)]
fn add_mul_acc2(acc0: &mut Vector, dr: &Vector, s: f32, w: &Vector, r: f32) {
    for k in 0..NDIM {
        acc0[k] += dr[k] * s + w[k] * r;
    }
}

/// Number of interactions processed per SIMD vector in the `simd` feature
/// kernels. The reduction is folded back into the running `f32` accumulators in
/// exact list order, so the lane width is irrelevant to byte-exactness and is
/// only a throughput knob.
#[cfg(feature = "simd")]
const SIMD_LANES: usize = 8;

/// SIMD force accumulation for the monopole (node) interaction list.
///
/// The per-interaction arithmetic (`sqrt`, divide, `dr*mr3i`) is an independent
/// *map* across the 8 lanes, so it is safe to vectorize. The fold back into
/// `phi0`/`acc0` is done one lane at a time in list order (`i = 0..LANES`),
/// exactly mirroring the scalar `sumnode` reduction — no reassociation, so the
/// result is bit-identical to the scalar path. FMA is avoided (separate mul/add)
/// to match the `-fma`-disabled SSE2 C reference.
#[cfg(feature = "simd")]
#[inline]
fn sumnode_simd(
    eps2: f32,
    interact: &[Interact],
    start: usize,
    finish: usize,
    pos0: Vector,
    phi0: &mut f32,
    acc0: &mut Vector,
) {
    const L: usize = SIMD_LANES;
    type F = f32x8;
    let px = pos0[0];
    let py = pos0[1];
    let pz = pos0[2];
    let mut base = start;
    while base + L <= finish {
        let mut mass = [0.0f32; L];
        let mut dx = [0.0f32; L];
        let mut dy = [0.0f32; L];
        let mut dz = [0.0f32; L];
        for i in 0..L {
            let c = &interact[base + i];
            mass[i] = c.mass;
            dx[i] = c.pos[0] - px;
            dy[i] = c.pos[1] - py;
            dz[i] = c.pos[2] - pz;
        }
        let vdx = F::from(dx);
        let vdy = F::from(dy);
        let vdz = F::from(dz);
        let vmass = F::from(mass);
        let mut vdr2 = F::splat(0.0);
        vdr2 += vdx * vdx;
        vdr2 += vdy * vdy;
        vdr2 += vdz * vdz;
        vdr2 += F::splat(eps2);
        let vdrab = vdr2.sqrt();
        let vphi = vmass / vdrab;
        let vmr3i = vphi / vdr2;
        let pa = vphi.to_array();
        let ma = vmr3i.to_array();
        let dxa = vdx.to_array();
        let dya = vdy.to_array();
        let dza = vdz.to_array();
        for i in 0..L {
            *phi0 -= pa[i];
            acc0[0] += dxa[i] * ma[i];
            acc0[1] += dya[i] * ma[i];
            acc0[2] += dza[i] * ma[i];
        }
        base += L;
    }
    // Remainder stays on the scalar, byte-exact path.
    sumnode(eps2, interact, base, finish, pos0, phi0, acc0);
}

/// SIMD force accumulation for the quadrupole (cell) interaction list.
///
/// Same strategy as [`sumnode_simd`]: vectorize the independent per-interaction
/// math (including the `quad_dot` matrix-vector product and the `drqdr`
/// contraction) and fold back in exact list order. The quadrupole tensor is
/// loaded lane-by-lane into `f32` arrays and reconstructed as vectors; no
/// reassociation is introduced, so output matches the scalar `sumcell`.
#[cfg(feature = "simd")]
#[inline]
fn sumcell_simd(
    eps2: f32,
    interact: &[Interact],
    start: usize,
    finish: usize,
    pos0: Vector,
    phi0: &mut f32,
    acc0: &mut Vector,
) {
    const L: usize = SIMD_LANES;
    type F = f32x8;
    let px = pos0[0];
    let py = pos0[1];
    let pz = pos0[2];
    let mut base = start;
    while base + L <= finish {
        let mut mass = [0.0f32; L];
        let mut dx = [0.0f32; L];
        let mut dy = [0.0f32; L];
        let mut dz = [0.0f32; L];
        let mut qxx = [0.0f32; L];
        let mut qxy = [0.0f32; L];
        let mut qxz = [0.0f32; L];
        let mut qyx = [0.0f32; L];
        let mut qyy = [0.0f32; L];
        let mut qyz = [0.0f32; L];
        let mut qzx = [0.0f32; L];
        let mut qzy = [0.0f32; L];
        let mut qzz = [0.0f32; L];
        for i in 0..L {
            let c = &interact[base + i];
            mass[i] = c.mass;
            dx[i] = c.pos[0] - px;
            dy[i] = c.pos[1] - py;
            dz[i] = c.pos[2] - pz;
            let q = c.quad;
            qxx[i] = q[0][0];
            qxy[i] = q[0][1];
            qxz[i] = q[0][2];
            qyx[i] = q[1][0];
            qyy[i] = q[1][1];
            qyz[i] = q[1][2];
            qzx[i] = q[2][0];
            qzy[i] = q[2][1];
            qzz[i] = q[2][2];
        }
        let vdx = F::from(dx);
        let vdy = F::from(dy);
        let vdz = F::from(dz);
        let vmass = F::from(mass);
        let vqxx = F::from(qxx);
        let vqxy = F::from(qxy);
        let vqxz = F::from(qxz);
        let vqyx = F::from(qyx);
        let vqyy = F::from(qyy);
        let vqyz = F::from(qyz);
        let vqzx = F::from(qzx);
        let vqzy = F::from(qzy);
        let vqzz = F::from(qzz);
        let mut vdr2 = F::splat(0.0);
        vdr2 += vdx * vdx;
        vdr2 += vdy * vdy;
        vdr2 += vdz * vdz;
        vdr2 += F::splat(eps2);
        let vdrab = vdr2.sqrt();
        let vphi = vmass / vdrab;
        let mut vmr3i = vphi / vdr2;
        let vqdrx = vqxx * vdx + vqxy * vdy + vqxz * vdz;
        let vqdry = vqyx * vdx + vqyy * vdy + vqyz * vdz;
        let vqdrz = vqzx * vdx + vqzy * vdy + vqzz * vdz;
        let vdrqdr = vqdrx * vdx + vqdry * vdy + vqdrz * vdz;
        let vdr5i = F::splat(1.0) / (vdr2 * vdr2 * vdrab);
        let vphi_q = F::splat(0.5) * vdr5i * vdrqdr;
        let vphi_tot = vphi + vphi_q;
        vmr3i += F::splat(5.0) * vphi_q / vdr2;
        let vnegdr5i = F::splat(0.0) - vdr5i;
        let pat = vphi_tot.to_array();
        let mat = vmr3i.to_array();
        let dxa = vdx.to_array();
        let dya = vdy.to_array();
        let dza = vdz.to_array();
        let qdrxa = vqdrx.to_array();
        let qdrya = vqdry.to_array();
        let qdrza = vqdrz.to_array();
        let nd5a = vnegdr5i.to_array();
        for i in 0..L {
            *phi0 -= pat[i];
            acc0[0] += dxa[i] * mat[i] + qdrxa[i] * nd5a[i];
            acc0[1] += dya[i] * mat[i] + qdrya[i] * nd5a[i];
            acc0[2] += dza[i] * mat[i] + qdrza[i] * nd5a[i];
        }
        base += L;
    }
    sumcell(eps2, interact, base, finish, pos0, phi0, acc0);
}

impl Tree {
    pub fn gravcalc(&mut self) -> Result<()> {
        let rmid: Vector = Vector::zero();

        self.actlen = self.estimate_active_length();
        let n = if self.actlen > 0 { self.actlen } else { 1 };
        let nb = self.nbody;
        let results = Mutex::new(vec![(0.0, Vector::zero()); nb]);
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

        let mut ctx = WalkContext {
            active: &mut active,
            interact: &mut interact,
            aptr: 0,
            nptr: 1,
            cptr: 0,
            bptr: n,
            np: 0,
            p: NodeRef::Cell(root),
            psize: self.rsize,
            pmid: rmid,
            results: &results,
            counters: &counters,
            parallel_depth: PARALLEL_FANOUT_DEPTH,
        };
        self.walktree(&mut ctx)?;

        self.cpuforce = (cputime()? - cpustart) as f32;
        let c = counters.into_inner().unwrap_or_default();
        self.actmax = c.actmax;
        self.nbbcalc = c.nbbcalc;
        self.nbccalc = c.nbccalc;

        let res = results.into_inner().unwrap_or_default();
        for (p, (phi, acc)) in self.bodytab.iter_mut().zip(res) {
            p.phi = phi;
            p.acc = acc;
        }
        Ok(())
    }

    fn estimate_active_length(&self) -> usize {
        let base = (FACTIVE * 216.0 * self.tdepth as f32) as usize;
        (base as f32 * self.theta_pow_m2_5) as usize
    }

    #[inline]
    fn walktree(&self, ctx: &mut WalkContext) -> Result<()> {
        let pnode = *self.node(ctx.p);
        if pnode.update == 0 {
            return Ok(());
        }
        let mut np = ctx.nptr;
        let actsafe = self.actlen - NSUB;
        let mut cptr = ctx.cptr;
        let mut bptr = ctx.bptr;
        let mut ap = ctx.aptr;
        while ap < ctx.nptr {
            let apnode = ctx.active[ap];
            let anode = *self.node(apnode);
            if anode.node_type == CELL {
                self.process_cell_node(ctx, apnode, actsafe, &mut np, &mut cptr)?;
            } else if apnode != ctx.p {
                self.process_body_node(ctx, apnode, &mut bptr);
            }
            ap += 1;
        }
        if let Ok(mut c) = ctx.counters.lock()
            && np > c.actmax
        {
            c.actmax = np;
        }
        if np != ctx.nptr {
            ctx.np = np;
            ctx.cptr = cptr;
            ctx.bptr = bptr;
            self.walksub(ctx)?;
        } else if pnode.node_type != BODY {
            return Err(TreeError::RecursionTerminated);
        } else {
            self.gravsum(ctx.interact, ctx.p, cptr, bptr, ctx.results, ctx.counters);
        }
        Ok(())
    }

    #[inline]
    fn process_cell_node(
        &self,
        ctx: &mut WalkContext,
        apnode: NodeRef,
        actsafe: usize,
        np: &mut usize,
        cptr: &mut usize,
    ) -> Result<()> {
        let cid = match apnode {
            NodeRef::Cell(c) => c,
            _ => unreachable!(),
        };
        if self.accept(apnode, ctx.psize, ctx.pmid) {
            let c = &self.cells[cid];
            let mass = c.cellnode.mass;
            let pos = c.cellnode.pos;
            let quad = if self.usequad != 0 {
                c.sorq.quad()
            } else {
                Matrix::zero()
            };
            ctx.interact[*cptr].mass = mass;
            ctx.interact[*cptr].pos = pos;
            ctx.interact[*cptr].quad = quad;
            *cptr += 1;
        } else {
            if *np >= actsafe {
                return Err(TreeError::ActiveListOverflow);
            }
            let pnext = self.node(apnode).next;
            let mut q = self.cells[cid].more;
            while q != pnext {
                if q.is_none() {
                    return Err(TreeError::TreeStructure);
                }
                ctx.active[*np] = q;
                *np += 1;
                q = self.node(q).next;
            }
        }
        Ok(())
    }

    #[inline]
    fn process_body_node(&self, ctx: &mut WalkContext, apnode: NodeRef, bptr: &mut usize) {
        *bptr -= 1;
        let (mass, pos) = match apnode {
            NodeRef::Body(b) => (self.bodytab[b].bodynode.mass, self.bodytab[b].bodynode.pos),
            _ => unreachable!(),
        };
        ctx.interact[*bptr].mass = mass;
        ctx.interact[*bptr].pos = pos;
    }

    #[inline]
    fn accept(&self, c: NodeRef, psize: f32, pmid: Vector) -> bool {
        let cn = self.node(c);
        let mut dmax = psize;
        let mut dsq: f32 = 0.0;
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

    fn walksub(&self, ctx: &mut WalkContext) -> Result<()> {
        let poff = ctx.psize / 4.0;
        let pnode = *self.node(ctx.p);
        if let NodeRef::Cell(pid) = ctx.p {
            let pnext = pnode.next;
            let mut q = self.cells[pid].more;
            let mut children = Vec::with_capacity(NSUB);
            while q != pnext {
                if q.is_none() {
                    return Err(TreeError::TreeStructure);
                }
                children.push(q);
                q = self.node(q).next;
            }
            if ctx.parallel_depth > 0 {
                self.fan_out_parallel(&children, ctx)?;
            } else {
                self.walk_sequential(&children, ctx)?;
            }
        } else {
            let nmid = next_midpoint(ctx.pmid, pnode.pos, poff);
            ctx.aptr = ctx.nptr;
            ctx.nptr = ctx.np;
            ctx.psize /= 2.0;
            ctx.pmid = nmid;
            ctx.parallel_depth = 0;
            self.walktree(ctx)?;
        }
        Ok(())
    }

    #[inline]
    fn walk_sequential(&self, children: &[NodeRef], ctx: &mut WalkContext) -> Result<()> {
        let poff = ctx.psize / 4.0;
        for qr in children {
            let nmid = next_midpoint(ctx.pmid, self.node(*qr).pos, poff);
            let mut child = WalkContext {
                active: &mut *ctx.active,
                interact: &mut *ctx.interact,
                aptr: ctx.nptr,
                nptr: ctx.np,
                cptr: ctx.cptr,
                bptr: ctx.bptr,
                np: ctx.np,
                p: *qr,
                psize: ctx.psize / 2.0,
                pmid: nmid,
                results: ctx.results,
                counters: ctx.counters,
                parallel_depth: 0,
            };
            self.walktree(&mut child)?;
        }
        Ok(())
    }

    #[inline]
    fn fan_out_parallel(&self, children: &[NodeRef], ctx: &WalkContext) -> Result<()> {
        self.walk_parallel(children, ctx)
    }

    #[inline]
    fn walk_parallel(&self, children: &[NodeRef], ctx: &WalkContext) -> Result<()> {
        let err_flag = Mutex::new(None::<TreeError>);
        rayon::scope(|s| {
            for qr in children {
                s.spawn(|_| self.run_child_walk(*qr, ctx, &err_flag));
            }
        });
        if let Some(e) = err_flag.into_inner().unwrap_or_default() {
            return Err(e);
        }
        Ok(())
    }

    #[inline]
    fn run_child_walk(&self, qr: NodeRef, ctx: &WalkContext, err_flag: &Mutex<Option<TreeError>>) {
        let actlen = self.actlen;
        let poff = ctx.psize / 4.0;
        let nmid = next_midpoint(ctx.pmid, self.node(qr).pos, poff);
        let mut active2 = Vec::new();
        if active2.try_reserve(actlen).is_err() {
            set_result(err_flag, Err(TreeError::OutOfMemory(actlen)));
            return;
        }
        active2.resize(actlen, NodeRef::Body(0));
        active2[0..(ctx.np - ctx.nptr)].copy_from_slice(&ctx.active[ctx.nptr..ctx.np]);
        let mut interact2 = Vec::new();
        if interact2.try_reserve(actlen).is_err() {
            set_result(err_flag, Err(TreeError::OutOfMemory(actlen)));
            return;
        }
        interact2.resize(actlen, Interact::default());
        interact2[0..ctx.cptr].copy_from_slice(&ctx.interact[0..ctx.cptr]);
        interact2[ctx.bptr..actlen].copy_from_slice(&ctx.interact[ctx.bptr..actlen]);
        let mut child = WalkContext {
            active: &mut active2,
            interact: &mut interact2,
            aptr: 0,
            nptr: ctx.np - ctx.nptr,
            cptr: ctx.cptr,
            bptr: ctx.bptr,
            np: ctx.np - ctx.nptr,
            p: qr,
            psize: ctx.psize / 2.0,
            pmid: nmid,
            results: ctx.results,
            counters: ctx.counters,
            parallel_depth: ctx.parallel_depth.saturating_sub(1),
        };
        set_result(err_flag, self.walktree(&mut child));
    }

    fn gravsum(
        &self,
        interact: &[Interact],
        p0: NodeRef,
        cptr: usize,
        bptr: usize,
        results: &Mutex<Vec<(f32, Vector)>>,
        counters: &Mutex<WalkCounters>,
    ) {
        let pos0 = match p0 {
            NodeRef::Body(b) => self.bodytab[b].bodynode.pos,
            _ => unreachable!(),
        };
        let eps2 = self.eps2;
        let usequad = self.usequad != 0;
        let actlen = self.actlen;
        let mut phi0: f32 = 0.0;
        let mut acc0: Vector = Vector::zero();
        #[cfg(feature = "simd")]
        {
            if usequad {
                sumcell_simd(eps2, interact, 0, cptr, pos0, &mut phi0, &mut acc0);
            } else {
                sumnode_simd(eps2, interact, 0, cptr, pos0, &mut phi0, &mut acc0);
            }
            sumnode_simd(eps2, interact, bptr, actlen, pos0, &mut phi0, &mut acc0);
        }
        #[cfg(not(feature = "simd"))]
        {
            if usequad {
                sumcell(eps2, interact, 0, cptr, pos0, &mut phi0, &mut acc0);
            } else {
                sumnode(eps2, interact, 0, cptr, pos0, &mut phi0, &mut acc0);
            }
            sumnode(eps2, interact, bptr, actlen, pos0, &mut phi0, &mut acc0);
        }
        if let NodeRef::Body(b) = p0 {
            if let Ok(mut g) = results.lock() {
                g[b] = (phi0, acc0);
            }
            if let Ok(mut c) = counters.lock() {
                c.nbbcalc += actlen - bptr;
                c.nbccalc += cptr;
            }
        }
    }

    pub fn force_max_active(&self) -> usize {
        self.actmax
    }

    pub fn force_bb_calc(&self) -> usize {
        self.nbbcalc
    }

    pub fn force_bc_calc(&self) -> usize {
        self.nbccalc
    }

    pub fn force_cpu_time(&self) -> f64 {
        self.cpuforce as f64
    }
}
