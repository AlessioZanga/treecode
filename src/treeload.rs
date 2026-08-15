use crate::{
    error::{Result, TreeError},
    mathfns,
    treecode::{MAXLEVEL, Tree},
    types::{
        BodyId, CELL, Cell, CellId, NDIM, NSUB, Node, NodeRef, Sorq, Vector, cputime, scanopt,
    },
    vecmath::{Matrix, matrix_zero, vector_zero},
};

fn set_center_of_mass(p: &Cell, cmpos: &mut Vector) {
    if p.cellnode.mass > 0.0 {
        for k in 0..NDIM {
            cmpos[k] /= p.cellnode.mass;
        }
    } else {
        for k in 0..NDIM {
            cmpos[k] = p.cellnode.pos[k];
        }
    }
}

fn verify_center(p: &Cell, cmpos: &Vector, psize: f32) -> Result<()> {
    for k in 0..NDIM {
        if cmpos[k] < p.cellnode.pos[k] - psize / 2.0 || p.cellnode.pos[k] + psize / 2.0 <= cmpos[k]
        {
            return Err(TreeError::TreeStructure);
        }
    }
    Ok(())
}

fn compute_rcrit_exact(rsize: f32) -> f32 {
    mathfns::rsqr(2.0 * rsize)
}

fn compute_rcrit_sw94(cmpos: &Vector, psize: f32, theta2: f32, p: &Cell) -> f32 {
    let mut bmax2: f32 = 0.0;
    for k in 0..NDIM {
        let d = cmpos[k] - p.cellnode.pos[k] + psize / 2.0;
        bmax2 += mathfns::rsqr(d.max(psize - d));
    }
    bmax2 / theta2
}

fn compute_rcrit_bh86(psize: f32, theta: f32) -> f32 {
    mathfns::rsqr(psize / theta)
}

fn compute_rcrit_default(cmpos: &Vector, psize: f32, theta: f32, p: &Cell) -> f32 {
    let mut d: f32 = 0.0;
    for k in 0..NDIM {
        let dk = cmpos[k] - p.cellnode.pos[k];
        d += dk * dk;
    }
    mathfns::rsqr(psize / theta + d.sqrt())
}

fn dot_product(dr: &Vector) -> f32 {
    let mut drsq: f32 = 0.0;
    for k in 0..NDIM {
        drsq += dr[k] * dr[k];
    }
    drsq
}

impl Tree {
    pub(crate) fn node(&self, r: NodeRef) -> &Node {
        match r {
            NodeRef::Body(b) => &self.bodytab[b].bodynode,
            NodeRef::Cell(c) => &self.cells[c].cellnode,
            NodeRef::None => unreachable!("node() called on NodeRef::None"),
        }
    }

    pub(crate) fn node_mut(&mut self, r: NodeRef) -> &mut Node {
        match r {
            NodeRef::Body(b) => &mut self.bodytab[b].bodynode,
            NodeRef::Cell(c) => &mut self.cells[c].cellnode,
            NodeRef::None => unreachable!("node_mut() called on NodeRef::None"),
        }
    }

    pub fn maketree(&mut self, nbody: usize) -> Result<()> {
        let cpustart = cputime()?;
        self.newtree();
        let root = self.makecell()?;
        self.root = Some(root);
        vector_zero(&mut self.cells[root].cellnode.pos);
        self.expandbox()?;
        self.load_all_bodies(nbody)?;
        self.parse_options()?;
        self.reset_stats();
        self.hackcofm(root, self.rsize, 0)?;
        self.threadtree(NodeRef::Cell(root), NodeRef::None)?;
        if self.usequad != 0 {
            self.hackquad(root)?;
        }
        self.cputree = (cputime()? - cpustart) as f32;
        Ok(())
    }

    fn load_all_bodies(&mut self, nbody: usize) -> Result<()> {
        for i in 0..nbody {
            self.loadbody(i)?;
        }
        Ok(())
    }

    fn parse_options(&mut self) -> Result<()> {
        self.bh86 = scanopt(&self.options, "bh86");
        self.sw94 = scanopt(&self.options, "sw94");
        if self.bh86 && self.sw94 {
            return Err(TreeError::IncompatibleOptions);
        }
        Ok(())
    }

    fn reset_stats(&mut self) {
        self.tdepth = 0;
        for i in 0..MAXLEVEL {
            self.cellhist[i] = 0;
            self.subnhist[i] = 0;
        }
    }

    fn newtree(&mut self) {
        if !self.firstcall {
            let mut p: NodeRef = self.root.map(NodeRef::Cell).unwrap_or(NodeRef::None);
            while !p.is_none() {
                if self.node(p).node_type == CELL {
                    let cid = p.index();
                    self.freecell.push(cid);
                    p = self.cells[cid].more;
                } else {
                    p = self.node(p).next;
                }
            }
        } else {
            self.firstcall = false;
        }
        self.root = None;
        self.ncell = 0;
    }

    fn makecell(&mut self) -> Result<CellId> {
        let id = if let Some(free_id) = self.freecell.pop() {
            free_id
        } else {
            self.cells.push(Cell::default());
            self.cells.len() - 1
        };
        self.cells[id].cellnode.node_type = CELL;
        self.cells[id].cellnode.update = 0;
        self.cells[id].more = NodeRef::None;
        self.cells[id].sorq = Sorq::default();
        self.ncell += 1;
        Ok(id)
    }

    fn expandbox(&mut self) -> Result<()> {
        let root = self.root.ok_or(TreeError::TreeStructure)?;
        let mut dmax: f32 = 0.0;
        for p in &self.bodytab {
            for k in 0..NDIM {
                let d = (p.bodynode.pos[k] - self.cells[root].cellnode.pos[k]).abs();
                if d > dmax {
                    dmax = d;
                }
            }
        }
        while self.rsize < 2.0 * dmax {
            self.rsize *= 2.0;
        }
        Ok(())
    }

    fn loadbody(&mut self, p: BodyId) -> Result<()> {
        let mut q: CellId = self.root.ok_or(TreeError::TreeStructure)?;
        let mut qind = self.subindex(p, q);
        let mut qsize = self.rsize;
        loop {
            let cur = self.cells[q].sorq.subp()[qind];
            if cur.is_none() {
                break;
            }
            if !cur.is_cell() {
                let other = cur.index();
                self.require_distinct_positions(p, other)?;
                let c = self.makecell()?;
                self.set_cell_midpoint(c, q, p, qsize);
                let sub = self.subindex(other, c);
                self.cells[c].sorq.subp_mut()[sub] = NodeRef::body(other);
                self.cells[q].sorq.subp_mut()[qind] = NodeRef::cell(c);
                q = c;
            } else {
                q = cur.index();
            }
            qind = self.subindex(p, q);
            qsize /= 2.0;
        }
        self.cells[q].sorq.subp_mut()[qind] = NodeRef::body(p);
        Ok(())
    }

    fn require_distinct_positions(&self, p: BodyId, other: BodyId) -> Result<()> {
        let mut dist2: f32 = 0.0;
        for k in 0..NDIM {
            let d = self.bodytab[p].bodynode.pos[k] - self.bodytab[other].bodynode.pos[k];
            dist2 += d * d;
        }
        if dist2 == 0.0 {
            return Err(TreeError::CoincidentBodies);
        }
        Ok(())
    }

    fn subindex(&self, p: BodyId, q: CellId) -> usize {
        let mut ind: usize = 0;
        for k in 0..NDIM {
            ind |= usize::from(self.cells[q].cellnode.pos[k] <= self.bodytab[p].bodynode.pos[k])
                << (NDIM - 1 - k);
        }
        ind
    }

    fn set_cell_midpoint(&mut self, c: CellId, q: CellId, p: BodyId, qsize: f32) {
        for k in 0..NDIM {
            let pk = self.bodytab[p].bodynode.pos[k];
            let qk = self.cells[q].cellnode.pos[k];
            let offset = (if pk < qk { -qsize } else { qsize }) / 4.0;
            self.cells[c].cellnode.pos[k] = qk + offset;
        }
    }

    fn hackcofm(&mut self, p: CellId, psize: f32, lev: usize) -> Result<()> {
        let mut cmpos: Vector = Vector::zero();
        self.update_depth(lev);
        self.cellhist[lev] += 1;
        self.cells[p].cellnode.mass = 0.0;
        self.accumulate_subnodes(p, psize, lev, &mut cmpos)?;
        set_center_of_mass(&self.cells[p], &mut cmpos);
        verify_center(&self.cells[p], &cmpos, psize)?;
        self.setrcrit(p, &cmpos, psize);
        for k in 0..NDIM {
            self.cells[p].cellnode.pos[k] = cmpos[k];
        }
        Ok(())
    }

    fn update_depth(&mut self, lev: usize) {
        if lev > self.tdepth {
            self.tdepth = lev;
        }
    }

    fn accumulate_subnodes(
        &mut self,
        p: CellId,
        psize: f32,
        lev: usize,
        cmpos: &mut Vector,
    ) -> Result<()> {
        let mut tmpv: Vector = Vector::zero();
        let subp = *self.cells[p].sorq.subp();
        for sub in &subp {
            if sub.is_none() {
                continue;
            }
            self.subnhist[lev] += 1;
            if sub.is_cell() {
                let cid = sub.index();
                self.hackcofm(cid, psize / 2.0, lev + 1)?;
            }
            let qnode = *self.node(*sub);
            self.cells[p].cellnode.update |= qnode.update;
            self.cells[p].cellnode.mass += qnode.mass;
            for k in 0..NDIM {
                tmpv[k] = qnode.pos[k] * qnode.mass;
                cmpos[k] += tmpv[k];
            }
        }
        Ok(())
    }

    fn setrcrit(&mut self, p: CellId, cmpos: &Vector, psize: f32) {
        let rcrit2 = if self.theta == 0.0 {
            compute_rcrit_exact(self.rsize)
        } else if self.sw94 {
            compute_rcrit_sw94(cmpos, psize, self.theta2, &self.cells[p])
        } else if self.bh86 {
            compute_rcrit_bh86(psize, self.theta)
        } else {
            compute_rcrit_default(cmpos, psize, self.theta, &self.cells[p])
        };
        self.cells[p].rcrit2 = rcrit2;
    }

    fn displacement(&self, q: NodeRef, p: CellId) -> Vector {
        let mut dr: Vector = Vector::zero();
        for k in 0..NDIM {
            dr[k] = self.node(q).pos[k] - self.cells[p].cellnode.pos[k];
        }
        dr
    }

    fn quadrupole_tensor(&self, q: NodeRef, dr: &Vector) -> Matrix {
        let drsq = dot_product(dr);
        let mut tmpm: Matrix = Matrix::zero();
        for j in 0..NDIM {
            for k in 0..NDIM {
                let id_rsq = if j == k { drsq } else { 0.0 };
                tmpm[j][k] = (3.0 * dr[j] * dr[k] - id_rsq) * self.node(q).mass;
            }
        }
        tmpm
    }

    fn collect_descendants(&self, c: CellId, desc: &mut [NodeRef]) -> usize {
        let mut ndesc: usize = 0;
        let subp = self.cells[c].sorq.subp();
        for s in subp.iter() {
            if s.is_none() {
                continue;
            }
            desc[ndesc] = *s;
            ndesc += 1;
        }
        ndesc
    }

    fn threadtree(&mut self, p: NodeRef, n: NodeRef) -> Result<()> {
        self.node_mut(p).next = n;
        if self.node(p).node_type == CELL {
            let cid = p.index();
            let mut desc: [NodeRef; NSUB + 1] = [NodeRef::None; NSUB + 1];
            let ndesc = self.collect_descendants(cid, &mut desc);
            self.cells[cid].more = desc[0];
            desc[ndesc] = n;
            for i in 0..ndesc {
                let child = desc[i];
                if child.is_none() {
                    return Err(TreeError::TreeStructure);
                }
                self.threadtree(child, desc[i + 1])?;
            }
        }
        Ok(())
    }

    fn hackquad(&mut self, p: CellId) -> Result<()> {
        let mut desc: [NodeRef; NSUB] = [NodeRef::None; NSUB];
        let ndesc = self.collect_descendants(p, &mut desc);
        matrix_zero(self.cells[p].sorq.quad_mut());
        for d in &desc[..ndesc] {
            let q = *d;
            if q.is_none() {
                return Err(TreeError::TreeStructure);
            }
            if q.is_cell() {
                self.hackquad(q.index())?;
            }
            self.accumulate_moment(p, q);
        }
        Ok(())
    }

    fn accumulate_moment(&mut self, p: CellId, q: NodeRef) {
        let dr = self.displacement(q, p);
        let mut tmpm = self.quadrupole_tensor(q, &dr);
        if let NodeRef::Cell(cid) = q {
            let qm = self.cells[cid].sorq.quad();
            for j in 0..NDIM {
                for k in 0..NDIM {
                    tmpm[j][k] += qm[j][k];
                }
            }
        }
        for j in 0..NDIM {
            for k in 0..NDIM {
                self.cells[p].sorq.quad_mut()[j][k] += tmpm[j][k];
            }
        }
    }

    pub fn tree_depth(&self) -> usize {
        self.tdepth
    }

    pub fn cell_count(&self) -> usize {
        self.ncell
    }

    pub fn tree_build_time(&self) -> f64 {
        self.cputree as f64
    }
}
