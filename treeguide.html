<HTML>

<HEAD>
<TITLE>Treecode Guide</TITLE>
</HEAD>

<BODY BGCOLOR=WHITE>

<H1>TREECODE GUIDE</H1>

<H2>Joshua E. Barnes</H2>

<H4><I>Institute for Astronomy, University of Hawai`i,<BR>
2680 Woodlawn Drive, Honolulu, Hawai`i 96822</I></H4>

<HR>

<P><B>treecode</B> is a program for self-consistent N-body simulation.
It is faster than previous codes, and generally provides better error
control.  This document describes <B>treecode</B> Version 1.4; it
provides an overview of the algorithm, instructions for compiling and
running the code, and links to the actual sources.</P>

<H2>0. Introduction</H2>

<P>Hierarchical force calculation algorithms (<I>e.g.</I> Greengard
1990) provide fast, general, and reasonably accurate approximations
for gravity and other inverse-square forces.  They fill the gap
between direct sum methods, which are accurate and general but require
O(<I>N</I><TT><SUP>2</SUP></TT>) operations for a complete force
calculation, and field methods, which have limited generality and
accuracy but require only O(<I>N</I>) operations.  All hierarchical
methods partition the mass distribution into a tree structure, where
each node of the tree provides a concise description of the matter
within some spatial volume.  This tree structure is used to introduce
explicit approximations into the force calculation.  Hierarchical
methods require either O(<I>N</I>) or O(<I>N</I> log <I>N</I>)
operations per force calculation, depending on the representation
employed.  The algorithm described here improves on an earlier
hierarchical O(<I>N</I> log <I>N</I>) method (Barnes &amp; Hut 1986,
hereafter BH86) which has been widely employed in astrophysical
simulations.</P>

<P>Below, <A HREF="#strategies">Section 1</A> outlines the general
strategy of the new code.  The data structures used in the tree
construction and force calculation are discussed in <A
HREF="#data">Section 2</A>.  Routines for tree construction, force
calculation, control &amp; integration, and input/output are described
in <A HREF="#routines">Section 3</A>.  Instructions for copying,
compiling, and running <B>treecode</B> are given in <A
HREF="#instructions">Section 4</A>.  Differences between Version 1.4
and earlier versions are discussed in <A HREF="#appendix_a">Appendix
A</A>.  Finally, system-level software used in the code is described
in the <A HREF="#appendix_b">Appendix B</A>.</P>

<A NAME="strategies"></A>
<H2>1. Strategies</H2>

<P><B>treecode</B> reduces the overhead of tree search by using the
fact that neighboring bodies have similar interaction lists.  This
idea was previously used to speed up hierarchical force calculation on
vector machines (Barnes 1990), but the earlier code simply performed a
tree search for some small volume defined by a single cell.  In the
new code this idea is applied to all levels of the tree.</P>

<P>Forces on bodies are computed during a single recursive scan of the
entire tree.  This scan maintains and continuously updates an
<I>interaction list</I>; at each level, the fact that a body
<TT>b</TT> lies somewhere within a cell <TT>c</TT> is used to winnow
the set of possible interactions that <TT>b</TT> might have.  (The
information used to define interaction lists is thus similar to that
used in an early parallel code (Barnes 1986), and a similar strategy
figures in the Fast Multipole Method (Greengard &amp; Rokhlin 1987).)
The cost of this winnowing is spread over all bodies within
<TT>c</TT>, promising a significant speed-up over codes which
construct a separate interaction list for each body (<I>e.g.</I>
BH86).  When the recursive scan arrives at a body <TT>b</TT>, the
interaction list on hand is used to obtain the gravitational force and
potential at <TT>b</TT>.</P>

<P>Typical hierarchical N-body codes spend about half their time
searching the tree (Hernquist 1987); thus one might hope that a factor
of two could be gained by quickly generating interaction lists for all
bodies in a single scan of the tree.  In practice, since the
interaction lists so generated are not precisely tailor-made for
individual bodies, they must be somewhat longer for a given accuracy.
Thus to obtain the expected speed-up requires faster evaluation of
interactions, achieved largely via tighter coding.  Preliminary timing
tests show that <B>treecode</B> is indeed about a factor of two faster
than the BH86 code <I>at a fixed level of accuracy</I>; it also
outperforms other versions of that algorithm (Barnes 1998).  Moreover,
it is significantly more robust when presented with `perverse' mass
distributions (Salmon &amp; Warren 1994).</P>

<P>The new algorithm should lend itself to efficient use of vector
processors.  On a scalar machine, it spends 90% to 95% of the time
summing up interaction lists; a large speed-up would be realized by
offloading this task to a vector unit.  This approach can make
efficient use of N-body force calculation hardware such as the GRAPE
(Gravity PIpe) processor (Sugimoto <I>et al.</I> 1990).</P>

<P>Parallel implementation should also be quite straight-forward
provided that each processor has enough memory to hold a copy of the
entire particle array and tree structure.  At a minimum, the task of
summing each interaction list can be distributed across a small number
of processors; the bottleneck is then that each processor must search
the entire tree.  A better approach is to order bodies as they would
be encountered in a tree-walk, estimate the computational work
required to calculate the force on each, and give each processor the
job of computing forces for a contiguous block of bodies.</P>

<A NAME="data"></A>
<H2>2. Data Structures</H2>

<P>The file <A HREF="treedefs.h.html"><TT>treedefs.h</TT></A> has
specifications for the data structures and global variables used in
tree construction and force calculation.  This file requires
definitions from several general-purpose C include files, which are
further described in the <A HREF="#appendix_b">Appendix B</A>.</P>

<H3>2.1. Tree Structure</H3>

<P>The primary data structure used in the code is an eight-way tree
(or `oct-tree'), composed of <I>bodies</I> and <I>cells</I>.  Bodies
are the leaves of the tree, while cells are the internal branches.
The entire tree can be reached starting from a single cell, known as
the <I>root</I>.  A simple tree is shown here:</P>

<IMG BORDER=1 ALIGN=MIDDLE SRC="fig1.gif">

<P><A HREF="treedefs.h.html#node"><B>node</B></A> structures contain
the information common to both bodies and cells.  In principle, each
element of the tree could be represented as a <TT>union</TT> of a body
and a cell, but this would be inefficient since bodies and cells
require different amounts of memory.  Instead, a node is a common
header for bodies and cells, and casts are used to convert a pointer
of arbitrary type into a pointer to a node, cell or a body.  All of
this ugly code is hidden in macros; the macros used to access the
components of nodes are:</P>

<UL>

<LI><B>Type(q)</B> is the actual type of node <TT>q</TT>; the only
possible values are the constants <TT>BODY</TT> and <TT>CELL</TT>.

<LI><B>Update(q)</B> is a boolean value which governs the scope of
force calculation.  Bodies with <TT>Update</TT> equal to <TT>TRUE</TT>
get updated forces; cells with <TT>Update</TT> equal to <TT>TRUE</TT>
are searched for bodies needing updated forces.

<LI><B>Next(q)</B> is a pointer to the next node encountered in a
recursive scan of the tree <I>after</I> <TT>q</TT>'s descendents (if
any) have been visited.

<LI><B>Mass(q)</B> is the mass of a body, or the total mass of all
bodies within a cell.

<LI><B>Pos(q)</B> is the position of a body, or the
position of the center of mass of all bodies within a cell.

</UL>

<P><A HREF="treedefs.h.html#body"><B>body</B></A> structures represent
particles.  The macros used to access the components of bodies,
besides those for nodes, are:</P>

<UL>

<LI><B>Vel(b)</B> is the velocity of body <TT>b</TT>.

<LI><B>Acc(b)</B> is the  acceleration of body <TT>b</TT>.

<LI><B>Phi(b)</B> is the potential of body <TT>b</TT>.

</UL>

<P><A HREF="treedefs.h.html#cell"><B>cell</B></A> structures represent
the eight-way internal branchings of the tree.  The macros used to
access the components of cells, besides those for nodes, are:</P>

<UL>

<LI><B>Subp(c)</B> is an array of pointers to the descendents of cell
<TT>c</TT>.

<LI><B>More(c)</B> is a pointer to the first of these descendents.

<LI><B>Quad(c)</B> is a matrix of quadrupole moments.

<LI><B>Rcrit2(c)</B> is the square of the critical radius beyond which
the cell <TT>c</TT> can safely appear in interaction lists.  This is
<I>not</I> defined if the <TT>QUICKSCAN</TT> version is compiled.

</UL>

<P>(Note that the <TT>Subp</TT> and <TT>Quad</TT> fields share memory;
thus only one of each is defined at any point in the calculation.)</P>

<P>The <TT>Next</TT> and <TT>More</TT> fields of the tree are
initialized in the second phase of tree construction.  In this phase,
the structure is `threaded' in such a way that a tree search can be
performed by a simple iterative procedure.  This transforms the tree
sketched above into the following:</P>

<IMG BORDER=1 ALIGN=MIDDLE SRC="fig2.gif">

<P>In essence, the <TT>More</TT> link is followed to get more detailed
information on the mass distribution, and the <TT>Next</TT> link is
followed to move on to the next part of the tree.  Threading was
originally done to speed up force calculation; now, however, it's used
to free up the memory previously used to store the <TT>Subp</TT>
array.</P>

<P>Before threading, the standard idiom to iterate over the
immediate descendents of a <TT>cellptr p</TT> is
<PRE>
    int i;
    for (i = 0; i &lt; NSUB; i++)
        if (Subp(p)[i] != NULL)
            &lt;process Subp(p)[i]&gt;;
</PRE>
where <TT>NSUB = 2<SUP>3</SUP></TT> is the maximum number of direct
descendents, while after threading, the idiom is
<PRE>
    nodeptr q;
    for (q = More(p); q != Next(p); q = Next(q))
        &lt;process q&gt;;
</PRE>
</P>

<H3>2.2. Global Variables</H3>

<P>Global variables are declared with the symbol <A
HREF="treedefs.h.html#global"><B>global</B></A>, introduced to deal
with ANSI C's requirement that the <TT>extern</TT> keyword be used in
all but one compilation module.  In compiling the tree construction,
force calculation, and input/output routines, <TT>global</TT> expands
to <TT>extern</TT>; in the main module <A
HREF="treecode.c.html"><TT>treecode.c</TT></A>, the <TT>global</TT>
symbol is predefined before <A
HREF="treedefs.h.html"><TT>treedefs.h</TT></A> is included:
<PRE>
    #define global  /* don't default to extern  */
</PRE>
This definition prevents <TT>global</TT> from expanding to
<TT>extern</TT> when compiling the main module.</P>

<P><B>a) <A HREF="treedefs.h.html#parameters">Input parameters</A></B>
for tree construction and force calculation are as follows:</P>

<UL>

<LI><B>theta</B> governs the accuracy of the force calculation.  This
parameter is <I>not</I> defined if the <TT>QUICKSCAN</TT> version is
compiled.

<LI><B>options</B> is a string listing various run-time control
options.

<LI><B>usequad</B> is a flag governing the use of quadrupole
corrections.

</UL>

<P><B>b) <A HREF="treedefs.h.html#tree">Tree construction</A></B>
assigns values to the following variables:</P>

<UL>

<LI><B>root</B> is a pointer to the root cell.

<LI><B>rsize</B> is the linear size of the root cell.

<LI><B>ncell</B> counts the number of cells used to build the tree.

<LI><B>tdepth</B> counts the number of levels in the tree below the
root.

<LI><B>cputree</B> is the CPU time required for tree construction.

</UL>

<P><B>c) <A HREF="treedefs.h.html#force">Force calculation</A></B>
assigns values to the following variables:</P>

<UL>

<LI><B>actmax</B> is the maximum length of the active list during
force calculation.

<LI><B>nbbcalc</B> is the total number of interactions between bodies
and bodies.

<LI><B>nbccalc</B> is the total number of interactions between bodies
and cells.

<LI><B>cpuforce</B> is the CPU time required for force calculation.

</UL>

<A NAME="routines"></A>
<H2>3. Routines</H2>

<P>The routines implementing <B>treecode</B> are grouped into four
categories: tree construction, force calculation, control &amp
integration, and input/output.  Each set of routines is defined in a
separate file.</P>

<H3>3.1. Tree Construction</H3>

<P>The tree construction task is handled by the routines in <A
HREF="treeload.c.html"><TT>treeload.c</TT></A>.  Construction begins
by fitting a cube, of linear dimension <TT>rsize</TT>, around the body
positions.  This cube is identified with the <TT>root</TT> cell of the
tree.  Bodies are then loaded into the tree one at a time.  A given
cell can hold up to eight bodies if each one happens to fall in a
different octant.  Whenever two bodies fall in the <I>same</I> octant,
a new cell is allocated and the tree is extended to another level.
Once all bodies are loaded, total mass and center-of-mass information
is propagated from the bodies towards the root.</P>

<P><A HREF="treeload.c.html#maketree"><B>maketree</B></A> supervises the
process of tree construction.  Its prototype is
<PRE>
    void maketree(bodyptr btab, int nbody);
</PRE>
where <TT>btab</TT> is an array of <TT>nbody</TT> body structures.
When <TT>maketree</TT> returns, the root of the tree is addressed by
the global cell pointer, <TT>root</TT>; also updated are the root cell
size, <TT>rsize</TT>, the number of cells in the tree, <TT>ncell</TT>,
the depth of the tree, <TT>tdepth</TT> (counting from <TT>0</TT>), and
the CPU time used, <TT>cputree</TT>.</P>

<P>After initializing an empty <TT>root</TT> cell, <TT>maketree</TT>
loops over <TT>btab</TT>, calling <TT>loadbody</TT> to install each
body in the tree.  It then propagates information from the leaves to
the root, threads the tree structure, and optionally computes
quadrupole moments.</P>

<P><A HREF="treeload.c.html#newtree"><B>newtree</B></A> scavenges the
cells in the existing tree and prepares to build a new one.  Its
prototype is
<PRE>
    void newtree(void);
</PRE>
</P>

<P><A HREF="treeload.c.html#makecell"><B>makecell</B></A> returns a
pointer to a free cell.  Its prototype is</P>
<PRE>
    cellptr makecell(void);
</PRE>
In addition, <TT>makecell</TT> also updates <TT>ncell</TT>.</P>

<P><A HREF="treeload.c.html#expandbox"><B>expandbox</B></A> fits a box
around the body distribution.  Its prototype is</P>
<PRE>
    void expandbox(bodyptr btab, int nbody);
</PRE>
where <TT>btab</TT> is an array of <TT>nbody</TT> body structures.
The result is stored in <TT>rsize</TT>.</P>

<P>To take advantage of the exact floating-point representation of
powers of two in binary computers, the box size is successively
doubled until it fits the bodies.  This insures that the corners and
midpoints of cells have exact binary representations.</P>

<P><A HREF="treeload.c.html#loadbody"><B>loadbody</B></A> inserts a body
in the tree.  Its prototype is</P>
<PRE>
    void loadbody(bodyptr p);
</PRE>
where <TT>p</TT> is the body in question.</P>

<P>To find the appropriate place to insert the body, <TT>loadbody</TT>
traces out a path from the <TT>root</TT> node toward the leaves of the
tree.  At each level, it calls <TT>subindex</TT> to decide which
branch to take.  Eventually, one of two things happens.  First, the
indicated subcell may be empty; the body <TT>p</TT> is then stuffed
into the empty slot.  Second, the indicated subcell may already hold a
body; then a new cell is allocated to extend the tree and both bodies
are installed within the new cell.</P>

<P><A HREF="treeload.c.html#subindex"><B>subindex</B></A> decides which
of a cell's octants a body falls in.  Its prototype is</P>
<PRE>
    int subindex(bodyptr p, cellptr q);
</PRE>
where <TT>p</TT> is the body in question and <TT>c</TT> is the cell.

<P>The <TT>subindex</TT> function uses the fact that during the first
phase of tree construction the <TT>Pos</TT> vector of a cell is its
geometric midpoint.  Note that it's assumed that <TT>p</TT> actually
lies within the volume represented by <TT>c</TT>; this assumption is
checked for the entire tree when <TT>hackcofm</TT> is called.</P>

<P><A HREF="treeload.c.html#hackcofm"><B>hackcofm</B></A> propagates
cumulative information towards the <TT>root</TT> cell.  Its
prototype is</P>
<PRE>
    void hackcofm(cellptr p, real psize, int lev);
</PRE>
where <TT>p</TT> is the current cell, <TT>psize</TT> is its linear
size, and <TT>lev</TT> is its level in the tree, counting the
<TT>root</TT> as level <TT>0</TT>.</P>

<P>In outline, <TT>hackcofm</TT> is a simple depth-first recursive
tree scan; it loops over the descendents of <TT>p</TT> and calls
itself on those which point to cells. It thereby accumulates total
masses, center-of-mass positions, and update flags.  This done,
<TT>hackcofm</TT> checks that the center-of-mass actually lies in the
volume represented by <TT>p</TT>; failure indicates a bug in tree
construction.  Optionally, <TT>setrcrit</TT> is called to set the
critical opening radius for this cell.</P>

<P><A HREF="treeload.c.html#setrcrit"><B>setrcrit</B></A> assigns each
cell a critical squared opening radius.  Its prototype is</P>
<PRE>
    void setrcrit(cellptr p, vector cmpos, real psize);
</PRE>
where <TT>p</TT> is the cell, <TT>cmpos</TT> is its center-of-mass
position, and <TT>psize</TT> is its size.</P>

<P>This routine is defined only if the <TT>QUICKSCAN</TT> version is
<I>not</I> compiled.  It offers several different methods to compute
the squared critical radius, <TT>Rcrit2</TT>, beyond which a force
calculation need not open <TT>p</TT>.  The default is a criterion
which is extra careful with cells having relatively off-center mass
distributions (Barnes 1995); also provided are Salmon &amp; Warren's
(1994) <I>bmax</I> criterion (option <TT>sw94</TT>) and BH86's
original criterion (option <TT>bh86</TT>).</P>

<P><A HREF="treeload.c.html#threadtree"><B>threadtree</B></A> relinks
the tree using the <TT>Next</TT> and <TT>More</TT> pointers.  Its
prototype is</P>
<PRE>
    void threadtree(nodeptr p, nodeptr n);
</PRE>
where <TT>p</TT> is any node and <TT>n</TT> will be its successor.</P>

<P>In outline, <TT>threadtree</TT> is also a depth-first recursive
tree walk.  First, it makes <TT>n</TT> the <TT>Next</TT> link of
<TT>p</TT>.  Then, if <TT>p</TT> is a cell, it makes a list its
descendents and sets <TT>p</TT>'s <TT>More</TT> link to be the first
of these; finally it calls itself on each member of the list, passing
the next member as the successor (note that the successor of the last
member is the successor of <TT>p</TT>).</P>

<P><A HREF="treeload.c.html#hackquad"><B>hackquad</B></A> propagates
quadrupole moments towards the <TT>root</TT> cell.  Its prototype
is</P>
<PRE>
    void hackquad(cellptr p);
</PRE>
where <TT>p</TT> is the current cell.</P>

<P>Again, this routine does a depth-first tree walk.  The only trick
is that the storage used for the <TT>Subp</TT> pointers is reused to
store the quadrupole matrix, so the latter are first copied to a
temporary array.  (Note: this routine could be simplified using the
post-threading idiom to iterate over descendents.)</P>

<H3>3.2. Force Calculation</H3>

<P>The force calculation code is implemented by routines provided in
<A HREF="treegrav.c.html"><TT>treegrav.c</TT></A>.  As described
above, forces are calculated during a single recursive scan of the
tree, which visits every body whose <TT>Update</TT> flag is set.
Gravitational forces and potentials are assigned to these bodies.</P>

<P><A HREF="treegrav.c.html#gravcalc"><B>gravcalc</B></A> supervises
the process of force calculation.  Its prototype is
<PRE>
    void gravcalc(void);
</PRE>
The tree structure to be used by <TT>gravcalc</TT> is addressed by
the global <TT>root</TT> pointer; also referenced are the tree depth
<TT>tdepth</TT> and root cell size <TT>rsize</TT>.</P>

<P>The main job of <TT>gravcalc</TT> is to set up the initial call to
the worker routine <TT>walktree</TT>.  It begins by allocating
temporary storage for several lists; the length of these lists is
estimated from the depth of the tree.  The <TT>interact</TT> pointer
addresses a linear array of cells which list all the interactions
acting on a body; body-cell interactions are listed from the front
toward the back, while body-body interactions are listed from the back
toward the front.  The <TT>active</TT> pointer addresses an array of
node pointers which will be examined when constructing interaction
lists.  With these arrays in place, <TT>gravcalc</TT> places the
<TT>root</TT> node on the active list and calls <TT>walktree</TT> to
scan the tree.</P>

<P><A HREF="treegrav.c.html#walktree"><B>walktree</B></A> is the main
recursive routine for force calculation.  Its prototype is
<PRE>
    void walktree(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                  nodeptr p, real psize, vector pmid);
</PRE>
The effect of <TT>walktree</TT> is to compute gravity on bodies within
node <TT>p</TT>.  This is accomplished via a recursive scan of
<TT>p</TT> and its descendents.  At each point in the scan,
information from levels between the root and <TT>p</TT> is contained
in a set of nodes which will appear on the final interaction lists of
all bodies within <TT>p</TT>; this set is split into separate lists of
cells and bodies, addressed by <TT>cptr</TT> and <TT>bptr</TT>,
respectively.  The rest of the tree is represented by a set of active
nodes which include node <TT>p</TT> and surround it in space; pointers
to these nodes are stored in an array between <TT>aptr</TT> and
<TT>nptr</TT>.  Node <TT>p</TT> has linear size <TT>psize</TT> and
geometric midpoint <TT>pmid</TT>.</P>

<P>In the main loop of the routine, <TT>walktree</TT> examines the
active nodes, deciding which may be appended to the interaction lists,
and which are so close that their descendents must be examined at the
next level of recursion.  Cells are tested by the function
<TT>accept</TT> to determine if they are sufficiently well-separated
from <TT>p</TT>.  If so, they are placed on the cell interaction list,
headed by <TT>cp</TT>; if not, their <I>descendents</I> are placed on
the new active list, headed by <TT>np</TT>.  Bodies, unless they
happen to be the body <TT>p</TT> itself, are placed directly on the
body interaction list, which is headed by <TT>bp</TT>.</P>

<P>(It may be worth testing the type of a node before placing it on
the new active list; that way, bodies can be copied once and for all
to the body interaction list, instead of being copied again for each
immediate descendent of <TT>p</TT>.  But this will complicate the
handling of self-interaction.)</P>

<P>If some active cells were rejected by <TT>accept</TT>, recursion
continues to the next level of the tree, taking as active the
descendents of the rejected cells.  The actual recursive descent is
performed by <TT>walksub</TT>.  Otherwise, <TT>p</TT> points to a
body, and the interaction lists are complete, so <TT>gravsum</TT> may
be called to sum up the interactions.</P>

<P><A HREF="treegrav.c.html#accept"><B>accept</B></A> determines if a
cell passes the separation test.  Its prototype is
<PRE>
    bool accept(nodeptr c, real psize, vector pmid);
</PRE>
where <TT>c</TT> is the cell under consideration, and <TT>psize</TT>
and <TT>pmid</TT> specify the volume represented by the current node
<TT>p</TT>.</P>

<P>Two versions of <TT>accept</TT> exist.  The default version insures
that the squared distance from <TT>c</TT> to the nearest point within
<TT>p</TT> is greater than <TT>Rcrit2(c)</TT>.  But if the
<TT>QUICKSCAN</TT> version is compiled, a minimal criterion is used;
no cell touching <TT>p</TT>, even at a single corner, is accepted.</P>

<P><A HREF="treegrav.c.html#walksub"><B>walksub</B></A> performs the
actual recursive call-back to <TT>walktree</TT>.  Its prototype is
<PRE>
    void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                 nodeptr p, real psize, vector pmid);
</PRE>
and all parameters have exactly the same values that they do at the
calling point in <TT>walktree</TT>.</P>

<P>Two possible cases arise in <TT>walksub</TT>.  Most often by far,
the node <TT>p</TT> is a cell.  In this case, <TT>walksub</TT> loops
over its descendents, invoking <TT>walktree</TT> for each, with the
appropriately shifted cell center.  Much more rarely, <TT>p</TT> is
actually a body.  In this case, the active list contains nodes which
can't be sorted into the interaction lists at the present level;
<TT>walksub</TT> calls <TT>walktree</TT> exactly once, continuing the
search to the next level by `virtually' extending the tree.</P>

<P><A HREF="treegrav.c.html#gravsum"><B>gravsum</B></A> supervises the
process of evaluating forces from the cell and body interaction
lists.  Its prototype is
<PRE>
    void gravsum(bodyptr p0, cellptr cptr, cellptr bptr);
</PRE>
where <TT>p0</TT> is the body requiring updated forces and
<TT>cptr</TT> and <TT>bptr</TT> are pointers to the cell and body
interaction lists, respectively.</P>

<P>The main job of <TT>gravsum</TT> is to call the worker routines
<TT>sumnode</TT> and <TT>sumcell</TT> to sum up interaction terms.
The latter routine is only invoked if quadrupole moments are
included.</P>

<P><A HREF="treegrav.c.html#sumnode"><B>sumnode</B></A> sums up
interactions without quadrupole corrections.  Its prototype is
<PRE>
    void sumnode(cellptr start, cellptr finish,
                 vector pos0, real *phi0, vector acc0);
</PRE>
where <TT>start</TT> and <TT>finish</TT> point to the front and back
of the interaction list, <TT>pos0</TT> is the place where the force is
evaluated, and <TT>phi0</TT> and <TT>acc0</TT> are the resulting
potential and acceleration.</P>

<P><A HREF="treegrav.c.html#sumcell"><B>sumcell</B></A> sums up
interactions, including quadrupole corrections.  Its prototype is
<PRE>
    void sumcell(cellptr start, cellptr finish,
                 vector pos0, real *phi0, vector acc0);
</PRE>
where the arguments have the same meanings they do for
<TT>sumnode</TT>.</P>

<P>This routine is similar to <TT>sumnode</TT>, but includes
quadrupole-moment corrections (Hernquist 1987) to improve the forces
and potentials generated by cells.</P>

<P><B>NOTE:</B> <TT>sumnode</TT> and <TT>sumcell</TT> together account
for upwards of 90% of the cycles in typical calculations.  Optimizing
the performance of these routines is critical.</P>

<H3>3.3. Control &amp; Integration</H3>

<P>The routines governing the N-body integration are collected in <A
HREF="treecode.c.html"><TT>treecode.c</TT></A>.  Parameters and state
variables associated with control and integration are defined in <A
HREF="treecode.h.html"><TT>treecode.h</TT></A>.</P>

<P><A HREF="treecode.c.html#main"><B>main</B></A> is the main controlling
routine.  Its prototype is
<PRE>
    void main(int argc, string argv[]);
</PRE>

<P><A HREF="treecode.c.html#treeforce"><B>treeforce</B></A> is the
supervisory routine for force calculation.  Its prototype is
<PRE>
    void treeforce(void);
</PRE>

<P><A HREF="treecode.c.html#stepsystem"><B>stepsystem</B></A> advances
the system by one time-step.  Its prototype is
<PRE>
    void stepsystem(void);
</PRE>
Integration is performed using a `synchronized leap-frog', which is
computationally equivalent to the traditional time-centered leap-frog
but retains the synchronization of positions and velocities.  The
formula to advance the positions <TT><B>r</B></TT> and velocities
<TT><B>v</B></TT> from step <TT>n</TT> to step <TT>n+1</TT> is:
<PRE><TT>
                <B>v</B><SUB>n+1/2</SUB> = <B>v</B><SUB>n    </SUB> + (</TT><I>h</I><TT>/2) <B>a</B><SUB>n    </SUB> ,

                <B>r</B><SUB>n+1  </SUB> = <B>r</B><SUB>n    </SUB> +   </TT><I>h</I><TT>   <B>v</B><SUB>n+1/2</SUB> ,

                <B>v</B><SUB>n+1  </SUB> = <B>v</B><SUB>n+1/2</SUB> + (</TT><I>h</I><TT>/2) <B>a</B><SUB>n+1  </SUB> ,
</TT></PRE>
where <TT><B>a</B> = <B>a</B>(<B>r</B>)</TT> is the acceleration
computed from positions at the corresponding step, and <I>h</I>
is the time-step.</P>

<P><A HREF="treecode.c.html#startrun"><B>startrun</B></A> initializes
parameters and data for the N-body run.  Its prototype is
<PRE>
    void startrun(void);
</PRE>

<P><A HREF="treecode.c.html#testdata"><B>testdata</B></A> sets up a
Plummer model (Aarseth, Henon, &amp; Wielen 1974).  Its prototype is
<PRE>
    void testdata(void);
</PRE>

<H3>3.4. Input/Output</H3>

<P>The routines responsible for input &amp; output of N-body data are
in <A HREF="treeio.c.html"><TT>treeio.c</TT></A>.</P>

<P><A HREF="treeio.c.html#inputdata"><B>inputdata</B></A> reads initial
conditions from an input file.  Its prototype is
<PRE>
    void inputdata(void);
</PRE></P>

<P><A HREF="treeio.c.html#startoutput"><B>startoutput</B></A> prints a
header describing the calculation.  Its prototype is
<PRE>
    void startoutput(void);
</PRE></P>

<P><A HREF="treeio.c.html#forcereport"><B>forcereport</B></A> prints
statistics on tree construction and force calculation.  Its prototype
is
<PRE>
    void forcereport(void);
</PRE></P>

<P><A HREF="treeio.c.html#output"><B>output</B></A> prints diagnostics
and determines if a data output is required.  Its prototype is
<PRE>
    void output(void);
</PRE></P>

<P><A HREF="treeio.c.html#outputdata"><B>outputdata</B></A> outputs the
actual N-body data.  Its prototype is
<PRE>
    void outputdata(void);
</PRE></P>

<P><A HREF="treeio.c.html#diagnostics"><B>diagnostics</B></A> computes
various dynamical diagnostics.  Its prototype is
<PRE>
    local void diagnostics(void);
</PRE></P>

<P><A HREF="treeio.c.html#savestate"><B>savestate</B></A> writes current
program state to a binary file.  Its prototype is
<PRE>
    void savestate(string pattern);
</PRE></P>

<P><A HREF="treeio.c.html#restorestate"><B>restorestate</B></A> restores
the program state from a binary file.  Its prototype is
<PRE>
    void restorestate(string file);
</PRE></P>

<A NAME="instructions">
<H2>4. Instructions</H2>

<H3>4.1. Distribution</H3>

<P>A complete version of <B>treecode</B> may be down-loaded from
this web site.  The entire code is available as a single gzipped tar
file:
<UL>
    <LI> <A HREF="treecode.tar.gz"><TT>treecode.tar.gz</TT></A>
</UL>
Alternatively, individual files may be saved as plain text (with the
appropriate <TT>.c</TT> and <TT>.h</TT> extensions).  The files
required are
<A HREF="treecode.h.html"><TT>treecode.h</TT></A>,
<A HREF="treedefs.h.html"><TT>treedefs.h</TT></A>,
<A HREF="treecode.c.html"><TT>treecode.c</TT></A>,
<A HREF="treegrav.c.html"><TT>treegrav.c</TT></A>,
<A HREF="treeio.c.html"><TT>treeio.c</TT></A>, and
<A HREF="treeload.c.html"><TT>treeload.c</TT></A>.
Some general-purpose include files and library routines are also
needed; these are found in
<A HREF="getparam.h.html"><TT>getparam.h</TT></A>,
<A HREF="mathfns.h.html"><TT>mathfns.h</TT></A>,
<A HREF="stdinc.h.html"><TT>stdinc.h</TT></A>,
<A HREF="vectdefs.h.html"><TT>vectdefs.h</TT></A>, 
<A HREF="vectmath.h.html"><TT>vectmath.h</TT></A>,
<A HREF="clib.c.html"><TT>clib.c</TT></A>,
<A HREF="getparam.c.html"><TT>getparam.c</TT></A>, and
<A HREF="mathfns.c.html"><TT>mathfns.c</TT></A>.
A <A HREF="Makefile"><TT>Makefile</TT></A> is also provided to help
organize the compilation process.</P>

<P>The <B>treecode</B> sources are provided as free software with no
restrictions on their use or modification.  You may redistribute this
software at will provided that you do so without restricting its use,
modification, or redistribution by others.  See the <A
HREF="http://www.gnu.org/copyleft/gpl.html">GNU General Public
License</A> for details.</P>

<P>This software is provided `as is'; no warranty is expressed or
implied.  Naturally, I have run extensive tests of this software
before releasing it; just as naturally, there are probably bugs or
limitations not identified in these tests.  If you encounter problems
in using this software, please let me know!</P>

<P>If you publish results obtained using <B>treecode</B>, I would
appreciate a brief acknowledgment (and a preprint or reprint).  I will
not be able to improve this software and continue making it available
without adequate research support; acknowledging the use of this
software is an investment in future developments as well as a basic
courtesy.</P>

<H3>4.2. Compilation</H3>

<P><B>treecode</B> is written in ANSI C; while some I/O operations may
depend on the operating system, the bulk of the code should work on
any system with an ANSI C compiler.  The instructions below assume
that a UNIX-like operating system is available.</P>

<P>Begin by placing the files in a single directory.  If you have
copied the <A HREF="treecode.tar.gz"><TT>treecode.tar.gz</TT></A>
file, this may be done by giving the following commands:
<PRE>
    % gunzip treecode.tar.gz
    % tar xvf treecode.tar
</PRE>
The directory should then contain all of the <TT>.h</TT> and
<TT>.c</TT> files listed above, as well as the <TT>Makefile</TT>.</P>

<P>Before actually compiling the code, you may need to edit the
<TT>Makefile</TT> and modify some of the parameters or switches.
Defaults are provided for SGI, Sun, and LINUX compilers; suggestions
for other operating system/compiler combinations are welcome.</P>

<P>Both single-precision and double-precision versions of the code can
be generated by changing the value of the <TT>PRECISION</TT> parameter
in the <TT>Makefile</TT>.  Single-precision is good enough for most
calculations, so that is the default.  Some run-time libraries,
however, do not provide single-precision versions of floating-point
functions.  On such systems, use the <TT>MIXEDPREC</TT> option for
best performance.  Note that if you change the <TT>PRECISION</TT>
parameter, you must delete any existing <TT>.o</TT> and <TT>.a</TT>
files before recompiling.</P>

<P>Once the <TT>Makefile</TT> has been edited, build a version of the
code by giving either of these commands:
<PRE>
    % make treecode
    % make treecode_q
</PRE>
The first of these builds the standard version, while the second
builds the <TT>QUICKSCAN</TT> version.</P>

<H3>4.3. Operation</H3>

<P>A test calculation using a self-generated Plummer model (Aarseth,
Henon, &amp; Wielen 1974) may be run by giving the either of these
commands
<PRE>
    % treecode
    % treecode_q
</PRE>
This calculation uses <TT>nbody=4096</TT> bodies, a timestep of
<TT>dtime=1/32</TT>, a smoothing length of <TT>eps=0.025</TT>, an
accuracy parameter <TT>theta=1.0</TT> (for the standard version), and
no quadrupole moments.  A log of the calculation is sent to standard
output; no other output is generated.  On a 500 MHz Pentium, this test
calculation takes about 1.2 minutes.</P>

<P><B>treecode</B> uses the <A HREF="#getparam"><TT>getparam</TT></A>
command-line package to get all input parameters, including the names
of any input and output files.  To see a full list of parameters,
along with their default values and descriptions, give either of these
commands:
<PRE>
    % treecode -help
    % treecode_q -help
</PRE>
In response, the standard <B>treecode</B> prints out the following:
<PRE>
  treecode                          Hierarchical N-body code (theta scan)
    in=                             Input file with initial conditions
    out=                            Output file of N-body frames
    dtime=1/32                      Leapfrog integration timestep
    eps=0.025                       Density smoothing length
    theta=1.0                       Force accuracy parameter
    usequad=false                   If true, use quad moments
    options=                        Various control options
    tstop=2.0                       Time to stop integration
    dtout=1/4                       Data output timestep
    nbody=4096                      Number of bodies for test run
    seed=123                        Random number seed for test run
    save=                           Write state file as code runs
    restore=                        Continue run from state file
    VERSION=1.4                     Joshua Barnes  February 21 2001
</PRE>
(the <TT>QUICKSCAN</TT> version produces a similar listing with one
less parameter).  This printout lists the name of each parameter, its
default value if any, and a brief explanation of its function.  The
<TT>getparam</TT> package accepts argument values on the command line,
and matches them to parameters either by position or by name.
Initially, positional matching is used: the first argument is matched
to the first parameter, and so on.  However, if any argument has the
form <I>name</I><TT>=</TT><I>value</I>, that argument and any that
follow it are matched by name.  An error results if a <I>name</I> does
not match any parameter or if more than one value is assigned to a
parameter.</P>

<P>At somewhat greater length than above, the parameters accepted by
<B>treecode</B> are interpreted as follows:
<UL>

    <LI> <B>in</B>, if given, names an input file containing initial
         conditions for the calculation.  The format of this file is
         described below.

    <LI> <B>out</B>, if given, is a pattern naming output files for
         N-body snapshots taken at regular intervals.  This pattern is
         used as an argument to <TT>sprintf</TT>, along with the
         integration step number, to generate an actual file name.  If
         the resulting file already exists, the new data is appended.
         The format used is similar to the format used for input
         files.

    <LI> <B>dtime</B> is the integration time-step.  It's convenient
         to use a timestep which has an exact representation as a
         floating-point number; for example, a value <I>n/d</I>, where
         <I>n</I> is an integer and <I>d</I> is a power of two.  To
         simplify specification of such timesteps, the code accepts
         fractions as well as ordinary floating-point values.  If
         <TT>dtime=0</TT>, <B>treecode</B> does a single force
         calculation, outputs the result, and exits.

    <LI> <B>eps</B> is the smoothing length used in the gravitational
         force calculation.  In effect, the mass distribution is
         smoothed by replacing each body by a Plummer sphere with
         scale length <TT>eps</TT>, and the gravitational field of
         this smoothed distribution is calculated.

    <LI> <B>theta</B> is the <I>opening angle</I> used to adjust the
         accuracy of the force calculation.  Values less than unity
         produce more accurate forces, albeit at greater computational
         expense.  The <TT>QUICKSCAN</TT> version does not use or
         accept this parameter.

    <LI> <B>usequad</B> determines if quadrupole corrections are used
         when calculating gravitational fields.  These corrections can
         significantly improve the accuracy of force calculation at a
         fairly modest computational cost.

    <LI> <B>options</B> is a list of comma-separated key words used to
         obtain various run-time options.  The options available are

         <UL>
             <LI> <B>reset-time</B>: set the time to zero,
                  regardless of the value in the input file; 
             <LI> <B>new-tout</B>: reschedule the first snapshot
                  output;
             <LI> <B>out-phi</B>: include potential values in the
                  output files;
             <LI> <B>out-acc</B>: include acceleration vectors in
                  the output files;
             <LI> <B>bh86</B>: use the BH86 criterion to set
                  critical cell radii;
             <LI> <B>sw94</B>: use Salmon &amp; Warren's (1994)
                  <I>bmax</I> criterion to set critical cell radii.
         </UL>

    <LI> <B>tstop</B> is the time at which the N-Body integration
         terminates.

    <LI> <B>dtout</B> is the time interval between output files.  To
         insure that outputs are performed when expected,
         <TT>dtout</TT> should be a multiple of <TT>dtime</TT>.  Like
         <TT>dtime</TT>, this parameter may be specified as a
         fraction.

    <LI> <B>nbody</B> is the number of bodies used to self-generate
         initial conditions.  This parameter is only used if
         no input or restore file is given.

    <LI> <B>seed</B> is the random number seed used to self-generate
         initial conditions.  This parameter is only used if
         no input or restore file is given.

    <LI> <B>save</B>, if given, is a pattern naming binary files used
         to record the state of the code after each step.  This
	 pattern is used as an argument to <TT>sprintf</TT>, along
	 with the lowest bit of the step number, to construct the
	 actual file name.

    <LI> <B>restore</B> names a binary file written using the
         <TT>save</TT> function above.  The calculation resumes from
         that point.  New values of the <TT>eps</TT>, <TT>theta</TT>,
         <TT>usequad</TT>, <TT>options</TT>, <TT>tstop</TT>, and
         <TT>dtout</TT> parameters may be supplied; if not, the values
         from the saved calculation are used.  If N-body data output
         is required, a new output file pattern <TT>out</TT> must be
         given.

</UL>

<P>For example, to run a test calculation using a Plummer model with
<TT>32768</TT> bodies and an opening angle of <TT>0.75</TT>, type
<PRE>
    % treecode nbody=32768 theta=0.75
</PRE>
To compute forces and potentials for a N-body configuration in the
input file <TT>input.data</TT>, writing the result to output file
<TT>forces.data</TT>, type
<PRE>
    % treecode input.data forces.data dtime=0 options=out-phi,out-acc
</PRE>
To run the initial conditions in <TT>input.data</TT> with a timestep of
<TT>1/64</TT>, writing results at intervals of <TT>1/16</TT> to
separate files <TT>run_000.data</TT>, <TT>run_004.data</TT>, ...,
type
<PRE>
    % treecode input.data run_%03d.data dtime=1/64 dtout=1/16
</PRE>
To perform the same calculation using the <TT>QUICKSCAN</TT> version,
while saving the program's state after each step, type
<PRE>
    % treecode_q input.data run_%03d.data dtime=1/64 dtout=1/16 save=state.data
</PRE>
To continue this calculation until time <TT>4</TT>, type
<PRE>
    % treecode_q restore=state.data out=run_%03d.data tstop=4
</PRE>
</P>

<H3>4.4. File Formats</H3>

<P>By default, the input and output files used by <B>treecode</B> are
written in ASCII.  Each file contains one or more snapshots.  Input
files have the following format:
<PRE>
    nbody
    NDIM
    time
    mass[1]
    ....
    mass[n]
    x[1] y[1] z[1]
    ....
    x[n] y[n] z[n]
    vx[1] vy[1] vz[1]
    ....
    vx[n] vy[n] vz[n]
</PRE>
Thus an input file contains a total of <TT>3+3*nbody</TT> lines.  A
similar format is used for output files.  If the <TT>out-phi</TT>
and/or <TT>out-acc</TT> options are set, velocities are followed by
potential values and/or acceleration vectors, respectively.</P>

<P>While ASCII output is easy to read, it is relatively inefficient.
If <A HREF="treeio.c.html"><TT>treeio.c</TT></A> is compiled with the
<TT>BINARYIO</TT> preprocessor symbol defined, input and output are
performed in binary.  This is recommended for production
calculations!</P>

<H2>Thanks</H2>

<P>The bucolic scenery of Leiden inspired the invention of this
algorithm and I thank Tim de Zeeuw for his hospitality.  I'm also
grateful to Jun Makino for hospitality while I developed the public
version described here.  Initial development of <B>treecode</B> was
supported by NASA grant NAG-2836.</P>

<A NAME="appendix_a"></A>
<H2>Appendix A. Changes From Version 1.3</H2>

<P>The main change is that the previous version expected the leapfrog
stepsize to be specified as a frequency; instead of the parameter
<TT>dtime</TT>, Version 1.3 accepted a parameter <TT>freq</TT> equal
to the inverse of the timestep.  I decided to change this because many
people find an integration timestep more intuitive than an integration
frequency.  It seemed convenient to use a frequency since power-of-two
timesteps could be easily specified; for example, <TT>freq=32</TT>.
However, the new version allows the user to give the timestep as a
fraction; for example, <TT>dtime=1/32</TT>.  The new version of the
code also uses <TT>dtout</TT> to specify the time between data
outputs, instead of <TT>freqout</TT> to specify the frequency of
outputs.</P>

<P>If you prefer to continue using <TT>freq</TT> and <TT>freqout</TT>,
you can compile the code with the <TT>USEFREQ</TT> preprocessor
symbol; see the <A HREF="Makefile"><TT>Makefile</TT></A> for details.
Note that if you want to resume calculations from state files
generated with Version 1.3, you <B>must</B> compile Version 1.4 with
<TT>USEFREQ</TT>.</P>

<P>In addition, the <A HREF="mathfns.h.html"><TT>mathfns.h</TT></A>
file in Version 1.4 has additional definitions to use alternate naming
conventions for single-precision math functions.  This is useful when
compiling under some versions of LINUX.</P>

<P>I no longer have access to an IBM AIX system for testing, so I have
removed a compilation option specific to that operating system.
Earlier tests indicated that <TT>sumnode</TT> and <TT>sumcell</TT> ran
more efficiently under IBM AIX if the local variables in these
routines were declared <TT>double</TT> instead of <TT>real</TT>.  I do
not know if this is still the case.</P>

<A NAME="appendix_b"></A>
<H2>Appendix B. Zeno System Software</H2>

<P>The <B>treecode</B> software depends on a number of include files
and function libraries, collectively known as the `Zeno System'.
These files define keywords, data types, and library routines useful
in C programs.</P>

<H3>A.1. Standard Include File</H3>

<P><A HREF="stdinc.h.html"><B>stdinc.h</B></A> is a standard include
file.  It defines some constants and data types, and provides
prototypes for a few functions from a general-purpose C library.  The
following constants and types are used in <B>treecode</B>.</P>

<UL>

<LI><B>NULL</B> indicates a pointer to no object.  This definition is
identical to the one in C's I/O include file <TT>stdio.h</TT>; it is
included here for completeness.

<LI><B>local</B> is a synonym for <TT>static</TT>, used to indicate
that a file-level data object or routine is defined only within that
file.

<LI><B>bool</B> is a storage type for logical values.  Also defined
are the constants <TT>TRUE</TT> and <TT>FALSE</TT>.

<LI><B>string</B> is a storage type for a pointer to a null-terminated
sequence of characters.

<LI><A HREF="stdinc.h.html#real"><B>real</B></A> is a storage type for
real-valued numbers.  The precision of <TT>real</TT> values is fixed
at compile time; <TT>real</TT> is synonymous to <TT>float</TT> if
either <TT>SINGLEPREC</TT> or <TT>MIXEDPREC</TT> is defined, and to
<TT>double</TT> if <TT>DOUBLEPREC</TT> is defined.  For more details,
see the description of <TT>mathfns.h</TT> below.

</UL>

<P>The following functions, used in <B>treecode</B>, have prototypes
in <A HREF="stdinc.h.html"><TT>stdinc.h</TT></A>; their source code is
given in <A HREF="clib.c.html">clib.c</A>.</P>

<P><B>allocate</B> is a memory-allocation function.  Its prototype is
<PRE>
    void *allocate(int nbyte);
</PRE>
where <TT>nbyte</TT> is the number of bytes to allocate.  The
<TT>allocate</TT> function exits via <TT>error</TT> if it can't get
the requested amount of memory.  The memory is cleared before
<TT>allocate</TT> returns.</P>

<P><B>cputime</B> returns the total process CPU time in minutes.  Its
prototype is
<PRE>
    double cputime(void);
</PRE>
</P>

<P><B>error</B> reports an error and exits.  Its prototype is
<PRE>
    void error(string fmt, ...);
</PRE>
The <TT>fmt</TT> string and any further arguments are interpreted
exactly like the system routine <TT>printf</TT>, but output to
<TT>stderr</TT> instead of <TT>stdout</TT>.</P>

<P><B>scanopt</B> scans an option string for a keyword.  Its prototype
is
<PRE>
    bool scanopt(string opts, string word);
</PRE>
where <TT>opts</TT> is a series of words separated by commas, and
<TT>word</TT> is a single keyword; it returns <TT>TRUE</TT> if
<TT>word</TT> appears in <TT>opts</TT>, and <TT>FALSE</TT>
otherwise.</P>

<P><B>stropen</B> opens a <TT>stdio.h</TT> stream for input/output.
Its prototype is
<PRE>
    stream stropen(string file, string mode);
</PRE></P>

<H3>A.2. Real Functions</H3>

<P><A HREF="mathfns.h.html"><B>mathfns.h</B></A> defines <TT>real</TT>
versions of the math functions in <TT>math.h</TT>, along with some
extensions. Most function name are derived from their counterparts in
<TT>math.h</TT> by prepending the letter `<TT>r</TT>' (for
<TT>real</TT>).  The following functions are used by <B>treecode</B>;
sources appear in <A HREF="mathfns.c.html"><TT>mathfns.c</TT></A>.</P>

<P><B>rsqrt</B> computes a real-valued square root.  Its prototype is
<PRE>
    real rsqrt(real x);
</PRE></P>

<P><B>rsqr</B> computes the inverse of <TT>rsqrt</TT>.  Its prototype
is
<PRE>
    real rsqr(real x);
</PRE></P>

<P><B>rpow</B> computes a real-valued power.  Its prototype is
<PRE>
    real rpow(real x, real y);
</PRE></P>

<P><B>rabs</B> computes a real absolute value.  Its prototype is
<PRE>
    real rabs(real x);
</PRE></P>

<P><B>xrandom</B> generates a random number between specified limits.
Its prototype is
<PRE>
    double xrandom(double x1, double x2);
</PRE></P>

<P><B>pickshell</B> picks a point on the surface of a shell in an
n-dimensional space.  Its prototype is
<PRE>
    void pickshell(real *point, int ndim, real radius);
</PRE></P>

<A NAME="getparam"></A>
<H3>A.3. Getparam Package</H3>

<P><A HREF="getparam.h.html"><B>getparam.h</B></A> defines a general
command-line argument processing package.  Source code is given in <A
HREF="getparam.c.html"><TT>getparam.c</TT></A>.  This package provides
the client program with a simple self-documenting user interface.  The
user may specify arguments by position or by keyword; defaults are
supplied for all arguments not specified.  Used in <B>treecode</B> are
the following.</P>

<P><B>initparam</B> initializes the command-line processing package.
Its prototype is
<PRE>
    void initparam(string argv[], string defv[]);
</PRE>
where <TT>argv</TT> is the <TT>NULL</TT>-terminated vector of
command-line arguments supplied to <TT>main</TT>, and <TT>defv</TT> is
a vector of parameter names, default values, and documentation
messages.</P>

<P><B>getparam</B> returns the value of a parameter.  Its prototype
is
<PRE>
    string getparam(string name);
</PRE>
where <TT>name</TT> is the name of the parameter as supplied to
<TT>initparam</TT>.  An error occurs if <TT>name</TT> is not found.</P>

<P><B>getiparam</B> returns the integer value of a parameter.  Its
prototype is
<PRE>
    int getiparam(string name);
</PRE></P>

<P><B>getdparam</B> returns the double-precision floating-point value
of a parameter.  Its prototype is
<PRE>
    double getdparam(string name);
</PRE></P>

<P><B>getbparam</B> returns the boolean value of a parameter.  Its
prototype is
<PRE>
    bool getbparam(string name);
</PRE></P>

<P><B>getparamstat</B> returns an indication of parameter status.  Its
prototype is
<PRE>
    int getparamstat(string name);
</PRE></P>

<H3>A.4. Vectors and Matrices</H3>

<P><A HREF="vectdefs.h.html"><B>vectdefs.h</B></A> defines vectors and
matrices.  The number of dimensions is a compile-time constant,
specified by the preprocessor symbols <TT>THREEDIM</TT>,
<TT>TWODIM</TT>, and <TT>NDIM</TT>; the default is <TT>THREEDIM</TT>.
The types defined are as follows.</P>

<UL>

<LI><B>vector</B> is a storage type for a vector of <TT>NDIM</TT>
<TT>real</TT> values.

<LI><B>matrix</B> is a storage type for a matrix of
<TT>NDIM<SUP>2</SUP></TT> <TT>real</TT> values.

</UL>

<P><A HREF="vectmath.h.html"><B>vectmath.h</B></A> defines a
collection of macros for operating on vectors and matrices.  Those
used by <B>treecode</B> are listed below.</P>

<UL>

<LI><B>ABSV(</B><I>s</I><B>,v)</B> computes the magnitude of a vector:

    <P>&nbsp;&nbsp;<I>s = sqrt(v<SUB>i</SUB> v<SUB>i</SUB>)</I></P>

<LI><B>ADDM(p,q,r)</B> adds matrices:

    <P>&nbsp;&nbsp;<I>p<SUB>ij</SUB> = q<SUB>ij</SUB> + r<SUB>ij</SUB></I></P>

<LI><B>ADDV(u,v,w)</B> adds vectors:

    <P>&nbsp;&nbsp;<I>u<SUB>i</SUB> = v<SUB>i</SUB> + w<SUB>i</SUB></I></P>

<LI><B>CLRM(p)</B> sets the elements of a matrix to zero:

    <P>&nbsp;&nbsp;<I>p<SUB>ij</SUB> = 0</I></P>

<LI><B>CLRV(v)</B> sets the elements of a vector to zero:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = 0</I></P>

<LI><B>DISTV(</B><I>s</I><B>,u,v)</B> computes the distance between
two vectors:

    <P>&nbsp;&nbsp;<I>s = sqrt((v<SUB>i</SUB> - u<SUB>i</SUB>)(v<SUB>i</SUB> -
                u<SUB>i</SUB>))</I></P>

<LI><B>DIVVS(v,u,</B><I>s</I><B>)</B> divides a vector by a scalar:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = u<SUB>i</SUB> / s</I></P>

<LI><B>DOTVP(</B><I>s</I><B>,v,u)</B> forms the dot product of two
vectors:

    <P>&nbsp;&nbsp;<I>s = v<SUB>i</SUB> u<SUB>i</SUB></I></P>

<LI><B>MULMS(p,q,</B><I>s</I><B>)</B> multiplies a matrix by a scalar:

    <P>&nbsp;&nbsp;<I>p<SUB>ij</SUB> = q<SUB>ij</SUB> s</I></P>

<LI><B>MULVS(v,u,</B><I>s</I><B>)</B> multiplies a vector by a scalar:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = u<SUB>i</SUB> s</I></P>

<LI><B>OUTVP(p,v,u)</B> forms the outer product of two vectors:

    <P>&nbsp;&nbsp;<I>p<SUB>ij</SUB> = v<SUB>i</SUB> u<SUB>j</SUB></I></P>

<LI><B>SETM(p,q)</B> sets one matrix equal to another:

    <P>&nbsp;&nbsp;<I>p<SUB>ij</SUB> = q<SUB>ij</SUB></I></P>

<LI><B>SETV(v,u)</B> sets vector equal to another:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = u<SUB>i</SUB></I></P>

<LI><B>SETVS(v,</B><I>s</I><B>)</B> sets all components of a vector
equal to a scalar:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = s</I></P>

<LI><B>SETMI(p)</B> sets a matrix to the identity:

    <P>&nbsp;&nbsp;<I>p<SUB>ij</SUB> = delta<SUB>ij</SUB></I></P>

<LI><B>SUBM(p,q,r)</B> subtracts matrices:

    <P>&nbsp;&nbsp;<I>p<SUB>ij</SUB> = q<SUB>ij</SUB> - r<SUB>ij</SUB></I></P>

<LI><B>SUBV(v,u,w)</B> subtracts vectors:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = u<SUB>i</SUB> - w<SUB>i</SUB></I></P>

</UL>

<P>There are also some macros defining compound operations specific to
<B>treecode</B>.  These are written so as to help the compiler produce
good code in the force summation loops.</P>

<UL>

<LI><B>ADDMULVS(v,u,</B><I>s</I><B>)</B> scales a vector, and adds the
result to another:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = v<SUB>i</SUB> + u<SUB>i</SUB> s</I></P>

<LI><B>ADDMULVS2(v,u,</B><I>s</I><B>,w,</B><I>r</I><B>)</B> scales two
vectors, and adds the result to a third:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = v<SUB>i</SUB> + u<SUB>i</SUB> s
    + w<SUB>i</SUB> r</I></P>

<LI><B>DOTPMULMV(</B><I>s</I><B>,v,p,u)</B> multiplies a vector by a
matrix, and also forms the dot product:

    <P>&nbsp;&nbsp;<I>v<SUB>j</SUB> = p<SUB>ij</SUB> u<SUB>i</SUB></I>
    &nbsp;&nbsp;and&nbsp;&nbsp; <I>s = v<SUB>j</SUB> v<SUB>j</SUB></I></P>

<LI><B>DOTPSUBV(</B><I>s</I><B>,v,u,w)</B> subtracts two vectors, and
forms the dot product:

    <P>&nbsp;&nbsp;<I>v<SUB>i</SUB> = u<SUB>i</SUB> - w<SUB>i</SUB></I>
    &nbsp;&nbsp;and&nbsp;&nbsp; <I>s = v<SUB>i</SUB> v<SUB>i</SUB></I></P>

</UL>

<HR>

<H2>References</H2>

<UL>

    <LI>Aarseth, S.J., Henon, M., &amp; Wielen, R. 1974 <I>Astr. &amp
        Ap.</I> <B>37</B>, 183.
    <LI>Barnes, J.E. 1986. In <I>The Use of Supercomputers in Stellar
        Dynamics</I>, eds. P. Hut and S. McMillan (Berlin:
        Springer-Verlag), p. 175.
    <LI>Barnes, J.E. 1990. <I>J. Comp. Phys.</I> <B>87</B>, 161.
    <LI>Barnes, J.E. 1995. In <I>The Formation of Galaxies</I>, eds. C.
        Munoz-Tunon &amp; F. Sanchez (Cambridge: Cambridge University
        Press), p. 399.
    <LI>Barnes, J.E. 1998. <I>Dynamics of Galaxy Interactions</I>, in
        <I>Galaxies: Interactions and Induced Star formation</I>, by
        R.C. Kennicutt, Jr., F. Schweizer, &amp; J.E. Barnes (Berlin:
        Springer).
    <LI>Barnes, J. &amp; Hut, P. 1986. <I>Nature</I> <B>324</B>, 446.
    <LI>Greengard, L. 1990. <I>Computers in Physics</I>, <B>4</B>, 142.
    <LI>Greengard, L. &amp; Rokhlin, V. 1987. <I>J. Comp. Phys.</I>
        <B>73</B>, 325.
    <LI>Hernquist, L. 1987. <I>Ap. J. Suppl.</I> <B>64</B>, 715.
    <LI>Salmon, J.K. &amp; Warren, M.S. 1994. <I>J. Comp. Phys.</I>
        <B>111</B>, 136 .
    <LI>Sugimoto, D. <I>et al.</I> 1990. <I>Nature</I> <B>345</B>, 33.
</UL>

<HR>

<P><ADDRESS><A HREF="../barnes.html">Joshua E. Barnes</A>
(<A HREF="mailto:barnes@ifa.hawaii.edu">barnes@ifa.hawaii.edu</A>)</ADDRESS>
Last modified: February 23, 2001</P>

<P><TT>http://www.ifa.hawaii.edu/~barnes/treecode/treeguide.html</TT></P>

</BODY>

</HTML>
