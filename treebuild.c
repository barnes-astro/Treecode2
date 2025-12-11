/*
 * treebuild.c: routines to create tree structure.
 * Copyright (c) 2025 by Joshua E. Barnes, Honolulu, Hawaii.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include <assert.h>

#include "treedefs.h"

//  GCPos1: geometric center of cell in tree coordinate system.
//  Here and below, vectors with a "1" suffix are in the tree CS.
//  _____________________________________________________________
 
#define GCPos1(x)  Pos(x)			// reuse cell pos. vector

//  Local routines and variables to perform tree construction.
//  __________________________________________________________
 
local void pickrotation(matrix rmat);
local void fitroot(treeptr tptr, bodyptr btab, int nbody);
local void loadbody(treeptr tptr, bodyptr p);
local int subindex(vector pos, cellptr q);
local void threadtree(nodeptr p, nodeptr q, int lev);
local void evalcells(cellptr p, real psize, treeptr tptr);
local void setrcrit(cellptr p, vector dr1, real psize, treeptr tptr);
local void tracequads(treeptr tptr);
local cellptr makecell(void);

local bool bh86;				// original opening criterion
local bool debug;				// flag messages, as needed

#define MAXLEVEL  50				// max height of tree

#define SCALERNG  0.5				// range of scale factors

//  treebuild: initialize tree structure for force calculation.
//  ___________________________________________________________
 
treeptr treebuild(bodyptr btab, int nbody, real eps, real theta, bool quad,
		  bool rotate, bool scale, real translate, string options)
{
  double cpustart = cputime();			// record time at start
  treeptr tptr;
  
  tptr = (treeptr) allocate(sizeof(tree));	// allocate new tree
  Theta(tptr) = theta;				// save key parameters
  Eps(tptr) = eps;
  Usequad(tptr) = quad;
  Options(tptr) = options;
  if (rotate)					// pick random rotation
    pickrotation(Rotn(tptr));
  else						// use default tree CS  
    SETMI(Rotn(tptr));
  if (scale) {					// pick scale factor
    Scale(tptr) = rpow(2.0, xrandom(-SCALERNG, SCALERNG));
    eprintf("[%s.treebuild: tree scale = %f]\n", getprog(), Scale(tptr));
  } else					// use default scale factor
    Scale(tptr) = 1.0;
  if (translate > 0.0) {			// pick random translation
    pickball(Trans(tptr), NDIM, translate);
    eprintf("[%s.treebuild: tree origin = %f %f %f]\n", getprog(),
	    Trans(tptr)[0], Trans(tptr)[1], Trans(tptr)[2]);
  } else					// use default translation
    CLRV(Trans(tptr));
  Nbody(tptr) = nbody;				// store number of bodies
  Root(tptr) = makecell();			// allocate the root cell
  Ncell(tptr) = 1;				// track cells used
  Nlevel(tptr) = 0;				// init max level
  CLRV(GCPos1(Root(tptr)));			// root at 0,0,0 in tree CS
  fitroot(tptr, btab, nbody);			// expand to fit bodies
  eprintf("[%s.treebuild: starting loadbody]\n", getprog());
  for (bodyptr p = btab; p < btab+nbody; p++) {	// loop over all bodies
    loadbody(tptr, p);				// insert each into tree
  }
  eprintf("[%s.treebuild: finished loadbody]\n", getprog());
  threadtree((nodeptr) Root(tptr), NULL, 0);	// convert to threaded tree
  if (scanopt(options, "rawdump"))
    treedump(tptr, stdout);
  bh86 = scanopt(options, "bh86");		// check original criterion
  evalcells(Root(tptr), Rsize(tptr), tptr);
  tracequads(tptr);
  CPUbuild(tptr) = cputime() - cpustart;	// store elapsed CPU time
  return (tptr);
}

//  treedestroy: reclaim threaded tree storage after force calculation.
//  ___________________________________________________________________

void treedestroy(treeptr tptr)
{
  nodeptr p = (nodeptr) Root(tptr), q;		// start at root of tree

  while (p != NULL)				// loop following thread
    if (Cell(p)) {				// if we reached a cell
      assert((Type(p) & ThreadFlag) != 0);	// make sure it's threaded
      q = More(p);				// get link to descendents
      free(p);					// return to free memory
      p = q;					// scan descendents
    } else					// skip over bodies
      p = Next(p);				// by going on to the next
  free(tptr);					// lastly, free structure
}

//  treedump: output text representation of tree structure.
//  _______________________________________________________

void treedump(treeptr tptr, stream ostr)
{
  matrix r;
  nodeptr p;

  fprintf(ostr, "# %s %s %s %s %s\n", "theta", "quad",
	  "rsize", "ncell", "nlevel");
  fprintf(ostr, "%f %d %f %d %d\n", Theta(tptr), Usequad(tptr),
	  Rsize(tptr), Ncell(tptr), Nlevel(tptr));
  SETM(r, Rotn(tptr));
  fprintf(ostr, "# %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
	  "rot00", "rot01", "rot02", "rot10", "rot11", "rot12",
	  "rot20", "rot21", "rot22", "scale", "trans0", "trans1", "trans2");
  fprintf(ostr, "%f %f %f %f %f %f %f %f %f %f %f %f %f\n",
	  r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2],
	  r[2][0], r[2][1], r[2][2], Scale(tptr),
	  Trans(tptr)[0], Trans(tptr)[1], Trans(tptr)[2]);
  fprintf(ostr, "# %s %s %s %s %s %s %s %s %s %s %s\n",
	  "addr", "type", "level", "active", "next", "more",
	  "rcrit2", "mass", "x", "y", "z");
  for (p = (nodeptr) Root(tptr); p != NULL; )	// first dump cells
    if (Cell(p)) {
      assert((Type(p) & ThreadFlag) != 0);
      fprintf(ostr, "%p %02x %02x %04x %14p %14p %f %f %f %f %f\n",
	      p, Type(p), Level(p), Active(p), Next(p), More(p),
	      Rcrit2(p), Mass(p), Pos(p)[0], Pos(p)[1], Pos(p)[2]);
      p = More(p);
    } else
      p = Next(p);
  for (p = (nodeptr) Root(tptr); p != NULL; )	// next dump bodies
    if (Body(p)) {
      fprintf(ostr, "%p %02x %02x %04x %14p %14p %f %f %f %f %f\n",
	      p, Type(p), Level(p), Active(p), Next(p), NULL, 0.0,
	      Mass(p), Pos(p)[0], Pos(p)[1], Pos(p)[2]);
      p = Next(p);
    } else
      p = More(p);
}

//  pickrotation: generate random rotation matrix.
//  ______________________________________________

local void pickrotation(matrix rotn)
{
  vector xhat, vtmp, yhat, zhat;
  real ymag;

  pickshell(xhat, NDIM, 1.0);			// 1st vector points anywhere
  pickshell(vtmp, NDIM, 1.0);
  CROSSVP(yhat, xhat, vtmp);			// 2nd is orthogonal to 1st
  ymag = absv(yhat);
  DIVVS(yhat, yhat, ymag);			// and of unit magnitude
  CROSSVP(zhat, xhat, yhat);			// 3rd is orthogonal to both
  eprintf("[%s.pickrotation: dot products =  %f %f %f]\n", getprog(),
	  dotvp(xhat,yhat), dotvp(yhat,zhat), dotvp(zhat,xhat));
  for (int i = 0; i < NDIM; i++) {
    rotn[0][i] = xhat[i];
    rotn[1][i] = yhat[i];
    rotn[2][i] = zhat[i];
  }
}

//  fitroot: expand root cell to fit particle configuration
//  in smallest enclosing power-of-2 box.
//  _______________________________________________________
 
local void fitroot(treeptr tptr, bodyptr btab, int nbody)
{
  real dmax = 0.0;				// find max dist from root
  vector rp1;					// positions in tree frame
 
  for (bodyptr p = btab; p < btab+nbody; p++) {	// loop over all bodies
    SUBV(rp1, Pos(p), Trans(tptr));		// offset by tree origin
    MULVS(rp1, rp1, Scale(tptr));
    MULMV1(rp1, Rotn(tptr), rp1);		// rotate to tree frame
    for (int k = 0; k < NDIM; k++)		// loop over dimensions
      dmax = MAX(rabs(rp1[k]), dmax);		// find max dist from center
  }
  Rsize(tptr) = 1.0;				// set minimum root size
  while (Rsize(tptr) < 2 * dmax)		// loop until bodies fit
    Rsize(tptr) = 2 * Rsize(tptr);		// doubling box each time
  eprintf("[%s.fitroot: dmax = %f  rsize = %f]\n", getprog(),
	  dmax, Rsize(tptr));
}

//  loadbody: descend tree and insert body p in appropriate place.
//  NOTE: should check depth of tree and enforce limit!
//  ______________________________________________________________
 
local void loadbody(treeptr tptr, bodyptr p)
{
  cellptr q = Root(tptr), c;			// start q at root
  vector rp1, rb1;
  bodyptr b;
  real qsize = Rsize(tptr);			// track cell size
  int ipq, ibc;					// subcell indicies

  SUBV(rp1, Pos(p), Trans(tptr));		// offset by tree origin
  MULVS(rp1, rp1, Scale(tptr));			// scale by tree scale
  MULMV1(rp1, Rotn(tptr), rp1);			// rotate to tree frame
  assert(ABS(rp1[0]) <= qsize/2 &&
	 ABS(rp1[1]) <= qsize/2 &&
	 ABS(rp1[2]) <= qsize/2);		// check body fits in root
  Next(p) = NULL;				// terminate body list
  ipq = subindex(rp1, q);			// get index of p in root
  while (Desc(q)[ipq] != NULL) {		// loop descending tree
    if (Body(Desc(q)[ipq])) {			// found a body in subcell?
      b = (bodyptr) Desc(q)[ipq];
      SUBV(rb1, Pos(b), Trans(tptr));		// offset by tree origin
      MULVS(rb1, rb1, Scale(tptr));		// scale by tree scale
      MULMV1(rb1, Rotn(tptr), rb1);		// rotate to tree frame
      if (EQUALV(rp1, rb1)) {			// check for collision
	// do this test in tree CS because two bodies with slightly different
	// sim. positions may map to same tree CS position due to roundoff
	eprintf("[%s.loadbody: warning: bodies %p, %p have same position]\n",
		getprog(), (void *) p, (void *) b);
	Next(p) = (nodeptr) b;			// link old body into list
	break;					// done with search
      }
      Ncell(tptr)++;				// track number of cells
      c = makecell();				// allocate new cell
      for (int k = 0; k < NDIM; k++)		// set offset from parent
	GCPos1(c)[k] = (rp1[k] < GCPos1(q)[k] ? GCPos1(q)[k] - qsize/4 :
			                        GCPos1(q)[k] + qsize/4);
      ibc = subindex(rb1, c);			// index old body in new cell
      Desc(c)[ibc] = (nodeptr) b;		// put old body in new cell
      Desc(q)[ipq] = (nodeptr) c;		// link new cell in tree
    }
    q = (cellptr) Desc(q)[ipq];			// advance to next level
    assert(Cell(q));
    qsize = qsize / 2;				// shrink current cell
    ipq = subindex(rp1, q);			// get index for next level
  }
  Desc(q)[ipq] = (nodeptr) p;			// found place, store p
}

//  subindex: compute subcell index for position pos in cell q.
//  ___________________________________________________________
 
local int subindex(vector pos, cellptr q)
{
  int ind = 0;					// accumulate subcell index

  for (int k = 0; k < NDIM; k++)		// loop over dimensions
    if (GCPos1(q)[k] <= pos[k])			// if beyond midpoint
      ind += NDESC >> (k + 1);			// then skip over subcells
  return (ind);
}

//  threadtree: do a recursive treewalk starting from node p,
//  with next stop q, installing Next and More links.
//  _________________________________________________________
 
local void threadtree(nodeptr p, nodeptr q, int lev)
{
  nodeptr subn[NDESC+1];			// vector of subnodes
  int nsubn = 0;				// count nodes in vector
 
  if (lev >= MAXLEVEL)				// limit depth of tree
    fatal("%s.threadtree: tree depth exceeds maximum!\n", getprog());
  if (Body(p)) {				// reached body (list)
    while (Next(p) != NULL)			// find end of list
      p = Next(p);
    Next(p) = q;				// link into tree
  } else {					// if descendents to thread
    Next(p) = q;				// set link to next node
    for (int i = 0; i < NDESC; i++)		// loop over all subcells
      if (Desc(p)[i] != NULL)			// if this one is occupied
	subn[nsubn++] = Desc(p)[i];		// then store it in vector
    More(p) = subn[0];				// set link to 1st one
    subn[nsubn] = q;				// thread last one to next
    for (int i = 0; i < nsubn; i++)		// loop over descendents
      threadtree(subn[i], subn[i+1], lev+1);	// and thread them together
    Type(p) = CellType | ThreadFlag;		// re-initialize node type
    Level(p) = lev;				// store level within tree
  }
}

//  evalcells: evaluate cell properties (mass, center of mass position,
//  critical radius, etc) by recursive descent of threaded tree.
//  ___________________________________________________________________
 
local void evalcells(cellptr p, real psize, treeptr tptr)
{
  vector rgeo1, rcom1, dr1;			// positions in tree CS
  real phalf = psize / 2.0;
  vector dr, trQ;
  matrix Qp, Qq2, Qq3;

  Nlevel(tptr) = MAX(Nlevel(tptr), Level(p));	// remember maximum level
  SETV(rgeo1, GCPos1(p));			// save cell's geo. center
  Mass(p) = 0.0;				// init cell's total mass
  CLRV(Pos(p));					// zero CM (shared w/ GCPos1)
  for (nodeptr q = More(p); q != Next(p); q = Next(q)) {
						// loop over descendents
    if (Cell(q))				// recurse on cells
      evalcells((cellptr) q, phalf, tptr);
    Mass(p) += Mass(q);				// accumulate total mass
    ADDMULVS(Pos(p), Pos(q), Mass(q));		// and CM position
    Active(p) = MIN(Active(p) + Active(q),	// propagate active count
		    (1 << (8 * sizeof(Active(p)))) - 1);
  }
  if (Mass(p) > 0.0) {				// usually, cell has mass
    DIVVS(Pos(p), Pos(p), Mass(p));		// finalize CM position
  } else {					// but if no mass inside
    eprintf("[%s.evalcells: warning: cell %p has no mass]\n", getprog(), p);
    SETV(Pos(p), Pos(More(p)));			// use any node within cell
  }
  SUBV(rcom1, Pos(p), Trans(tptr));		// map CM to tree CS
  MULVS(rcom1, rcom1, Scale(tptr));
  MULMV1(rcom1, Rotn(tptr), rcom1);
  SUBV(dr1, rcom1, rgeo1);			// get offset in tree CS
  if (ABS(dr1[0])>phalf || ABS(dr1[1])>phalf || ABS(dr1[2])>phalf) {
    eprintf("[%s.evalcells: WARNING: center of mass is outside cell:\n"
	    "  lev = %d  phalf = %a  dr1 = (%a,%a,%a)\n"
	    "  rcom1 = (%a,%a,%a)\n  rgeo1 = (%a,%a,%a)]\n",
	    getprog(), Level(p), phalf, dr1[0], dr1[1], dr1[2],
	    rcom1[0], rcom1[1], rcom1[2], rgeo1[0], rgeo1[1], rgeo1[2]);
    if (scanopt(Options(tptr), "fatal-centmass"))
      fatal("%s.evalcells: center of mass is outside cell\n", getprog());
  }
  setrcrit(p, dr1, psize, tptr);		// crit. radius uses offset
  if (Usequad(tptr)) {
    CLRM(Qp);					// accummulate quad moment
    for (nodeptr q = More(p); q != Next(p); q = Next(q)) {
						// loop over descendents
      SUBV(dr, Pos(q), Pos(p));			// find displacement vect
      OUTVP(Qq2, dr, dr);			// form outer product
      MULMS(Qq2, Qq2, 3 * Mass(q));		// scale by 3*mass of node
      if (Cell(q)) {				// descendent is a cell?
	SETMSM(Qq3, Quad(q));			// get its moment
	ADDM(Qq2, Qq2, Qq3);			// and add to total
      }
      ADDM(Qp, Qp, Qq2);			// sum each quad moments
    }
    SETSMM(Quad(p), Qp);			// convert to symmetric form
  }
}

//  setrcrit: assign critical radius for cell p.
//  ____________________________________________

local void setrcrit(cellptr p, vector dr1, real psize, treeptr tptr)
{
  real rcrit;

  if (Theta(tptr) > 0.0) {			// do hierarchical algorithm
    if (! bh86)					// use default criterion
      rcrit = psize/Theta(tptr) + absv(dr1);	// add offset from center
    else
      rcrit = psize/Theta(tptr);		// use original criterion
  } else
    rcrit = 2 * Rsize(tptr);			// force all cells open
  rcrit = rcrit / Scale(tptr);			// convert to sim. coords
#ifndef MODTREECODE
  Rcrit2(p) = rcrit * rcrit;			// store square of radius
#else
  Rcrit(p) = rcrit;				// store radius itself
#endif
}

//  tracequads: make quad moments traceless, saving original trace.
//  _______________________________________________________________

local void tracequads(treeptr tptr)
{
  nodeptr p = (nodeptr) Root(tptr);
  matrix I, Qp, It3;
  real Qtr;

  SETMI(I);
  while (p != NULL) {
    if (Cell(p)) {
      SETMSM(Qp, Quad(p));
      TRACEM(Qtr, Qp);
      MULMS(It3, I, Qtr / 3.0);
      SUBM(Qp, Qp, It3);
      SETSMM(Quad(p), Qp);
      QuadEps2(p) = Qtr * Eps(tptr) * Eps(tptr) / 3.0;
      p = More(p);
    } else {
      p = Next(p);
    }
  }
}

//  makecell: return pointer to free oct-tree cell.
//  _______________________________________________
 
local cellptr makecell()
{
  cellptr c = (cellptr) allocate(sizeof(cell));	// allocate a new cell

  Type(c) = CellType | OctFlag;			// initialize node type
  Active(c) = 0;				// and force active flag
  for (int i = 0; i < NDESC; i++)		// loop over subcells
    Desc(c)[i] = NULL;				// and empty each one
  return (c);					// return pointer to cell
}
