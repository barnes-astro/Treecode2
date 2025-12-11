/*
 * treegrav_mod.c: hierarchical gravity calculation; implements a
 * variant of the "Modified Tree Code" (Barnes 1990).
 * Assumes that forces are required on ALL bodies.  
 * Copyright (c) 2025 by Joshua E. Barnes, Honolulu, Hawaii.
 */
 
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include <assert.h>
#include "treedefs.h"
#include "treewalk.c"

//  Local routines and variables to perform force calculations.
//  ___________________________________________________________

local void listintern(cellptr p);
local bool listextern(cellptr p, treeptr tptr);
local bool accept(nodeptr q, vector pmid, real rmax);
local void sumgrav(cellptr p, treeptr tptr);
local void sumgravB(bodyptr bi, bodyptr barr, int nb, real eps2);
local void sumgravQ(bodyptr bi, cellptr carr, int nc, real eps2);
local void sumgravC(bodyptr bi, cellptr carr, int nc, real eps2);
local bool treewalkQ(bodyptr p, treeptr tptr);	// code in treewalk.c
local bool treewalkM(bodyptr p, treeptr tptr);	// code in treewalk.c

local void *objarr;			// storage for bodies and cells
local int objlen = 0;			// length in bytes of objarr
local bodyptr bptr;			// pointer++ to bodies in objarr
local cellptr cptr;			// --pointer to cells in objarr
local bodyptr *intptr;			// pointers to internal bodies
local int nbint;			// number of internal bodies
local int nbext;			// number of external bodies
local int ncext;			// number of external cells
local int maxbody, maxcell;		// track source mem used
local int nsingle;			// count descents for single body
local int narrovf;			// count array overflows

//  Max. length of interaction list.  This should be a parameter!
//  _____________________________________________________________

#if !defined(MAXOBJ)
#define MAXOBJ  20000
#endif

//  treegrav: compute gravitational forces on all bodies in tree.
//  _____________________________________________________________

void treegrav(treeptr tptr, int nshare)
{
  double cpustart;
  nodeptr p;

  if (Active(Root(tptr)) <= nshare)		// fatal usr error (unlikely)
    error("%s.treegrav: nshare too large (Active(root) = %d)\n",
	  getprog(), Active(Root(tptr)));
  if (objlen == 0) {				// first time called?
    objlen = MAXOBJ * sizeof(cell);		// just a guess, for now
    objarr = allocate(objlen);			// alloc. mem. for sources
    intptr = (bodyptr *) allocate(nshare * sizeof(bodyptr));
  }
  cpustart = cputime();				// record time w/o alloc
  Ntrwalk(tptr) = 0;				// init various counters
  Nbbcalc(tptr) = 0;
  Nbccalc(tptr) = 0;
  maxbody = maxcell = 0;			// track max bodies and cells
  nsingle = narrovf = 0;			// count tree descents
  p = (nodeptr) Root(tptr);			// start at root of tree
  while (p != NULL)				// loop scanning tree
    if (Cell(p)) {				// reached a cell... 
      if (Active(p) <= nshare) {		// with nshare or less bodies
	listintern((cellptr) p);		// list bodies needing forces
	if (listextern((cellptr) p, tptr)) {	// list nodes generating field
	  sumgrav((cellptr) p, tptr);		// perform force summations
	  Ntrwalk(tptr)++;			// count a tree walk
	  maxbody = MAX(maxbody, nbext+nbint);	// track max counts
	  maxcell = MAX(maxcell, ncext);
	  p = Next(p);				// move on to next node
	} else
	  p = More(p);				// overflow; examine desc.
      } else					// cell has too many bodies
	p = More(p);				// so examine descendents
    } else {					// reached a single body...
      if (Usequad(tptr))			// quad terms requested?
	assert(treewalkQ((bodyptr) p, tptr));	// use quad moments
      else
	assert(treewalkM((bodyptr) p, tptr));	// use monopoles only
      Ntrwalk(tptr)++;				// count one tree walk
      nsingle++;				// for a single body
      p = Next(p);				// on to next node
    }
  CPUgrav(tptr) = cputime() - cpustart;		// store elapsed CPU time
  eprintf("[%s.treegrav: maxbody = %d  maxcell = %d"
	  "  nsingle = %d  narrovf = %d]\n",
	  getprog(), maxbody, maxcell, nsingle, narrovf);
  if (narrovf > 0)				// warn about overflows
    eprintf("[%s.treegrav: warning: %d array overflows (%d treewalks)]\n",
	    getprog(), narrovf, Ntrwalk(tptr));
}

//  treegrav1: compute gravitational force on one body.
//  ___________________________________________________

void treegrav1(bodyptr p, treeptr tptr, bool selfcheck)
{
  double cpustart = cputime();			// record time at start
  bool selfskip;

  selfskip = (Usequad(tptr) ? treewalkQ(p, tptr) : treewalkM(p, tptr));
  if (selfcheck && ! selfskip) {
    if (! scanopt(Options(tptr), "warn-selfint"))
      error("%s.treegrav1: skipped self-interaction  pos = %f,%f,%f\n",
	    getprog(), Pos(p)[0], Pos(p)[1], Pos(p)[2]);
    eprintf("[%s.treegrav1: WARNING: missed self-interaction  pos ="
	    " %f,%f,%f]\n", getprog(), Pos(p)[0], Pos(p)[1], Pos(p)[2]);
  }
  CPUgrav(tptr) += cputime() - cpustart;	// add elapsed CPU time
}

//  listintern: list all bodies internal to cell p.
//  _______________________________________________

local void listintern(cellptr p)
{
  nodeptr q;

  bptr = (bodyptr) objarr;			// start at front of array
  nbint = 0;					// count internal bodies
  q = More(p);					// begin with 1st subnode
  while (q != Next(p)) {			// loop over p's descendents
    if (Body(q)) {				// found a body within cell p
      Mass(bptr) = Mass(q);			// copy data for grav field
      SETV(Pos(bptr), Pos(q));			// into source array
      bptr++;					// advance insertion ptr
      intptr[nbint++] = (bodyptr) q;		// save address for update
      q = Next(q);				// and on to next node
    } else					// found a cell...
      q = More(q);				// examine its descendents
  }
  assert(nbint == Active(p));			// check all were listed
}

//  listextern: list nodes outside cell p with which bodies inside can
//  safely interact.  Returns FALSE if objarr overflows.
//  __________________________________________________________________

local bool listextern(cellptr p, treeptr tptr)
{
  bodyptr bint = (bodyptr) objarr;
  vector pmin, pmax, pmid;
  real prad = 0.0;
  nodeptr q;
  bool skipself = FALSE;

  SETV(pmin, Pos(bint));			// estimate internal midpoint 
  SETV(pmax, Pos(bint));			// start with first body
  for (bodyptr s = bint+1; s < bint+nbint; s++)	// loop over remaining bodies
    for (int k = 0; k < NDIM; k++) {		// and over coordinate index
      pmin[k] = MIN(pmin[k], Pos(s)[k]);	// find minimum of coordinates
      pmax[k] = MAX(pmax[k], Pos(s)[k]);	// and maximum as well
    }
  ADDV(pmid, pmin, pmax);			// average to get midpoint
  DIVVS(pmid, pmid, 2.0);
  for (bodyptr s = bint; s < nbint+bint; s++)	// loop over internal bodies
    prad = MAX(prad, distv(Pos(s), pmid));	// find bounding radius 
  cptr = (cellptr) (objarr + objlen);		// start at end of cell array
  q = (nodeptr) Root(tptr);			// scan from root of tree
  while (q != NULL) {				// loop scanning tree
    if (Body(q)) {				// found a body?
      Mass(bptr) = Mass(q);			// copy data for grav field
      SETV(Pos(bptr), Pos(q));			// into source array
      bptr++;					// advance insertion ptr
      q = Next(q);				// on to next node
    } else if (q == (nodeptr) p) {		// reached cell p itself?
      skipself = TRUE;				// note self-int skipped
      q = Next(q);				// on to next node
    } else if (accept(q, pmid, prad)) {		// found acceptable cell?
      cptr--;					// point to open slot
      Mass(cptr) = Mass(q);			// copy data for grav field
      SETV(Pos(cptr), Pos(q));			// into source array
      SETSM(Quad(cptr), Quad(q));		// including quad field
      QuadEps2(cptr) = QuadEps2(q);
      q = Next(q);				// on to next node
    } else {					// cell too close to accept
      q = More(q);				// examine descendents
    }
    if ((void *) bptr >= (void *) cptr) {	// check for over-write
      narrovf++;				// count array overflows
      return (FALSE);				// indicate failure
    }
  }
  assert(skipself);				// check self-int skipped
  nbext = (bptr - bint) - nbint;		// count external bodies
  ncext = (cellptr) (objarr + objlen) - cptr;	// and external cells
  return (TRUE);				// indicate success
}

//  accept: determine if bodies in p, all within prad of pmid, can
//  safely interact with cell q.
//  ______________________________________________________________

local bool accept(nodeptr q, vector pmid, real prad)
{
  real dqp2;

  DISTSQV(dqp2, Pos(q), pmid);
  return (dqp2 > (Rcrit(q) + prad) * (Rcrit(q) + prad));
}

//  sumgrav: compute gravity on bodies in cell p due to nodes in objarr.
//  ____________________________________________________________________

local void sumgrav(cellptr p, treeptr tptr)
{
  bodyptr bint = (bodyptr) objarr, bext = ((bodyptr) objarr) + nbint, bi, bj;
  real eps2 = Eps(tptr) * Eps(tptr), dr2, dr2i, dr1i, phiM, mdr3i;
  vector dr, acc;

  for (bi = bint; bi < bint+nbint; bi++) {	// loop over int. body copies
    Phi(bi) = 0.0;				// initialize phi, acc
    CLRV(Acc(bi));
  }
  for (bi = bint+1; bi < bint+nbint; bi++)	// do O(N^2) force calculation
    for (bj = bint; bj < bi; bj++) {		// for set of internal bodies
      DOTPSUBV(dr2, dr, Pos(bj), Pos(bi));
      dr2i = ((real) 1.0) / (dr2 + eps2);
      dr1i = rsqrt(dr2i);
      phiM = Mass(bj) * dr1i;
      mdr3i = phiM * dr2i;
      MULVS(acc, dr, mdr3i);
      Phi(bi) -= phiM;				// sum potentials
      Phi(bj) -= phiM;
      ADDV(Acc(bi), Acc(bi), acc);		// sum accelerations
      SUBV(Acc(bj), Acc(bj), acc);
    }
  Nbbcalc(tptr) += nbint * (nbint - 1);		// count interactions per body
  for (bi = bint; bi < bint+nbint; bi++) {
    sumgravB(bi, bext, nbext, eps2);		// do body-body interactions
    if (Usequad(tptr))
      sumgravQ(bi, cptr, ncext, eps2);		// do body-cell interactions
    else
      sumgravC(bi, cptr, ncext, eps2);
  }
  Nbbcalc(tptr) += nbint * nbext;		// count external interactions
  Nbccalc(tptr) += nbint * ncext;
  for (int i = 0; i < nbint; i++) {
    Phi(intptr[i]) = Phi(bint + i);		// copy results back to
    SETV(Acc(intptr[i]), Acc(bint + i));	// original bodies
  }
}

//  sumgravB: sum gravity on body bi due to bodies in barr[nb].
//  ___________________________________________________________

local void sumgravB(bodyptr bi, bodyptr barr, int nb, real eps2)
{
  real dr2, dr2i, dr1i, phiM, mdr3i;
  vector dr;

  for (bodyptr bj = barr; bj < barr + nb; bj++) {
    DOTPSUBV(dr2, dr, Pos(bj), Pos(bi));
    dr2i = ((real) 1.0) / (dr2 + eps2);
    dr1i = rsqrt(dr2i);
    phiM = Mass(bj) * dr1i;
    mdr3i = phiM * dr2i;
    Phi(bi) -= phiM;
    ADDMULVS(Acc(bi), dr, mdr3i);
  }
}

//  sumgravQ: sum quad-order gravity on body bi due to cells in carr[nc].
//  _____________________________________________________________________

local void sumgravQ(bodyptr bi, cellptr carr, int nc, real eps2)
{
  real dr2, dr2i, dr1i, phiM, mdr3i, qdr2, dr5i, phiQ;
  vector dr, qdr;

  for (cellptr cj = carr; cj < carr + nc; cj++) {
    DOTPSUBV(dr2, dr, Pos(cj), Pos(bi));
    dr2i = ((real) 1.0) / (dr2 + eps2);
    dr1i = rsqrt(dr2i);
    phiM = Mass(cj) * dr1i;
    mdr3i = phiM * dr2i;
    DOTPMULSMV(qdr2, qdr, Quad(cj), dr);
#if !defined(NOSOFTCORR)
    qdr2 -= QuadEps2(cj);			// apply Keigo's correction
#endif
    dr5i = dr2i * dr2i * dr1i;			// form inverse 5th power
    phiQ = ((real) 0.5) * dr5i * qdr2;		// form quad potential
    Phi(bi) -= phiM + phiQ;			// add mono and quad pot
    mdr3i += ((real) 5.0) * phiQ * dr2i;	// adjust radial term
    ADDMULVS2(Acc(bi), dr, mdr3i, qdr, - dr5i);	// mono and quad acc
  }
}

//  sumgravC: sum gravity on body bi due to cells in carr[nc].
//  __________________________________________________________

local void sumgravC(bodyptr bi, cellptr carr, int nc, real eps2)
{
  real dr2, dr2i, dr1i, phiM, mdr3i;
  vector dr;

  for (cellptr cj = carr; cj < carr + nc; cj++) {
    DOTPSUBV(dr2, dr, Pos(cj), Pos(bi));
    dr2i = ((real) 1.0) / (dr2 + eps2);
    dr1i = rsqrt(dr2i);
    phiM = Mass(cj) * dr1i;
    mdr3i = phiM * dr2i;
    Phi(bi) -= phiM;
    ADDMULVS(Acc(bi), dr, mdr3i);
  }
}
