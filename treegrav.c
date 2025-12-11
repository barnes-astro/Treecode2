/*
 * treegrav.c: baseline hierarchical gravity calculation.
 * Copyright (c) 2025 by Joshua E. Barnes, Honolulu, Hawaii.
 */
 
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include <assert.h>
#include "treedefs.h"
#include "treewalk.c"

//  Local routines to perform force calculations.
//  _____________________________________________

local bool treewalkQ(bodyptr p, treeptr tptr);	// code in treewalk.c
local bool treewalkM(bodyptr p, treeptr tptr);	// code in treewalk.c

//  treegrav: compute gravitational forces on bodies in tree.
//  _________________________________________________________

void treegrav(treeptr tptr)
{
  double cpustart = cputime();			// record time at start
  bool selfskip;
  nodeptr q = (nodeptr) Root(tptr);		// start at root of tree

  Ntrwalk(tptr) = 0;				// initialize counters
  Nbbcalc(tptr) = 0;
  Nbccalc(tptr) = 0;
  while (q != NULL)				// loop scanning tree
    if (Body(q)) {				// reached a body...
      if (Active(q) != 0) {			// does it need updating?
	selfskip = (Usequad(tptr) ? treewalkQ((bodyptr) q, tptr)
		                  : treewalkM((bodyptr) q, tptr));
						// compute force on body
	if (! selfskip) {			// fail to skip self-int?
	  if (! scanopt(Options(tptr), "warn-selfskip"))
	    error("%s.treegrav: skipped self-interaction  pos = %f,%f,%f\n",
		  getprog(), Pos(q)[0], Pos(q)[1], Pos(q)[2]);
	  eprintf("[%s.treegrav: WARNING: skipped self-interaction  pos = "
		  "%f,%f,%f]\n",getprog(), Pos(q)[0], Pos(q)[1], Pos(q)[2]);
	}
      }
      q = Next(q);				// move on to next node
    } else {					// reached a cell...
      if (Active(q) != 0)			// active bodies in cell?
	q = More(q);				// yes, so scan descendents
      else
	q = Next(q);				// no, done with this cell
    }
  CPUgrav(tptr) = cputime() - cpustart;		// store elapsed CPU time
}

//  treegrav1: compute gravitational force on one body.
//  ___________________________________________________

void treegrav1(bodyptr p, treeptr tptr, bool selfcheck)
{
  double cpustart = cputime();			// record time at start
  bool selfskip;

  selfskip = (Usequad(tptr) ? treewalkQ(p, tptr) : treewalkM(p, tptr));
  if (selfcheck && ! selfskip) {
    if (! scanopt(Options(tptr), "warn-selfskip"))
      error("%s.treegrav1: skipped self-interaction  pos = %f,%f,%f\n",
	    getprog(), Pos(p)[0], Pos(p)[1], Pos(p)[2]);
    eprintf("[%s.treegrav1: WARNING: missed self-interaction  pos ="
	    " %f,%f,%f]\n", getprog(), Pos(p)[0], Pos(p)[1], Pos(p)[2]);
  }
  CPUgrav(tptr) += cputime() - cpustart;	// add elapsed CPU time
}
