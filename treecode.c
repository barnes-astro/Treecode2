/*
 * treecode.c: hierarchical N-body simulation code.
 * Copyright (c) 2025 by Joshua E. Barnes, Honolulu, Hawaii.
 */

#define global					// allocate global variables
#include "treecode.h"

//  Input parameters and default values.
//  ____________________________________

string defv[] = {
#if !defined(MODTREECODE)
				";Hierarchical N-body code:",
#else
				";Modified hierarchical N-body code:",
#endif
#if !defined(NOSOFTCORR)
                                ";forces corrected for softening.",
#else
                                ";no softening correction.",
#endif
  "in=",			";Initial conditions: input file name;",
				";if null, make Plummer model test data.",
  "out=",			";Simulation output: pattern for file name;",
				";directive (eg, %04x) formats step number.",
  "dtime=1/128",		";Timestep for leapfrog integration;",
				";if zero, does force calculation and exits.",
  "eps=0.01",			";Smoothing length for force calculation",
  "theta=0.8",			";Hierarchical force accuracy parameter",
  "quad=true",			";Enable use of quadupole moments",
  "rotate=true",		";Enable random tree orientation",
  "scale=true",			";Enable random tree scaling",
  "translate=1.0",		";Enable random tree translation;",
				";value is max translation distance.",
#if defined(MODTREECODE)
  "nshare=64",			";Max number of bodies sharing interactions",
#endif
  "options=",			";Comma-separated list of program options;",
				";bh86: original opening criteria,",
				";warn-selfint: allow self-interaction,",
				";new-tout: reschedule output times,",
  "dtout=1/4",			";Time interval between data outputs",
  "tstop=2.0",			";Time to stop simulation",
  "seed=123",			";Random seed for averaging and test data",
  "nbody=65536",		";Number of bodies for test data",
  "log=-",			";Log file name (default is stdout)",
  "cpumax=",			";If given, sets max CPU time in minutes",
  "VERSION=2.1exp",		";Joshua Barnes  24 September 2025",
  NULL,
};

//  Prototypes for local procedures.
//  ________________________________

local void startrun(void);			// initialize system state
local void newrun(void);			// start new simulation
local void testdata(void);			// generate test data
local void stepsystem(void);			// advance by one time-step
local void gravforce(void);			// do force calculation

//  main: toplevel routine for hierarchical N-body code.
//  ____________________________________________________

int main(int argc, string argv[])
{
  initparam(argv, defv);			// initialize param access
  startrun();					// get params & input data
  startoutput(defv);				// activate output code
  if (nstep == 0) {				// if data just initialized
    gravforce();				// calculate initial forces
    output();					// generate initial output
  }
  while (dtime > 0.0 && tnow < tstop &&		// evolve until terminated
	 (strnull(getparam("cpumax")) || cputime() < getdparam("cpumax"))) {
    stepsystem();				// advance step by step
    output();					// output results each time
  }
  return (0);					// exit with proper status
}

//  startrun: startup hierarchical N-body code.
//  ___________________________________________

local void startrun(void)
{
  eprintf("[%s: sizeof(node) = %d  sizeof(body) = %d  sizeof(cell) = %d]\n",
	  getprog(), sizeof(node), sizeof(body), sizeof(cell));
  infile = getparam("in");			// set I/O file names
  outfile = getparam("out");
  logfile = getparam("log");
  if (! checkfmt(outfile, "%0?x"))
    error("%s: bad directive in out=%s\n", getprog(), outfile);
  newrun();					// then start up new run
}

//  newrun: initialize new run from file or test data.
//  __________________________________________________

local void newrun(void)
{
  eps = getdparam("eps");			// set input parameters
  dtime = getdparam("dtime");
  theta = getdparam("theta");
  quad = getbparam("quad");
  rotate = getbparam("rotate");
  scale = getbparam("scale");
  translate = getdparam("translate");
  nshare = 0;					// default if not modified
#if defined(MODTREECODE)
  nshare = getiparam("nshare");			// set for modified code
#endif
  tstop = getdparam("tstop");
  dtout = getdparam("dtout");
  options = getparam("options");
  init_random(getiparam("seed"));		// set random number gen.
  if (! strnull(infile))			// if data file was given
    inputdata();				// then read inital data
  else {					// else make initial data
    nbody = getiparam("nbody");			// get number of bodies
    if (nbody < 1)				// check for silly values
      error("%s: absurd value for nbody\n", getargv0());
    testdata();					// and make plummer model
  }
  nstep = 0;					// begin counting steps
  tout = tnow;					// schedule first output
}

//  testdata: generate Plummer model initial conditions for test runs,
//  scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
//  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
//  __________________________________________________________________

local void testdata(void)
{
  real rsc, vsc, r, v, x, y, MCUT = 0.999;

  bodytab = (bodyptr) allocate(nbody * sizeof(body));
						// alloc space for bodies
  rsc = (3 * M_PI) / 16;			// set length scale factor
  vsc = rsqrt(1.0 / rsc);			// find speed scale factor
  for (bodyptr p = bodytab; p < bodytab+nbody; p++) {	// loop over bodies
    Type(p) = BodyType;				// tag as a body
    Mass(p) = 1.0 / nbody;			// set masses equal
    x = xrandom(0.0, MCUT);			// pick value for enc. mass
    r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1);	// find corresponding radius
    pickshell(Pos(p), NDIM, rsc * r);		// pick position vector
    do {					// select from fn g(x)
      x = xrandom(0.0, 1.0);			// for x in range 0:1
      y = xrandom(0.0, 0.1);			// max of g(x) is 0.092
    } while (y > x*x * rpow(1 - x*x, 3.5));	// using von Neumann tech
    v = x * rsqrt(2.0 / rsqrt(1 + r*r));	// find resulting speed
    pickshell(Vel(p), NDIM, vsc * v);		// pick velocity vector
  }
  tnow = 0.0;					// set elapsed model time
}

//  stepsystem: advance N-body system using simple leap-frog.
//  _________________________________________________________

local void stepsystem(void)
{
  for (bodyptr p = bodytab; p < bodytab + nbody; p++) {
    ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);	// advance v by 1/2 step
    ADDMULVS(Pos(p), Vel(p), dtime);		// advance r by 1 step
  }
  gravforce();					// compute forces
  for (bodyptr p = bodytab; p < bodytab + nbody; p++) {
    ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);	// advance v by 1/2 step
  }
  nstep++;					// count another time step
  tnow = tnow + dtime;				// finally, advance time
}

//  gravforce: supervise force calculation.
//  _______________________________________

local void gravforce(void)
{
  treeptr tptr;

  for (bodyptr p = bodytab; p < bodytab + nbody; p++)
    Active(p) = 1;				// update forces on all bodies
  tptr = treebuild(bodytab, nbody, eps, theta, quad,
		   rotate, scale, translate, options);
#ifndef MODTREECODE
  treegrav(tptr);				// compute current forces
#else
  treegrav(tptr, nshare);			// compute current forces
#endif
#if defined(TESTBODY)
  for (bodyptr p = bodytab; p < bodytab+nbody; p++)
    if (Mass(p) == 0.0)				// select bodies not in tree
      treegrav1(p, tptr, FALSE);		// compute forces individually
#endif
  forcereport(tptr);				// print force statistics
  treedestroy(tptr);

}
