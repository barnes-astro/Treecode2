/*
 * treecode.h: define various things for treecode.c and treeio.c.
 * Copyright (c) 2025 by Joshua E. Barnes, Honolulu, Hawaii.
 */

#ifndef _treecode_h
#define _treecode_h

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "treedefs.h"

//  Parameters, state variables, and diagnostics for N-body integration.
//  ____________________________________________________________________

global string infile;			// file name for snapshot input
global string outfile;			// file pattern for snapshot output
global string logfile;			// file name for log file
global string savefile;			// file pattern for state output
global string restorefile;		// file name for state input
global string streamfile;		// file pattern for frame stream
global real dtime;			// basic integration timestep
global real eps;			// force softening length
global real theta;			// force accuracy parameter
global bool quad;			// flag for quadrupole moments
global bool rotate;			// flag for random tree orientation
global bool scale;			// flag for random tree scaling
global real translate;			// max distance to shift tree origin
global int nshare;			// max bodies sharing interactions
global string options;			// comma-separated keywords
global real dtout;			// data output timestep
global string outputs;			// list of data to output
global real tstop;			// time to stop calculation
global real tnow;			// current value of time
global real tout;			// time of next output
global int nstep;			// number of time-steps
global int nbody;			// number of bodies in system
global bodyptr bodytab;			// points to array of bodies

//  Prototypes for I/O routines.
//  ____________________________

void inputdata(void);			// read initial data file
void startoutput(string defv[]);	// open files for output
void forcereport(treeptr tptr);		// report on force calculation
void output(void);			// perform output operation
void savestate(void);			// save system state		
void restorestate(void);		// restore system state

#endif // ! _treecode_h
