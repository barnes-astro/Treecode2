/*
 * treeio.c: I/O routines for hierarchical N-body code.
 * Copyright (c) 2026 by Joshua E. Barnes, Honolulu, Hawaii.
 */

#define global extern				// allocate globals elsewhere

#include "treecode.h"
#include <strings.h>
#include <sys/types.h>
#include <sys/stat.h>

//  Prototypes for local routines and variables.
//  ____________________________________________

local stream logstr;				// log stream, or stdout
local bool csv = FALSE;				// if TRUE, use csv format

local real rsize, cpubuild, cpugrav;		// data saved by forcereport
local int ncell, nlevel, ntrwalk;
local long nbbcalc, nbccalc;

local void diagnostics(void);			// eval N-body diagnostics:
local real mtot;		                // total mass of system
local real etot[3];				// Etot, KE, PE of system
local vector Rbulk;				// center of mass position
local vector Vbulk;				// center of mass velocity
local vector Jbulk;				// angular momentum vector
local vector Fbulk;				// bulk force on system
local vector Tbulk;				// bulk torque on system

//  inputdata: read initial conditions from input file
//  in csv format: mass, x, y, z, vx, vy, vz.
//  __________________________________________________

void inputdata(void)
{
  stream instr = stropen(infile, "r");		// open input stream
  int nxtchr, nline = 0;
  double m, x, y, z, vx, vy, vz;
  bodyptr p = bodytab; 
  
  do {						// loop reading line by line
    ungetc(nxtchr = fgetc(instr), instr);	// look ahead in input
    nline++;
    if (nxtchr != '#') {			// a line of data?
      if (fscanf(instr, "%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
		 &m, &x, &y, &z, &vx, &vy, &vz) != 7)
	error("%s.inputdata: read error on line %d\n", getprog(), nline);
      if (p >= bodytab + nbody)			// make sure there's room
	error("%s.inputdata: too much data in file\n", getprog());
      Mass(p) = m;				// store body data
      Pos(p)[0] = x;  Pos(p)[1] = y;  Pos(p)[2] = z;
      Vel(p)[0] = vx; Vel(p)[1] = vy; Vel(p)[2] = vz;
      p++;					// advance body pointer
    } else {					// a line of comments
      while (fgetc(instr) != '\n')		// skip to end of line
	if (feof(instr))			// not likely!
	  error("%s.inputdata: error skipping line %d\n", getprog(), nline);
    }
  } while (! feof(instr));			// until end of file
  if (p != bodytab + nbody)			// read exactly number
    error("%s.inputdata: unexpected end of file\n", getprog());
}

//  startoutput: begin output to log and state files.
//  _________________________________________________

void startoutput(string defv[])
{
  string ext, lead;

  ext = rindex(logfile, '.');
  if (ext != NULL && streq(ext, ".csv")) {
    eprintf("[%s: writing .csv log file]\n", getprog());
    csv = TRUE;
  }
  lead = csv ? "#" : "";
  logstr = strnull(logfile) ? NULL : stropen(logfile, "w!");
  if (logstr != NULL) {				// announce start of run
    fprintf(logstr, "%s%s [VERS. %s]: ", lead, getprog(), getversion());
    for (int i = 0; *(defv[i]) == ';'; i++)
      fprintf(logstr, "%s ", defv[i] + 1);
    fprintf(logstr, "\n%s\n%s %11s %11s %11s %11s %11s %11s\n", lead, lead,
	    "nbody", "dtime", "eps", "dtout", "tnow", "tstop");
    fprintf(logstr, "%s %11d %11.5f %11.4f %11.5f %11.4f %11.4f\n", lead,
	    nbody, dtime, eps, dtout, tnow, tstop);
#if !defined(MODTREECODE)
    fprintf(logstr, "%s\n%s %11s %11s %11s %11s %11s\n", lead, lead,
	    "theta", "quad", "rotate", "scale", "translate");
    fprintf(logstr, "%s %11.2f %11s %11s %11s %11.2f\n%s\n", lead,
	    theta, quad ? "true" : "false", rotate ? "true" : "false",
	    scale ? "true" : "false", translate, lead);
#else
    fprintf(logstr, "%s\n%s %11s %11s %11s %11s %11s %11s\n", lead, lead,
	    "theta", "quad", "rotate", "scale", "translate", "nshare");
    fprintf(logstr, "%s %11.2f %11s %11s %11s %11.2f %11d\n%s\n", lead,
	    theta, quad ? "true" : "false", rotate ? "true" : "false",
	    scale ? "true" : "false", translate, nshare, lead);
#endif
    if (! strnull(options))			// print options, if any
      fprintf(logstr, "%s\toptions: %s\n%s\n", options, lead, lead);
    if (csv)
      fprintf(logstr, "#%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s," \
	      "%s\n",	"time", "Rsize", "Ncell", "Nlevel", "CPUbuild",
	      "Ntrwalk", "Nbbcalc", "Nbccalc", "CPUgrav", "|T+U|", "-T/U",
	      "|Vbulk|", "|Jbulk|", "|Fbulk|", "|Tbulk|", "CPUtot", "File");

  }
  fflush(NULL);					// empty all output buffers
}

//  forcereport: store statistics on tree construction and force calculation.
//  _________________________________________________________________________

void forcereport(treeptr tptr)
{
  rsize = Rsize(tptr);
  ncell = Ncell(tptr);
  nlevel = Nlevel(tptr);
  cpubuild = CPUbuild(tptr);
  ntrwalk = Ntrwalk(tptr);
  nbbcalc = Nbbcalc(tptr);
  nbccalc = Nbccalc(tptr);
  cpugrav = CPUgrav(tptr);
}

//  output: output particle data and simulation diagnostics.  Particles
//  are written using 'a' format to insure full precision; to read data
//  output files in Python, do:
//      m,x,y,z,vx,vy,vz = loadtxt('snap0000.csv', delimiter=',',
//          converters=float.fromhex, unpack=True, encoding=None)
//  ___________________________________________________________________

#define SEP ","
#if !defined(DOUBLEPREC)
#  define FMT  "%.6a"
#else
#  define FMT  "%.13a"
#endif

void output(void)
{
  char namebuf[256] = { CNULL };
  stream outstr;

  if (! strnull(outfile) && ABS(tout-tnow) < dtime/2) {	// time for output?
    snprintf(namebuf, sizeof(namebuf), outfile, nstep);	// construct file name
    outstr = stropen(namebuf, "w");		// create & open for output
    fprintf(outstr, "# %s: time = %f  nstep = %d\n", getprog(), tnow, nstep);
    fprintf(outstr, "# mass,x,y,z,vx,vy,vz\n");
    for (bodyptr p = bodytab; p < bodytab+nbody; p++)
      fprintf(outstr, FMT SEP FMT SEP FMT SEP FMT SEP FMT SEP FMT SEP FMT
	      "\n", Mass(p), Pos(p)[0], Pos(p)[1], Pos(p)[2],
	      Vel(p)[0], Vel(p)[1], Vel(p)[2]);
    fclose(outstr);				// close up output file
    tout += dtout;				// schedule next output
  }
  if (logstr != NULL) {				// is log file open?
    diagnostics();				// compute std diagnostics
    if (! csv) {				// write formatted text log
      fprintf(logstr, "\n %7s %9s %7s %9s %9s %11s %11s %9s\n", "Rsize",
	      "Ncell", "Nlevel", "CPUbuild", "Ntrwalk", "Nbbcalc", "Nbccalc",
	      "CPUgrav");
      fprintf(logstr, " %7.0f %9d %7d %9.5f %9d %11ld %11ld %9.5f\n", rsize,
	      ncell, nlevel, cpubuild, ntrwalk, nbbcalc, nbccalc, cpugrav);
      fprintf(logstr, "\n%12s %8s %8s %9s %9s %9s %9s %8s\n", "time", "|T+U|",
	      "-T/U", "|Vbulk|", "|Jbulk|", "|Fbulk|", "|Tbulk|", "CPUtot");
      fprintf(logstr, "%12.5f %8.5f %8.5f %9.2e %9.2e %9.2e %9.2e %8.2f\n",
	      tnow, ABS(etot[0]), -etot[1]/etot[2], absv(Vbulk), absv(Jbulk),
	      absv(Fbulk), absv(Tbulk), cputime());
      if (! strnull(namebuf))
	fprintf(logstr, "\n\tdata output to file %s at time %f\n",
		namebuf, tnow);
    } else					// write csv log
      fprintf(logstr, "%.6e,%.6e,%d,%d,%.6e,%d,%ld,%ld,%.6e,%.6e,%.6e,"
	      "%.6e,%.6e,%.6e,%.6e,%.6e,%s\n", tnow, rsize, ncell, nlevel,
	      cpubuild, ntrwalk, nbbcalc, nbccalc, cpugrav, ABS(etot[0]),
	      -etot[1]/etot[2], absv(Vbulk), absv(Jbulk), absv(Fbulk),
	      absv(Tbulk), cputime(), namebuf);
  }
  fflush(NULL);					// empty all output buffers
}

//  diagnostics: compute set of dynamical diagnostics.
//  __________________________________________________

local void diagnostics(void)
{
  real velsq;
  vector tmpv;

  mtot = 0.0;					// zero total mass
  etot[1] = etot[2] = 0.0;			// zero total KE and PE
  CLRV(Rbulk);					// zero c. of m. position
  CLRV(Vbulk);					// zero c. of m. velocity
  CLRV(Jbulk);					// zero am vector
  CLRV(Fbulk);
  CLRV(Tbulk);
  for (bodyptr p = bodytab; p < bodytab+nbody; p++) {
    mtot += Mass(p);				// sum particle masses
    DOTVP(velsq, Vel(p), Vel(p));		// square vel vector
    etot[1] += 0.5 * Mass(p) * velsq;		// sum current KE
    etot[2] += 0.5 * Mass(p) * Phi(p);		// and PE, weighted right
    MULVS(tmpv, Vel(p), 0.5 * Mass(p));		// sum 0.5 m v_i v_j
    MULVS(tmpv, Pos(p), Mass(p));		// sum cm position
    ADDV(Rbulk, Rbulk, tmpv);
    MULVS(tmpv, Vel(p), Mass(p));		// sum cm momentum
    ADDV(Vbulk, Vbulk, tmpv);
    CROSSVP(tmpv, Vel(p), Pos(p));		// sum angular momentum
    MULVS(tmpv, tmpv, Mass(p));
    ADDV(Jbulk, Jbulk, tmpv);
    MULVS(tmpv, Acc(p), Mass(p));		// sum bulk force
    ADDV(Fbulk, Fbulk, tmpv);
    CROSSVP(tmpv, Acc(p), Pos(p));		// sum bulk torque
    MULVS(tmpv, tmpv, Mass(p));
    ADDV(Tbulk, Tbulk, tmpv);
  }
  etot[0] = etot[1] + etot[2];			// sum KE and PE
  DIVVS(Rbulk, Rbulk, mtot);        		// normalize cm coords
  DIVVS(Vbulk, Vbulk, mtot);
  CROSSVP(tmpv, Vbulk, Rbulk);			// correct for CM offset
  SUBV(Jbulk, Jbulk, tmpv);
  CROSSVP(tmpv, Fbulk, Rbulk);
  SUBV(Tbulk, Tbulk, tmpv);
}
