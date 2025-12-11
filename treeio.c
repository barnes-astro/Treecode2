/*
 * treeio.c: I/O routines for hierarchical N-body code.
 * Copyright (c) 2025 by Joshua E. Barnes, Honolulu, Hawaii.
 */

#define global extern				// allocate globals elsewhere

#include "treecode.h"
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

//  inputdata: read initial conditions from input file.
//  ___________________________________________________

void inputdata(void)
{
  error("%s.inputdata: not yet implemented!\n", getprog());
#if 0
  stream instr = stropen(infile, "r");		// open input stream
  string intags[MaxBodyFields];

  get_history(instr);				// read file history
  bodytab = NULL;				// force body allocation
  if (! get_snap(instr, &bodytab, &nbody, &tnow, intags, FALSE, NULL))
    error("%s.inputdata: no data in input file\n", getprog());
  strclose(instr);				// close input stream
  if (! set_member(intags, MassTag))
    error("%s.inputdata: %s data missing\n", getprog(), MassTag);
  if (! set_member(intags, PosTag))
    error("%s.inputdata: %s data missing\n", getprog(), PosTag);
  if (! set_member(intags, VelTag))
    error("%s.inputdata: %s data missing\n", getprog(), VelTag);
  if (! set_member(intags, TypeTag)) {		// initialize body type
    eprintf("[%s.inputdata: setting body type = 0x%x]\n",
	    getprog(), BodyType);
    for (bodyptr p = bodytab; p < bodytab+nbody; p++)
      Type(p) = BodyType;
  } else {					// check types provided
    for (bodyptr p = bodytab; p < bodytab+nbody; p++)
      if (Cell(p) || ! Body(p))
	error("%s.inputdata: bad body type 0x%x\n", getprog(), Type(p));
  }
#endif
}

//  startoutput: begin output to log and state files.
//  _________________________________________________

void startoutput(string defv[])
{
  string lead;

  if (strlen(logfile) > 4 && streq(logfile + strlen(logfile) - 4, ".csv")) {
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
    fprintf(logstr, "%s %11.2f %11s %11s %11s %11.2f\n", lead,
	    theta, quad ? "true" : "false", rotate ? "true" : "false",
	    scale ? "true" : "false", translate);
#else
    fprintf(logstr, "%s\n%s %11s %11s %11s %11s %11s %11s\n", lead, lead,
	    "theta", "quad", "rotate", "scale", "translate", "nshare");
    fprintf(logstr, "%s %11.2f %11s %11s %11s %11.2f %11d\n", lead,
	    theta, quad ? "true" : "false", rotate ? "true" : "false",
	    scale ? "true" : "false", translate, nshare);
#endif
    if (! strnull(options))			// print options, if any
      fprintf(logstr, "%s\n%s\toptions: %s\n", lead, lead, options);
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

//  output: compute diagnostics and output binary data.
//  ___________________________________________________

void output(void)
{
  static bool csvhead = TRUE;
  char namebuf[256] = { CNULL };
  struct stat buf;
  stream outstr;

  diagnostics();				// compute std diagnostics
  if (!strnull(outfile) && ABS(tout-tnow) < dtime/2)	// time for output?
    snprintf(namebuf, sizeof(namebuf), outfile, nstep);	// construct file name
  if (logstr != NULL) {				// is log file open?
    if (! csv) {				// write human-readable log
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
    } else {					// write csv log
      if (csvhead) 
	fprintf(logstr, "#%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s," \
		"%s\n",	"time", "Rsize", "Ncell", "Nlevel", "CPUbuild",
		"Ntrwalk", "Nbbcalc", "Nbccalc", "CPUgrav", "|T+U|", "-T/U",
		"|Vbulk|", "|Jbulk|", "|Fbulk|", "|Tbulk|", "CPUtot", "File");
      csvhead = FALSE;
      fprintf(logstr, "%.6e,%.6e,%d,%d,%.6e,%d,%ld,%ld,%.6e,%.6e,%.6e,"
	      "%.6e,%.6e,%.6e,%.6e,%.6e,%s\n", tnow, rsize, ncell, nlevel,
	      cpubuild, ntrwalk, nbbcalc, nbccalc, cpugrav, ABS(etot[0]),
	      -etot[1]/etot[2], absv(Vbulk), absv(Jbulk), absv(Fbulk),
	      absv(Tbulk), cputime(), namebuf);
    }
  }
#if 0
  if (!strnull(outfile) && ABS(tout-tnow) < dtime/2) {	// time for output?
    if (stat(namebuf, &buf) != 0) {		// no output file exists?
      outstr = stropen(namebuf, "w");		// create & open for output
      put_history(outstr);			// write out history data
    } else					// else file already exists
      outstr = stropen(namebuf, "a");		// reopen and append output
    ntag = 0;					// list fields to output
    for (int k = 0; alltags[k] != NULL; k++)
      if (scanopt(outputs, alltags[k]) || (nstep == 0 && k < 2))
	outtags[ntag++] = alltags[k];		// add Type, Mass 1st time
    outtags[ntag] = NULL;
    put_snap(outstr, &bodytab, &nbody, &tnow, outtags);
    strclose(outstr);				// close up output file
    if (logstr != NULL && !csv)
      fprintf(logstr, "\n\tdata output to file %s at time %f\n",
	      namebuf, tnow);
    tout += dtout;				// schedule next output
  }
#endif
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
