/* 
 * error.c: routines to report errors, etc.
 */

#include "stdinc.h"
#include <stdarg.h>
#include <string.h>

#define ERRORFILE  "fatal_error.log"		// default error file

local stream errstr = NULL;			// use to report errors

//  error: announce error and exit cleanly.
//  _______________________________________

void error(string fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);			// printf msg to std error
  va_end(ap);
  fflush(stderr);
  if (errstr != NULL) {				// error logging enabled?
    va_start(ap, fmt);
    vfprintf(errstr, fmt, ap);
    va_end(ap);
    fflush(errstr);
  }
  exit(1);					// quit with error status
}

//  fatal: announce error, write to file, and leave core image.
//  ___________________________________________________________

void fatal(string fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);			// printf msg to std error
  va_end(ap);
  fflush(stderr);
  if (errstr == NULL)				// default to error log
    errstr = fopen(ERRORFILE, "a");
  if (errstr != NULL) {				// have valid log file?
    va_start(ap, fmt);
    vfprintf(errstr, fmt, ap);			// printf msg to log file
    va_end(ap);
    fflush(errstr);
  }
  abort();					// quit, leave core image
}

//  set_error_stream: direct error message to given stream.
//  _______________________________________________________

void set_error_stream(stream str)
{
  errstr = str;
}

//  eprintf: print messages and warnings.  Uses "ZENO_MSG_OPTION" env. var.
//  to control printing; "warn" or "none" suppress some or all output.
//  _______________________________________________________________________


void eprintf(string fmt, ...)
{
  static string msgopt = NULL;
  va_list ap;

  if (msgopt == NULL)				// if NULL, get env. value
    msgopt = getenv("ZENO_MSG_OPTION");
  if ((msgopt == NULL) ||
      (strne(msgopt, "none") &&
       (strne(msgopt, "warn") || strcasestr(fmt, "warn") != NULL))) {
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);			// printf msg to std error
    va_end(ap);
    fflush(stderr);				// drain std error buffer
  }
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
  "logfile=",
  "fatal=false",
  NULL,
};

int main(int argc, string argv[])
{
  initparam(argv, defv);
  if (!strnull(getparam("logfile")))
    set_error_stream(fopen(getparam("logfile"), "w"));
  eprintf("[%s: foo=%f  bar=%d  fum=\"%s\"]\n",
	  getprog(), 2.7183, 16384, "bilbo");
  eprintf("[%s: WARNING: foo=%f  bar=%d  fum=\"%s\"]\n",
	  getprog(), 3.1415, 32768, "frodo");
  eprintf("[%s: warning: foo=%f  bar=%d  fum=\"%s\"]\n",
	  getprog(), 3.1415, 32768, "sam");
  eprintf("[%s: Warning: foo=%f  bar=%d  fum=\"%s\"]\n",
	  getprog(), 3.1415, 32768, "gollum");
  if (! getbparam("fatal"))
    error("error: foo=%f  bar=%d  fum=\"%s\"\n", 0.5772, 65536, "gimli");
  else
    fatal("fatal: foo=%f  bar=%d  fum=\"%s\"\n", 0.5772, 65536, "legolas");
  return (0);
}

#endif
