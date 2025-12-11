/*
 * cputime.c: compute total process CPU time in minutes.
 */

#include "stdinc.h"
#include "getparam.h"

#include <sys/time.h>
#include <sys/resource.h>

double cputime(void)
{
  struct rusage buf;

  if (getrusage(RUSAGE_SELF, &buf) == -1)
    error("%s.cputime: getrusage() call failed\n", getprog());
  eprintf("[%s.cputime: utime.sec = %u  utime.usec = %u  stime.sec = %u"
	  "  stime.usec = %u  maxrss = %lu  ixrss = %lu  idrss = %lu"
	  "  isrss = %lu  minflt = %lu  majflt = %lu  nswap = %lu"
	  "  inblock = %lu  oublock = %lu  msgsnd = %lu  msgrcv = %lu"
	  "  nsignals = %lu  nvcsw = %lu  nivcsw = %lu]\n", getprog(),
	  buf.ru_utime.tv_sec, buf.ru_utime.tv_usec, buf.ru_stime.tv_sec,
	  buf.ru_stime.tv_usec, buf.ru_maxrss, buf.ru_ixrss, buf.ru_idrss,
	  buf.ru_isrss, buf.ru_minflt, buf.ru_majflt, buf.ru_nswap,
	  buf.ru_inblock, buf.ru_oublock, buf.ru_msgsnd, buf.ru_msgrcv,
	  buf.ru_nsignals, buf.ru_nvcsw, buf.ru_nivcsw);
  return ((buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1000000.0 +
	   buf.ru_stime.tv_sec + buf.ru_stime.tv_usec / 1000000.0) / 60.0);
}

#if defined(OBSOLETE_CODE)

#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>

double cputime(void)
{
    struct tms buffer;

    if (times(&buffer) == -1)
	error("%s.cputime: times() call failed\n", getprog());
    return ((buffer.tms_utime + buffer.tms_stime) / (60.0 * HZ));
}

#endif
