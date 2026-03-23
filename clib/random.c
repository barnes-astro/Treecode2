/*
 * random.c: useful functions for random numbers.  This version uses the
 * "better" generator random(3) from libc, for applications w/o GSL.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include <string.h>

local unsigned long rng_min, rng_max;

//  init_random: initialize random number generator.
//  ________________________________________________

void init_random(unsigned long seed)
{
  rng_min = 0;
  rng_max = 0x7fffffff;
  eprintf("[%s.init_random: generator random(3), range %lu:%lu]\n",
	  getprog(), rng_min, rng_max);
  srandom((unsigned) seed);
}

//  xrandom: floating-point random number in range [xl,xh], inclusive;
//  computed by scaling unsigned long, so at most 32 bits are random.
//  __________________________________________________________________

double xrandom(double xl, double xh)
{
  if (rng_min == rng_max) {			// uninitialized?
    eprintf("[%s.xrandom: WARNING: using default seed]\n", getprog());
    init_random(0);
  }
  return (xl + (xh - xl) *
	  (random() - rng_min) / ((double) (rng_max - rng_min)));
}

//  grandom: normally distributed random number (polar method).
//  Reference: Knuth, vol. 2, p. 104.
//  ___________________________________________________________

double grandom(double mean, double sdev)
{
  double v1, v2, s;

  do {
    v1 = xrandom(-1.0, 1.0);
    v2 = xrandom(-1.0, 1.0);
    s = v1*v1 + v2*v2;
  } while (s >= 1.0);
  return (mean + sdev * v1 * sqrt(-2.0 * log(s) / s));
}

#if defined(TESTBED)

string defv[] = {
  "seed=123",
  "count=10",
  NULL,
};

int main(int argc, string argv[])
{
  int n;
  double x;

  initparam(argv, defv);
  if (! strnull(getparam("seed")))
    init_random(getiparam("seed"));
  for (n = getiparam("count"); n > 0; n--) {
    x = xrandom(0.0, 1.0);
    printf("xrandom(0,1) -> %.15f (%a)\n", x, x);
  }
  return (0);
}

#endif
