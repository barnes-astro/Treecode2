/*
 * checkfmt.c: test legality of printf format string.
 */

#include "stdinc.h"
#include "mathfns.h"
#include <ctype.h>

#define COMMA  ','
#define COLON  ':'

//  checkfmt: check user-supplied format string fmt for legal directive;
//  must match ref letter-by-letter (case-blind), with '?' as wild card.
//  ____________________________________________________________________

bool checkfmt(string fmt, string ref)
{
  int ndir = 0;
  char *p = fmt, *q, *r;

  while (*p != CNULL) {				// scan format string
    if (*p == '%') {				// at start of directive?
      ndir++;					// count directives found
      for (q = p, r = ref; *r != CNULL; q++, r++) {
	if (tolower(*q) != tolower(*r) && *r != '?')
	  return (FALSE);			// mismatch, not wild card
      }
      p = q;					// past end of directive
    } else {
      p++;					// on to next char
    }
  }
  return (ndir < 2);				// allow 0 or 1 directive(s)
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
    "fmt=snap_%04x.dat",
    "ref=%0?x",
    NULL,
};

int main(int argc, string argv[])
{
  initparam(argv, defv);
  printf("checkfmt returns %s\n",
	 checkfmt(getparam("fmt"), getparam("ref")) ? "TRUE" : "FALSE");
  return 0;
}

#endif
