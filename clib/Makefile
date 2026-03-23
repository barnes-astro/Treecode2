##  .../clib/Makefile: build and test C library package for export treecode.
##  Copyright (c) 2026, Joshua E. Barnesm Honolulu, Hawai'i.
##  ________________________________________________________________________

OBJFILES = allocate.o checkfmt.o cputime.o error.o getparam.o\
	   mathfns_d.o mathfns_f.o pickpnt_d.o pickpnt_f.o\
	   random.o scanopt.o stropen.o vectmath_d.o vectmath_f.o

TESTCODE = checkfmt_test error_test getparam_test\
	   mathfns_test_d mathfns_test_f random_test scanopt_test\
	   stropen_test vectmath_test_d vectmath_test_f

##  Define compilation and loading options.
##  _______________________________________

ZCC = gcc
ZCCFLAGS = -std=gnu17 -DLINUX
ZLDFLAGS = 

##  Optimization for CPU-intensive packages: -O2 keeps code size down.
##  __________________________________________________________________

OPT = -O2

# Assemble object library.

libClib.a: $(OBJFILES)
	ar ruv libClib.a $(OBJFILES)

tidy:
	rm -f $(OBJFILES) $(TESTCODE)

# Create individual object files.

allocate.o: allocate.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c allocate.c

checkfmt.o: checkfmt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c checkfmt.c

cputime.o: cputime.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c cputime.c

error.o: error.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c error.c

getparam.o: getparam.c getparam.h stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c getparam.c

mathfns_d.o: mathfns.h mathfns.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c mathfns.c -DDOUBLEPREC -o mathfns_d.o

mathfns_f.o: mathfns.h mathfns.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c mathfns.c -DSINGLEPREC -o mathfns_f.o

pickpnt_d.o: pickpnt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c pickpnt.c -DDOUBLEPREC -o pickpnt_d.o

pickpnt_f.o: pickpnt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c pickpnt.c -DSINGLEPREC -o pickpnt_f.o

random.o: random.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c random.c

scanopt.o: scanopt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c scanopt.c

stropen.o: stropen.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c stropen.c

vectmath_d.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c vectmath.c -DDOUBLEPREC -o vectmath_d.o

vectmath_f.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c vectmath.c -DSINGLEPREC -o vectmath_f.o

# Test library routines.

tests: libClib.a
	$(MAKE) $(TESTCODE)

checkfmt_test: checkfmt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o checkfmt_test -DTESTBED checkfmt.c \
	       libClib.a

error_test: error.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o error_test -DTESTBED error.c libClib.a

getparam_test: getparam.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o getparam_test -DTESTBED getparam.c \
	       libClib.a

mathfns_test_d: mathfns.c mathfns.h stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o mathfns_test_d -DDOUBLEPREC -DTESTBED \
	       mathfns.c libClib.a -lm

mathfns_test_f: mathfns.c mathfns.h stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o mathfns_test_f -DSINGLEPREC -DTESTBED \
	       mathfns.c libClib.a -lm

random_test: random.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o random_test -DTESTBED random.c \
	       libClib.a -lm

scanopt_test: scanopt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o scanopt_test -DTESTBED scanopt.c \
	       libClib.a

stropen_test: stropen.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o stropen_test -DTESTBED stropen.c \
	       libClib.a

vectmath_test_d: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o vectmath_test_d -DDOUBLEPREC \
	       vectmath_test.c libClib.a -lm

vectmath_test_f: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -o vectmath_test_f -DSINGLEPREC \
	       vectmath_test.c libClib.a -lm

# Generate listing of source code.

zeno_clib.pdf: allocate.c checkfmt.c cputime.c error.c getparam.c getparam.h \
	       mathfns.c mathfns.h pickpnt.c random.c scanopt.c stdinc.h \
	       stropen.c vectdefs.h vectmath.c vectmath.h
	enscript -M Letterdj -r2 \
	      allocate.c checkfmt.c cputime.c error.c getparam.c getparam.h \
	      mathfns.c mathfns.h pickpnt.c random.c scanopt.c stdinc.h \
	      stropen.c vectdefs.h vectmath.c vectmath.h -o | \
	  ps2pdf - zeno_clib.pdf
