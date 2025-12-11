##  .../clib/Makefile: build and test C library package for public treecode.
##  Copyright (c) 2025  Joshua E. Barnes  Honolulu, Hawai'i.

OBJFILES = allocate.o checkfmt.o cputime.o error.o getparam.o\
	   mathfns_d.o mathfns_f.o pickpnt_d.o pickpnt_f.o\
	   random.o scanopt.o vectmath_d.o vectmath_f.o

TESTCODE = checkfmt_test error_test getparam_test\
	   mathfns_test_d mathfns_test_f random_test scanopt_test\
	   vectmath_test_d vectmath_test_f

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
	$(ZCC) $(ZCCFLAGS) -c allocate.c

checkfmt.o: checkfmt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c checkfmt.c

cputime.o: cputime.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c cputime.c

error.o: error.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c error.c

getparam.o: getparam.c getparam.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -c getparam.c

mathfns_d.o: mathfns.h mathfns.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c mathfns.c -DDOUBLEPREC -o mathfns_d.o

mathfns_f.o: mathfns.h mathfns.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c mathfns.c -DSINGLEPREC -o mathfns_f.o

pickpnt_d.o: pickpnt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c pickpnt.c -DDOUBLEPREC -o pickpnt_d.o

pickpnt_f.o: pickpnt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c pickpnt.c -DSINGLEPREC -o pickpnt_f.o

random.o: random.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c random.c

scanopt.o: scanopt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c scanopt.c

vectmath_d.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c vectmath.c -DDOUBLEPREC -o vectmath_d.o

vectmath_f.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c vectmath.c -DSINGLEPREC -o vectmath_f.o

# Test library routines.

tests: libClib.a
	$(MAKE) $(TESTCODE)

checkfmt_test: checkfmt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o checkfmt_test -DTESTBED checkfmt.c libClib.a

error_test: error.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o error_test -DTESTBED error.c libClib.a

getparam_test: getparam.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o getparam_test -DTESTBED getparam.c libClib.a

mathfns_test_d: mathfns.c mathfns.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -o mathfns_test_d -DDOUBLEPREC -DTESTBED mathfns.c \
	   libClib.a -lm

mathfns_test_f: mathfns.c mathfns.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -o mathfns_test_f -DSINGLEPREC -DTESTBED mathfns.c \
	   libClib.a -lm

random_test: random.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o random_test -DTESTBED random.c libClib.a \
	    -lgsl -lgslcblas -lm

scanopt_test: scanopt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o scanopt_test -DTESTBED scanopt.c libClib.a

vectmath_test_d: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o vectmath_test_d -DDOUBLEPREC vectmath_test.c \
	   libClib.a -lm

vectmath_test_f: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o vectmath_test_f -DSINGLEPREC vectmath_test.c \
	   libClib.a -lm

# Generate listing of source code.

zeno_clib.pdf: stdinc.h allocate.c bessel.c burststring.c checkfmt.c cputime.c \
	      datatypes.h datatypes.c error.c extstring.c filestruct.h \
	      filestruct.c getparam.h getparam.c \
	      history.c mathfns.h mathfns.c phatstruct.h phatstruct.c \
	      pickpnt.c random.c scanopt.c setrange.c stropen.c \
	      strset.h strset.c vectdefs.h vectmath.h vectmath.c \
	      vectmath_test.c within.c
	enscript -M Letterdj -r2 \
	      stdinc.h allocate.c bessel.c burststring.c checkfmt.c cputime.c \
	      datatypes.h datatypes.c error.c extstring.c filestruct.h \
	      filestruct.c getparam.h getparam.c \
	      history.c mathfns.h mathfns.c phatstruct.h phatstruct.c \
	      pickpnt.c random.c scanopt.c setrange.c stropen.c \
	      strset.h strset.c vectdefs.h vectmath.h vectmath.c \
	      vectmath_test.c within.c -o - | \
	  ps2pdf - zeno_clib.pdf
