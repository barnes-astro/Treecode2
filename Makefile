##  zeno/src/nbody/tree/treecode2.1/export/Makefile: compile tree-code.
##  Copyright (c) 2025, Joshua E. Barnes, Honolulu, Hawaii.
##  ____________________________________________________________________

BINFILES = treecode$X treecode_mod$X

##  Compilation options are selected by changing OPTIONS (and X).
##  To compile BINFILES in double precision, with "_dp" suffixes:
##	make make_all OPTIONS="-DDOUBLEPREC" X=_dp
##  To compile treecode without softening correction:
##      make treecode_nc OPTIONS="-DNOSOFTCORR" X=_nc
##  _____________________________________________________________

OPTIONS =

##  Optimization flags.
##  ___________________

OPT = -O3

ZCCFLAGS = -std=gnu17 -DLINUX -Iclib

make_all: $(BINFILES)

##  Build vanilla treecode; by default, includes softening correction.
##  __________________________________________________________________

treecode$X: treecode$X.o treeio$X.o treebuild$X.o treegrav$X.o
	$(ZCC) $(ZLDFLAGS) -o treecode$X treecode$X.o \
	    treeio$X.o treebuild$X.o treegrav$X.o \
	    -lNBody -lClib -lgsl -lgslcblas -lm

treecode$X.o: treecode.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -o treecode$X.o -c treecode.c

treeio$X.o: treeio.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -o treeio$X.o -c treeio.c

treebuild$X.o: treebuild.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -o treebuild$X.o -c treebuild.c

treegrav$X.o: treegrav.c treedefs.h treewalk.c
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -o treegrav$X.o -c treegrav.c

##  Build modified treecode.
##  ________________________

treecode_mod$X: treecode_mod$X.o treeio_mod$X.o treebuild_mod$X.o \
		treegrav_mod$X.o
	$(ZCC) $(ZLDFLAGS) -o treecode_mod$X treecode_mod$X.o \
	    treeio_mod$X.o treebuild_mod$X.o treegrav_mod$X.o \
	    -lNBody -lClib -lgsl -lgslcblas -lm

treecode_mod$X.o: treecode.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -DMODTREECODE \
	    -o treecode_mod$X.o -c treecode.c

treeio_mod$X.o: treeio.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -DMODTREECODE \
	    -o treeio_mod$X.o -c treeio.c

treebuild_mod$X.o: treebuild.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -DMODTREECODE \
	    -o treebuild_mod$X.o -c treebuild.c

treegrav_mod$X.o: treegrav_mod.c treedefs.h treewalk.c
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -DMODTREECODE \
	    -o treegrav_mod$X.o -c treegrav_mod.c

##  Generate listings.
##  _________________

treecode.pdf: treecode.c treeio.c treecode.h treedefs.h treebuild.c \
	      treewalk.c treegrav.c treegrav_mod.c Makefile
	enscript -r2 -M Letterdj treecode.c treeio.c treecode.h treedefs.h \
		 treebuild.c treewalk.c treegrav.c treegrav_mod.c Makefile \
		 -o - | \
	  ps2pdf - treecode.pdf
