This is a stand-alone implementation of Treecode2 (Barnes 2026; MNRAS, accepted Feb 2, 2026), a gravitational N-body simulation code.  Like previous oct-tree codes, it constructs a hierarchical representation of the mass distribution by enclosing the particles in a root cube, which is then recursively subdivided until each particle is isolated in a unique sub-cell.  Treecode2 can randomly reorient, reposition, and resize the root cube at each timestep.  This improves force accuracy and global conservation of energy, linear momentum, and angular momentum.

This repository consists of a top-level directory containing the actual treecode sources, and a sub-directory named 'clib' containing a small library of useful routines, adapted from the Zeno N-body software environment.

To build the treecode, first change to the 'clib' directory and type 'make' to build 'libClib.a'.  Then change to the top-level directory and type 'make' again to build 'treecode' and 'treecode_mod' (which is usually somewhat faster and more accurate).  You can also build a double-precision version, as well as a version w/o corrected softening; see the top-level 'Makefile' for details.

Before running the treecode, you may want to set the shell environment variable 'ZENO_MSG_OPTION=warn'.  This will suppress debugging messages, but not warnings and actual error messages.

To run the treecode, just type './treecode' (or './treecode_mod').  This will start a short calculation using Plummer model initial conditions and reasonable default parameters.  Type './treecode -help' to see the default values and brief documentation ('-clue' and '-explain' yield shorter and longer versions, respectively).  Command-line arguments can be used to change these defaults, matching either by position or by name=value pairs.

The present version of the code reads and writes simulation data using CSV files.  Each particle is represented by a single line giving its mass, position, and velocity, in that order.  Output files use 'a' formatting to preserve full precision, while input files may have any floating-point format.  Lines starting with '#' are comments.

Please send questions and feedback to barnes@hawaii.edu.
Josh Barnes
