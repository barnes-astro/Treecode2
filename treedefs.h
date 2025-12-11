/*
 * treedefs.h: include file for hierarchical force calculation routines.
 * These definitions are needed for treebuild.c and treegrav.c.
 * Copyright (c) 2025 by Joshua E. Barnes, Honolulu, Hawaii.
 */
 
#ifndef _treedefs_h
#define _treedefs_h

//  Body and cell data structures together represent the tree.  During
//  oct-tree construction, pointers are stored in the desc arrays:
// 
//           +-----------------------------------------------+
//  root --> | CELL: next, mass, pos, desc:[/,o,/,/,/,/,o,/] |
//           +--------------------------------|---------|----+
//                                            |         |
//      +-------------------------------------+         |
//      |                                               |
//      |    +--------------------------------------+   |
//      +--> | BODY: next, mass, pos, vel, acc, phi |   |
//           +--------------------------------------+   |
//                                                      |
//      +-----------------------------------------------+
//      |
//      |    +-----------------------------------------------+
//      +--> | CELL: next, mass, pos, desc:[o,/,/,o,/,/,o,/] |
//           +------------------------------|-----|-----|----+
//                                         etc   etc   etc
//
//  After the oct-tree is built, it's threaded using the next and more
//  pointers, and the desc array is reused for rcrit2, more, and quad:
// 
//           +-----------------------------------------------+
//  root --> | CELL: next:/, mass, pos, rcrit2, more:o, quad |
//           +---------------------------------------|-------+
//                                                   |
//      +--------------------------------------------+
//      |
//      |    +----------------------------------------+
//      +--> | BODY: next:o, mass, pos, vel, acc, phi |
//           +------------|---------------------------+
//                        |
//      +-----------------+
//      |
//      |    +-----------------------------------------------+
//      +--> | CELL: next:/, mass, pos, rcrit2, more:o, quad |
//           +---------------------------------------|-------+
//                                                  etc

//  node: data common to body and cell structures.
//  ______________________________________________

typedef struct _node {
  byte type;                    // code for node type
  byte level;			// for cells, depth in tree
  unsigned short active;        // status in force calc
  struct _node *next;		// link to next force calc
  real mass;                    // total mass of node
  vector pos;                   // position of node
} node, *nodeptr;
 
#define Type(x)    (((nodeptr) (x))->type)
#define Level(x)   (((nodeptr) (x))->level)
#define Active(x)  (((nodeptr) (x))->active)
#define Next(x)    (((nodeptr) (x))->next)
#define Mass(x)    (((nodeptr) (x))->mass)
#define Pos(x)     (((nodeptr) (x))->pos)

#define BodyType    0x40        // type code for bodies
#define CellType    0x80        // type code for cells
#define OctFlag     0x01	// cell in oct tree
#define ThreadFlag  0x02	// cell in threaded tree

#define Cell(x)   ((Type(x) & CellType) != 0)
#define Body(x)   ((Type(x) & BodyType) != 0)

//  body: data structure used to represent simulation bodies.
//  _________________________________________________________
 
typedef struct {
  node bodynode;                // data common to all nodes
  vector vel;                   // velocity of body
  vector acc;                   // acceleration of body
  real phi;                     // potential at body
} body, *bodyptr;

#define Vel(x)    (((bodyptr) (x))->vel)
#define Acc(x)    (((bodyptr) (x))->acc)
#define Phi(x)    (((bodyptr) (x))->phi)

//  cell: structure used to represent internal nodes of tree.
//  _________________________________________________________
  
#define NDESC (1 << NDIM)       // descendents of cell
 
typedef struct {
  node cellnode;                // data common to all nodes
  union {
    nodeptr desc[NDESC];	// links to descendents of cell
    struct {
#ifndef MODTREECODE
      real rcrit2;              // squared crit. dist. from cell's c-of-m
#else
      real rcrit;               // modified code uses actual distance
#endif
      nodeptr more;		// link to first descendent
      symatrix quad;		// traceless, symmetric quad. moment
      real quadEps2;		// softening corr: TraceQ * eps^2 / 3
    } cellprop;
  } dorp;
} cell, *cellptr;
 
#define Desc(x)   (((cellptr) (x))->dorp.desc)
#ifndef MODTREECODE
#define Rcrit2(x) (((cellptr) (x))->dorp.cellprop.rcrit2)
#else
#define Rcrit(x)  (((cellptr) (x))->dorp.cellprop.rcrit)
#define Rcrit2(x) (Rcrit(x) * Rcrit(x))
#endif
#define More(x)   (((cellptr) (x))->dorp.cellprop.more)
#define Quad(x)   (((cellptr) (x))->dorp.cellprop.quad)
#define QuadEps2(x)  (((cellptr) (x))->dorp.cellprop.quadEps2)

//  tree: structure used to hold tree and properties.
//  _________________________________________________

typedef struct {
  real eps;			// Plummer softening length
  real theta;			// force accuracy parameter
  bool usequad;			// apply quadrupole corrections
  string options;		// option string for internal reference
  matrix rotn;			// rotation for tree coordinates
  real scale;			// scale for tree coordinates
  vector trans;			// offset for tree coordinates
  cellptr root;			// pointer to root cell
  real rsize;			// side-length of root cell
  int nbody;			// number of bodies in tree
  int ncell;			// number of cells in tree
  int nlevel;			// number of levels in tree
  real cpubuild;		// CPU time for tree build
  int ntrwalk;			// number of tree walks
  long nbbcalc;			// number of body-body interactions
  long nbccalc;			// number of body-cell interactions
  real cpugrav;			// CPU time for force calculation
  void *treemisc;		// miscellaneous storage for extensions
} tree, *treeptr;

#define Eps(x)         ((x)->eps)
#define Theta(x)       ((x)->theta)
#define Usequad(x)     ((x)->usequad)
#define Options(x)     ((x)->options)
#define Rotn(x)        ((x)->rotn)
#define Scale(x)       ((x)->scale)
#define Trans(x)       ((x)->trans)
#define Root(x)        ((x)->root)
#define Rsize(x)       ((x)->rsize)
#define Nbody(x)       ((x)->nbody)
#define Ncell(x)       ((x)->ncell)
#define Nlevel(x)      ((x)->nlevel)
#define CPUbuild(x)    ((x)->cpubuild)
#define Ntrwalk(x)     ((x)->ntrwalk)
#define Nbbcalc(x)     ((x)->nbbcalc)
#define Nbccalc(x)     ((x)->nbccalc)
#define CPUgrav(x)     ((x)->cpugrav)
#define Treemisc(x)    ((x)->treemisc)

//  treebuild: construct tree from body array (treebuild.c).
//  ________________________________________________________

treeptr treebuild(bodyptr btab, int nbody, real eps, real theta, bool quad,
		  bool rotate, bool scale, real translate, string options);

//  treedestroy: collect memory used for tree (treebuild.c).
//  ________________________________________________________

void treedestroy(treeptr tptr);

//  treedump: dump tree structure in text form (treebuild.c).
//  _________________________________________________________

void treedump(treeptr tptr, stream ostr);

//  treegrav: compute gravitational forces (treegrav.c or treegrav_mod.c).
//  ______________________________________________________________________

#ifndef MODTREECODE
void treegrav(treeptr tptr);
#else
void treegrav(treeptr tptr, int nshare);
#endif

//  treegrav1: compute gravitational force on one body (treegrav.c).
//  ________________________________________________________________

void treegrav1(bodyptr p, treeptr tptr, bool selfcheck);

//  gravtrace: report diagnostics of gravity calculation (treegrav.c).
//  __________________________________________________________________

#if defined(GRAVTRACE)
void gravtrace(bodyptr p, nodeptr q, real phiM, real phiQ, vector dr,
	       real mdr3i, vector qdr, real dr5i, treeptr tptr);
#endif

#endif // ! _treedefs_h */
