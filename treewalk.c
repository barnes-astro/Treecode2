/*
 * treewalk.c: tree-walk code common to treegrav.c and treegrav_mod.c.
 */

//  treewalkQ: compute acceleration on body p including quadrupole terms.
//  Returns TRUE if self-interaction detected and excluded.
//  _____________________________________________________________________

local bool treewalkQ(bodyptr p, treeptr tptr)
{
  nodeptr q = (nodeptr) Root(tptr);		// start at root of tree
  vector pos_p, acc_p, dr, qdr;
  real phi_p, eps2, dr2, dr2i, dr1i, phiM, mdr3i, qdr2, dr5i, phiQ;
  bool selfskip = FALSE;

  SETV(pos_p, Pos(p));				// use local vars for speed
  phi_p = 0.0;
  CLRV(acc_p);
  eps2 = Eps(tptr) * Eps(tptr);
  while (q != NULL) {				// loop scanning tree
    DOTPSUBV(dr2, dr, Pos(q), pos_p);
    if (Body(q) || Rcrit2(q) < dr2) {		// check interaction rules
      if (p == (bodyptr) q) {
	// eprintf("[%s.treewalkQ: skipping self-interaction]\n", getprog());
	selfskip = TRUE;
      } else {					// compute force of q on p
	dr2i = ((real) 1.0) / (dr2 + eps2);
	dr1i = rsqrt(dr2i);
	phiM = Mass(q) * dr1i;
	mdr3i = phiM * dr2i;
	if (Body(q)) {
	  phi_p -= phiM;
	  ADDMULVS(acc_p, dr, mdr3i);
	  Nbbcalc(tptr)++;
	} else {				// compute quad. corrections
	  DOTPMULSMV(qdr2, qdr, Quad(q), dr);
#if !defined(NOSOFTCORR)
	  qdr2 -= QuadEps2(q);			// apply Keigo's correction
#endif
	  dr5i = dr2i * dr2i * dr1i;		// form inverse 5th power
	  phiQ = ((real) 0.5) * dr5i * qdr2;	// get quad potential
	  phi_p -= phiM + phiQ;			// add mono and quad pot
	  mdr3i += ((real) 5.0) * phiQ * dr2i;	// adjust radial term
	  ADDMULVS2(acc_p, dr, mdr3i, qdr, - dr5i);	// mono and quad acc
	  Nbccalc(tptr)++;
	}
      }
      q = Next(q);				// done, onto next node
    } else {
      q = More(q);				// need to dig deeper
    }
  }
  Phi(p) = phi_p;				// store results in body
  SETV(Acc(p), acc_p);
  Ntrwalk(tptr)++;				// count tree walk
  return (selfskip);
}

//  treewalkM: compute acceleration on body p using monopole terms.
//  Returns TRUE if self-interaction detected and excluded.
//  _______________________________________________________________

local bool treewalkM(bodyptr p, treeptr tptr)
{
  nodeptr q = (nodeptr) Root(tptr);		// start at root of tree
  vector pos_p, acc_p, dr;
  real phi_p, eps2, dr2, dr2i, dr1i, phiM, mdr3i;
  bool selfskip = FALSE;

  SETV(pos_p, Pos(p));				// use local vars for speed
  phi_p = 0.0;
  CLRV(acc_p);
  eps2 = Eps(tptr) * Eps(tptr);
  while (q != NULL) {				// loop scanning tree
    DOTPSUBV(dr2, dr, Pos(q), pos_p);
    if (Body(q) || Rcrit2(q) < dr2) {		// check interaction rules
      if (p == (bodyptr) q) {
	// eprintf("[%s.treewalkM: skipping self-interaction]\n", getprog());
	selfskip = TRUE;
      } else {					// compute force of q on p
	dr2i = ((real) 1.0) / (dr2 + eps2);
	dr1i = rsqrt(dr2i);
	phiM = Mass(q) * dr1i;
	mdr3i = phiM * dr2i;
	phi_p -= phiM;
	ADDMULVS(acc_p, dr, mdr3i);
	if (Body(q)) {
	  Nbbcalc(tptr)++;			// count body-body terms
	} else
	  Nbccalc(tptr)++;			// count body-cell terms
      }
      q = Next(q);				// done, go on to next node
    } else {
      q = More(q);				// need to dig deeper
    }
  }
  Phi(p) = phi_p;				// store results in body
  SETV(Acc(p), acc_p);
  Ntrwalk(tptr)++;				// count tree walk
  return (selfskip);
}
