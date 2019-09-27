#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

int ncomb(int a, int b)
{
  int bi, res = 0;
  if (a == 1) return b;
  for (bi = 0; bi < b; bi++) {
    res += ncomb(a - 1, bi + 1);
  }
  return res;
}

/******************************************************************************/
/*                                                                            */ 
/*                             SINGLE LOCUS                                   */ 
/*                for one sample, calculate all combinations                  */
/*      for each combination, ALL minus-one probabilities, return matrix      */
/*        grab logp (already filtered), logj (capped), plast (selected)       */
/* ALTERNATIVELY, can filter/cap/select inside the function, pass ppmat, etc  */
/*    RETURNS matrix; first row is a full set, subsequent rows are minus 1    */
/*              CAN NOT use void combProb() with .Call method                 */
/*                                                                            */ 
/******************************************************************************/

// based on mneqc2 (*** check) - timing the same as ... but incrementing is
// easier for the big program (no need to keep matrix, subset etc)
// basically make better documentation and decisions

SEXP combProb1(SEXP Rlogp, SEXP Rlogj, SEXP Rpplast)//, SEXP ncomb)
{  
  int a, b, nlast, i, jc, nr, nc;
  double *logp, *logj, *pplast, *res, spp, factn, factn1;  // sum of partial probs

  a = length(Rlogp) - 1;                 /* a = nux - 1 */
  b = length(Rlogj);                     /* nx = a + b  */
  nr = a + 2;                            /* 1 + nux (1 for full set) */
  nc = ncomb(a, b);
  logp   = REAL(Rlogp);
  logj   = REAL(Rlogj);
  pplast = REAL(Rpplast);
  //  SEXP Rres = PROTECT(allocMatrix(REALSXP, a + 2, INTEGER(ncomb)[0]));
  SEXP Rres = PROTECT(allocMatrix(REALSXP, nr, nc));
  res = REAL(Rres);

  int v[a];
  double pp[a];                         /* partial probs */

  /* log(nx!) and log((nx - 1)!) for numerators */
  factn1 = lgamma(a + b);  //  recast to double (function call conversion)
  factn  = factn1 + log(a + b);

  /* initialize vx, nxlast, ppx with nonexistent "pre-first" combination */
  v[0]  = 0;
  pp[0] = 0;
  for (i = 1; i < a; i++) {
    v[i]  = 1;
    pp[i] = logp[i];
  }
  nlast = b + 1;

  /* subsequent combinations */
  for (jc = 0; jc < nc; jc++) {  /* combinations */
    //  while (v[a - 1] < b) {              
    if (nlast > 1) {
      v[0]++;
      nlast--;
      pp[0] += logp[0] - logj[v[0] - 1];
    } else {
      i = 0;
      while (v[i] == 1) {  // since no == 1, will hit not 1 at some point before m
	i++;
      }
      nlast = v[i] - 1;  /* was 1; v[i] will change to 1; v[i + 1] incremented */
      v[i] = 1;
      v[i + 1]++;
      pp[i] = logp[i];
      pp[i + 1] += logp[i + 1] - logj[v[i + 1] - 1]; 
    }
    spp = pplast[nlast - 1]; 
    for (i = 0; i < a; i++) {
      spp += pp[i];
    }
    res[jc*nr] = spp + factn;         /* probability for the "full set" */
    for (i = 1; i < a + 1; i++) {     /* "minus one" prob start w/second row */
      res[i + jc*nr] = spp - logp[i - 1] + logj[v[i - 1] - 1] + factn1;
    }
    res[a + 1 + jc*nr] = spp - logp[a] + logj[nlast - 1] + factn1;
  }
  UNPROTECT(1);
  return Rres;  // do we need some kind of setAttrib() or something?
}

/* list: vector of probs for the "full set" nux and matrix for "minus one"  */
/* USE THIS ONE - a little faster and easier to read/use */
// to recompile
SEXP combProb2(SEXP Rlogp, SEXP Rlogj, SEXP Rpplast)
{  
  int a, b, nlast, i, jc, nux, nc;
  double *logp, *logj, *pplast, *vec, *mat, spp, factn, factn1;  // sum of partial probs
  SEXP Rres, Rvec, Rmat;

  nux = length(Rlogp);
  a   = nux - 1; 
  b   = length(Rlogj);                     /* nx = a + b  */
  nc  = ncomb(a, b);
  logp   = REAL(Rlogp);
  logj   = REAL(Rlogj);
  pplast = REAL(Rpplast);

  Rres = PROTECT(allocVector(VECSXP, 2));
  Rvec = PROTECT(allocVector(REALSXP, nc));
  Rmat = PROTECT(allocMatrix(REALSXP, a + 1, nc));
  vec = REAL(Rvec);
  mat = REAL(Rmat);

  int v[a];
  double pp[a];                         /* partial probs */

  /* log(nx!) and log((nx - 1)!) for numerators */
  factn1 = lgammafn(a + b);  //  recast to double (function call conversion)
  factn  = factn1 + log(a + b);

  /* initialize vx, nxlast, ppx with nonexistent "pre-first" combination */
  v[0]  = 0;
  pp[0] = 0;
  for (i = 1; i < a; i++) {
    v[i]  = 1;
    pp[i] = logp[i];
  }
  nlast = b + 1;

  /* subsequent combinations */
  for (jc = 0; jc < nc; jc++) {  /* combinations */
    //  while (v[a - 1] < b) {              
    if (nlast > 1) {
      v[0]++;
      nlast--;
      pp[0] += logp[0] - logj[v[0] - 1];
    } else {
      i = 0;
      while (v[i] == 1) {  // since no == 1, will hit not 1 at some point before m
	i++;
      }
      nlast = v[i] - 1;  /* was 1; v[i] will change to 1; v[i + 1] incremented */
      v[i] = 1;
      v[i + 1]++;
      pp[i] = logp[i];
      pp[i + 1] += logp[i + 1] - logj[v[i + 1] - 1]; 
    }
    spp = pplast[nlast - 1]; 
    for (i = 0; i < a; i++) {
      spp += pp[i];
    }
    vec[jc] = spp + factn;      /* probability for the "full set" */
    for (i = 0; i < a; i++) {    /* "minus one" probabilities      */
      mat[i + jc*nux] = spp - logp[i] + logj[v[i] - 1] + factn1;
    }
    mat[a + jc*nux] = spp - logp[a] + logj[nlast - 1] + factn1;
  }

  SET_VECTOR_ELT(Rres, 0, Rvec);
  SET_VECTOR_ELT(Rres, 1, Rmat);
  UNPROTECT(3);
  return Rres;  // do we need some kind of setAttrib() or something?
}











