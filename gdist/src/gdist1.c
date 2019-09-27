#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

#define MAX(a, b) ((a) > (b)  ? (a) : (b))

int ncomb(int a, int b)
{
  int bi, res = 0;
  if (a == 1) return b;
  for (bi = 0; bi < b; bi++) {
    res += ncomb(a - 1, bi + 1);
  }
  return res;
}

/* nrow(Rcomb1x) = nrow(Rcomb1y) = length(logp) = nuxy */
SEXP llik(SEXP Rcombx, SEXP Rcomby, SEXP Rcomb1x, SEXP Rcomb1y, SEXP Rlogp,
	   SEXP Rr)
{
  int nuxy, ncombx, ncomby, nr, i, jx, jy;
  double pIBD0, pIBD1;
  double *combx, *comby, *comb1x, *comb1y, *logp, *r, *lik;
  SEXP Rlik;

  nuxy   = length(Rlogp);
  ncombx = length(Rcombx);
  ncomby = length(Rcomby);
  nr     = length(Rr);
  combx  = REAL(Rcombx);
  comby  = REAL(Rcomby);
  comb1x = REAL(Rcomb1x);
  comb1y = REAL(Rcomb1y);
  if (!isNull(Rlogp)) {  // if NULL, no iterations in the loop
    logp   = REAL(Rlogp);
  }
  r      = REAL(Rr); 

  Rlik = PROTECT(allocVector(REALSXP, nr));
  lik    = REAL(Rlik); 

  for (i = 0; i < nr; i++) {
    lik[i] = 0;
  }
 
  for (jx = 0; jx < ncombx; jx++) {
    for (jy = 0; jy < ncomby; jy++) {
      pIBD0 = exp(combx[jx] + comby[jy]);
      pIBD1 = 0;
      for (i = 0; i < nuxy; i++) {
	pIBD1 += exp(logp[i] + comb1x[i + jx*nuxy] + comb1y[i + jy*nuxy]); // exp (!)
      }
      for (i = 0; i < nr; i++) {
	lik[i] += r[i]*pIBD1 + (1 - r[i])*pIBD0;
      }
    }
  }

  for (i = 0; i < nr; i++) {
    lik[i] = log(lik[i]);
  }
  UNPROTECT(1);
  return Rlik;
}

/******************************************************************************/
/*                                                                            */ 
/*               SINGLE LOCUS (alternatively, process all loci for a pair     */ 
/*                SCALABLE version: no large objects saved                    */
/*                      likelihood for a pair of samples                      */
/*                                                                            */ 
/******************************************************************************/

// version 1: for (ix = 0; ix < ncombx; ix++)
// version 2 is removed ( "while (vx[ax - 1] < bx || vx[0] == 0) {" )
/* ixy - indices of which ux are in uy (or uxy), sorted; likewise for iyx */
SEXP lliks(SEXP Rux, SEXP Ruy, SEXP Rixy, SEXP Riyx, SEXP Rnx, SEXP Rny, 
	  SEXP Rprob, SEXP Rr) {  // s for scalable
  int nux, nuy, nuxy, nr, K, nx, ny, ax, ay, bx, by, nj, i, ncombx, ncomby, ix, iy;
  int *ux, *uy, *ixy, *iyx; 
  double factnx, factny, factn1x, factn1y, combx, comby, pIBD0, pIBD1;
  double *p, *r, *lik;
  SEXP Rlik;

  nux  = length(Rux);         
  nuy  = length(Ruy);
  nuxy = length(Rixy);   // length(Rixy) = length(Riyx)
  K    = length(Rprob);
  nr   = length(Rr);
  ux   = INTEGER(Rux);   // one-based indices (passed from R)
  uy   = INTEGER(Ruy);   // one-based indices (passed from R)
  ixy  = INTEGER(Rixy);  // one-based indices (passed from R)
  iyx  = INTEGER(Riyx);  // one-based indices (passed from R)
  nx   = INTEGER(Rnx)[0];
  ny   = INTEGER(Rny)[0];
  p    = REAL(Rprob);   
  r    = REAL(Rr);

  Rlik  = PROTECT(allocVector(REALSXP, nr));
  lik   = REAL(Rlik); 
  for (i = 0; i < nr; i++) {
    lik[i] = 0;
  }

  /* rewrite ux, uy, ixy, and iyx to zero-based indices for convenience */
  for (i = 0; i < nux;  i++) ux[i]--;
  for (i = 0; i < nuy;  i++) uy[i]--;
  for (i = 0; i < nuxy; i++) {
    ixy[i]--;
    iyx[i]--;
  }

  ax = nux - 1;
  ay = nuy - 1;
  bx = nx - ax;
  by = ny - ay;
  nj = MAX(bx, by);
  ncombx = ncomb(ax, bx);
  ncomby = ncomb(ay, by); 

  int vx[nux], vy[nuy]; 
  double logj[nj], logp[K], logpx[nux], logpy[nuy], logpxy[nuxy],
    ppx[ax], ppy[ay], pplastx[bx], pplasty[by], comb1x[nuxy], comb1y[nuxy];
  
  /* log(nx!) and log((nx - 1)!) for numerators */
  factn1x = lgammafn(nx);  
  factnx  = factn1x + log(nx);
  factn1y = lgammafn(ny);  
  factny  = factn1y + log(ny);
  /* logj, logp, logpx, logpy, logpxy */
  for (i = 0; i < nj;   i++) logj[i]   = log(i + 1);
  for (i = 0; i < K;    i++) logp[i]   = log(p[i]);
  for (i = 0; i < nux;  i++) logpx[ i] = logp[ux[i]];
  for (i = 0; i < nuy;  i++) logpy[ i] = logp[uy[i]];
  for (i = 0; i < nuxy; i++) logpxy[i] = logp[ux[ixy[i]]]; 

  /* pplastx, pplasty - partial probs for last category */
  pplastx[0] = logpx[ax];
  for (i = 1; i < bx; i++) {
    pplastx[i] = pplastx[i - 1] + pplastx[0] - logj[i];
  }
  pplasty[0] = logpy[ay];
  for (i = 1; i < by; i++) {
    pplasty[i] = pplasty[i - 1] + pplasty[0] - logj[i];
  }

  /* initialize vx, nxlast, ppx with nonexistent "pre-first" combination */
  vx[0]  = 0;
  ppx[0] = 0;
  for (i = 1; i < ax; i++) {
    vx[i]  = 1;
    ppx[i] = logpx[i];
  }
  vx[ax] = bx + 1;

  /* subsequent combinations */
  for (ix = 0; ix < ncombx; ix++) {
    if (vx[ax] > 1) {
      vx[0]++;
      vx[ax]--;
      ppx[0] += logpx[0] - logj[vx[0] - 1];
    } else {
      i = 0;
      while (vx[i] == 1) {  // since no == 1, will hit != 1 at some point before m
	i++;
      }
      vx[ax] = vx[i] - 1;  /* was 1; v[i] will change to 1; v[i + 1] incremented */
      vx[i] = 1;
      vx[i + 1]++;
      ppx[i] = logpx[i];
      ppx[i + 1] += logpx[i + 1] - logj[vx[i + 1] - 1];
    }
    combx = pplastx[vx[ax] - 1]; 
    for (i = 0; i < ax; i++) {
      combx += ppx[i];
    }
    for (i = 0; i < nuxy; i++) {
      comb1x[i] = combx - logp[ux[ixy[i]]] + logj[vx[ixy[i]] - 1] + factn1x;
    }
    combx += factnx;

    /* initialize vy, nylast, ppy */
    vy[0]  = 0;
    ppy[0] = 0;
    for (i = 1; i < ay; i++) {
      vy[i]  = 1;
      ppy[i] = logpy[i];
    }
    vy[ay] = by + 1;
    /* subsequent combinations for y */ 
    for (iy = 0; iy < ncomby; iy++) {
      if (vy[ay] > 1) {
        vy[0]++;
	vy[ay]--;
        ppy[0] += logpy[0] - logj[vy[0] - 1];
      } else {
        i = 0;
        while (vy[i] == 1) {  // since no == 1, will hit not 1 at some point before m
	  i++;
        }
	vy[ay] = vy[i] - 1;
	vy[i] = 1;
        vy[i + 1]++;
        ppy[i] = logpy[i];
        ppy[i + 1] += logpy[i + 1] - logj[vy[i + 1] - 1];
      }
      comby = pplasty[vy[ay] - 1];
      for (i = 0; i < ay; i++) {
        comby += ppy[i];
      }
      for (i = 0; i < nuxy; i++) {
        comb1y[i] = comby - logpxy[i] + logj[vy[iyx[i]] - 1] + factn1y;
      }
      comby += factny;
      
      /* calculate likelihood for a combination X and combination Y */
      pIBD0 = exp(combx + comby);
      pIBD1 = 0;
      for (i = 0; i < nuxy; i++) {
        pIBD1 += exp(logpxy[i] + comb1x[i] + comb1y[i]);
      }
      for (i = 0; i < nr; i++) {
        lik[i] += r[i]*pIBD1 + (1 - r[i])*pIBD0;
      }
    }
  }

  for (i = 0; i < nr; i++) {
    lik[i] = log(lik[i]);
  }
  
  UNPROTECT(1);
  return Rlik;
}
