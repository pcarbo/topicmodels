#include <R.h>
#include <Rinternals.h>

// Implements core calculations for cost(X,A,B,e,version = "R") when X
// is a sparse matrix.
SEXP rcost (SEXP np, SEXP mp, SEXP kp, SEXP ip, SEXP jp, SEXP xp, SEXP ep) {
  int     n;
  int     m;
  int     k;
  int*    i;
  int*    j;
  double  ab;
  double  e;
  double* x;
  double* y;
  SEXP    yp;

  // Set up access to the inputs.
  np = PROTECT(coerceVector(np,INTSXP));
  mp = PROTECT(coerceVector(mp,INTSXP));
  kp = PROTECT(coerceVector(kp,INTSXP));
  ip = PROTECT(coerceVector(ip,INTSXP));
  jp = PROTECT(coerceVector(jp,INTSXP));
  xp = PROTECT(coerceVector(xp,REALSXP));
  ep = PROTECT(coerceVector(ep,REALSXP));
  n  = *INTEGER(np);
  m  = *INTEGER(mp);
  k  = *INTEGER(kp);
  i  = INTEGER(ip);
  j  = INTEGER(jp);
  x  = REAL(xp);
  e  = *REAL(ep);

  // Initialize the output.
  yp = PROTECT(allocVector(REALSXP,1));
  y  = REAL(yp);

  // This is equivalent to the following R code:
  //
  //   y = -X[i,j] * log(ab + e)
  //
  // where 
  // 
  //   ab = sum(A[i,] * B[,j])
  //
  n = length(xp);
  *y = 0;
  for (int i = 0; i < n; i++) {
    // ab = A * B.col(j);
    ab = 1;
    *y -= x[i] * log(ab + e);
  }

  UNPROTECT(8);
  return yp;
}

SEXP convolve2 (SEXP a, SEXP b) {
  int na, nb, nab;
  double *xa, *xb, *xab;
  SEXP ab;

  a = PROTECT(coerceVector(a,REALSXP));
  b = PROTECT(coerceVector(b,REALSXP));
  na = length(a); 
  nb = length(b); 
  nab = na + nb - 1;
  ab = PROTECT(allocVector(REALSXP,nab));
  xa = REAL(a); 
  xb = REAL(b); 
  xab = REAL(ab);

  for (int i = 0; i < nab; i++) 
    xab[i] = 0.0;
  for (int i = 0; i < na; i++)
    for (int j = 0; j < nb; j++) 
      xab[i + j] += xa[i] * xb[j];

  UNPROTECT(3);
  return ab;
}
