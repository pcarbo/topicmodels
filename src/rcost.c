#include <R.h>
#include <Rinternals.h>

// Implements core calculations for cost(X,A,B,e,version = "R") when X
// is a sparse matrix.
SEXP rcost (SEXP i, SEXP j, SEXP x, SEXP e) {
  int     n;
  int*    xi;
  int*    xj;
  double  ab;
  double  xe;
  double* xx;
  double* xy;
  SEXP    y;

  // Set up access to the inputs.
  i = PROTECT(coerceVector(i,INTSXP));
  j = PROTECT(coerceVector(j,INTSXP));
  x = PROTECT(coerceVector(x,REALSXP));
  e = PROTECT(coerceVector(e,REALSXP));
  xi = INTEGER(i);
  xj = INTEGER(j);
  xx = REAL(x);
  xe = *REAL(e);

  // Initialize the output.
  y  = PROTECT(allocVector(REALSXP,1));
  xy = REAL(y);

  // This is equivalent to the following R code:
  //
  //   y = -X[i,j] * log(ab + e)
  //
  // where 
  // 
  //   ab = sum(A[i,] * B[,j])
  //
  n = length(x);
  *xy = 0;
  for (int i = 0; i < n; i++) {
    ab = 1;
    *xy -= xx[i] * log(ab + xe);
  }

  UNPROTECT(5);
  return y;
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
