#include <R.h>
#include <Rinternals.h>

// Implements core calculations for cost(X,A,B,e,version = "R") when X
// is a sparse matrix.
SEXP rcost (SEXP np, SEXP mp, SEXP kp, SEXP ip, SEXP jp, SEXP xp, 
	    SEXP ap, SEXP bp, SEXP ep) {
  int     n;
  int     m;
  int     k;
  int     t;
  int     N;
  int     K;
  int     it;
  int     jt;
  int*    i;
  int*    j;
  double  ab;
  double  e;
  double* x;
  double* a;
  double* b;
  double* ai;
  double* bj;
  double* y;
  SEXP    yp;

  // Set up access to the inputs.
  np = PROTECT(coerceVector(np,INTSXP));
  mp = PROTECT(coerceVector(mp,INTSXP));
  kp = PROTECT(coerceVector(kp,INTSXP));
  ip = PROTECT(coerceVector(ip,INTSXP));
  jp = PROTECT(coerceVector(jp,INTSXP));
  xp = PROTECT(coerceVector(xp,REALSXP));
  ap = PROTECT(coerceVector(ap,REALSXP));
  bp = PROTECT(coerceVector(bp,REALSXP));
  ep = PROTECT(coerceVector(ep,REALSXP));
  n  = *INTEGER(np);
  m  = *INTEGER(mp);
  K  = *INTEGER(kp);
  i  = INTEGER(ip);
  j  = INTEGER(jp);
  x  = REAL(xp);
  a  = REAL(ap);
  b  = REAL(bp);
  e  = *REAL(ep);

  // Initialize the output.
  yp = PROTECT(allocVector(REALSXP,n));
  y  = REAL(yp);

  // This is equivalent to the following R code:
  //
  //   y = -X[i,j] * log(ab + e)
  //
  // where 
  // 
  //   ab = sum(A[i,] * B[,j])
  //
  N = length(xp);
  for (t = 0; t < n; t++)
    y[t] = 0;
  for (t = 0; t < N; t++) {
    it = i[t];
    jt = j[t];
    ab = 0;
    ai = a + it*K;
    bj = b + jt*K;
    for (k = 0; k < K; k++)
      ab += ai[k] * bj[k];
    y[it] -= x[t] * log(ab + e);
  }

  UNPROTECT(10);
  return yp;
}
