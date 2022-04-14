# Compute the log-posterior = log-likelihood + log-prior.
logposterior_multinom <- function (x, F, L, alpha, delta, e = 1e-8) {
  m <- nrow(F)
  k <- ncol(F)
  x <- sparseMatrix(i = x$i,j = x$j,x = x$v,dims = c(x$nrow,x$ncol))
  return(sum(loglik_multinom_const(x)
             - cost(x,L,t(F),e)
             + ddiri(L,rep(alpha,k),logged = TRUE))
         + sum(ddiri(t(F),rep(delta,m),logged = TRUE)))
}

# Computes loss function equal to negative log-likelihood, ignoring
# terms that do not depend on F and L.
cost <- function (X, A, B, e = 1e-8, version = c("R","cpp")) {
  version <- match.arg(version)
  if (is.matrix(X) | version == "R")
    y <- -rowSums(X*log(A %*% B + e))
  else {
    x <- summary(X)
    y <- .Call(C_rcost,as.integer(x$i),as.integer(x$j),as.numeric(x$x),
               as.numeric(e))
  }
  return(y)
}

# Compute the constant terms in the multinomial log-likelihoods.
loglik_multinom_const <- function (x)
  lgamma(rowSums(x) + 1) + loglik_poisson_const(x)

# Compute the constant terms in the Poisson log-likelihoods.
loglik_poisson_const <- function (x) {
  if (is.matrix(x))
    return(-rowSums(lgamma(x + 1)))
  else
    return(-rowSums(apply.nonzeros(x,function (x) lgamma(x + 1))))
}

# Apply operation f to all nonzeros of a sparse matrix.
apply.nonzeros <- function (x, f) {
  d <- summary(x)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(x)))
}
