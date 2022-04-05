# Compute the log-posterior = log-likelihood + log-prior.
logposterior_multinom <- function (x, F, L, alpha, delta, e = 1e-15) {
  x <- sparseMatrix(i = x$i,j = x$j,x = x$v,dims = c(x$nrow,x$ncol))
  return(loglik_multinom_const(x) - cost(x,L,t(F),e))
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

