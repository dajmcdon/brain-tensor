require(glmnet)
require(tidyverse)
require(GpGp)


big_lrt <- function(Xunfold, y, theta, dict, p, 
                    start_vals = NULL, maxit = 10, eps = 1e-6,
                    rhomin = 1e-6, rhomax = 1e6) {
  
  stopifnot(is.vector(y))
  Dmat <- fields::rdist(theta)
  n <- length(y)
  nd <- nrow(Dmat)
  y <- matrix(y, nrow = nd)
  nvox <- ncol(y)
  X <- dense_sp_mult(dict, Xunfold, n, p)
  
  big_lm <- function(pars) {
    dmat <- GpGp::exponential_isotropic(pars, theta)
    # not very stable
    E <- eigen(dmat, TRUE)
    Sighalf <- E$vectors %*% Matrix::diag(1 / sqrt(E$values), nd, nd)
    ytilde <- crossprod(Sighalf, y)
    dim(ytilde) <- c(length(ytilde), 1)
    SigDict <- crossprod(Sighalf, dict)
    xtilde <- dense_sp_mult(SigDict, Xunfold, n, p)
    fit <- glmnet::bigGlm(xtilde, ytilde, lower.limits = 0)
    bhat <- coef(fit)
    bhat
  }
  
  negll <- function(logpars, penalized = TRUE) {
    # per observation, times -1
    # input on logscale since used by optim
    pars <- exp(logpars)
    dmat <- GpGp::exponential_isotropic(pars, theta)
    Sinv <- solve(dmat)
    tr <- sum(Sinv * Scatter) / (2 * n)
    logdet <- determinant(dmat)$modulus * nvox / (2 * n)
    ll <- tr + logdet + log(2*pi) / 2 + sum(pen(pars))*(penalized)
    attributes(ll) <- NULL
    return(ll)
  }
  
  null_bhat <- big_lm(c(1, rhomin, 0))
  resids <- y - drop(X %*% null_bhat[-1]) - null_bhat[1]
  vy <- mean(resids^2)
  Scatter <- tcrossprod(resids)
  nulldev <- negll(log(c(vy, rhomin, 0)), FALSE)
  
  pen <- penalty(vy)
  if (is.null(start_vals)) start_vals <- c(vy, mean(Dmat / 4), 0.1)
  stopifnot(length(start_vals) == 3)
  pars <- start_vals
  fn <- nulldev
  
  for (iter in 1:maxit) {
    print(str_glue("On iteration {iter}, pars are: ",
                   "{paste(round(pars,3), collapse = ', ')}.\n"))
    bhat <- big_lm(pars)
    resids <- y - drop(X %*% bhat[-1]) - bhat[1]
    Scatter <- tcrossprod(resids) 
    logpars <- log(pars) # convert to real valued for optim
    out <- optim(logpars, negll)
    if (sum(abs(logpars - out$par)) < eps) break
    if (abs(out$value - fn) < eps) break
    pars <- exp(out$par)
    fn <- out$value
  }
  pars <- exp(out$par)
  mledev <- negll(out$par, FALSE)
  
  return(list(null_bhat = null_model$bhat,
         mle_bhat = mle$bhat,
         mle_pars = pars,
         nulldev = nulldev,
         mledev = mledev,
         lrt = 2*(nulldev - mledev)))
}


intexpit <- function(x) log(1 + exp(x))
pen_loglo <- function(x, tt, aa) {
  if (x == 0) return(0.0) 
  else return(pen_lo(log(x), tt, aa))
}
pen_hi <- function(x, tt, aa) tt * intexpit(x - aa)
pen_lo <- function(x, tt, aa) tt * intexpit(aa - x)

penalty <- function(vv) {
  # variance, range, nugget; from GpGp exponential_isotropic
  function(x) {
    c(pen_hi(x[1] / vv, 1, 6), x[2], pen_loglo(x[3], .1, log(0.001)))
  }
}

link_fun <- function(logpars) exp(logpars)
invlink_fun <- function(exppars) log(pars)
