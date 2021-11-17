library(glmnet)
big_lrt <- function(Xunfold, y, theta, dict, p, 
                    rhoinit = 1, maxit = 10, eps = 1e-6,
                    rhomin = 1e-6, rhomax = 1e6) {
  
  Dmat <- acos(crossprod(theta)) # theta should be unit norm vectors
  y <- matrix(y, nrow = nrow(theta))
  X <- dense_sp_mult(dict, Xunfold, p)
  nd <- nrow(Dmat)
  
  big_lm <- function(rhoinit, returnall = FALSE) {
    dmat <- exp( -Dmat / rhoinit) # covariance
    E <- eigen(dmat, TRUE)
    Sighalf <- E$vectors %*% diag(1 / sqrt(E$values), nd, nd)  
    ytilde <- Sighalf %*% y
    dim(ytilde) <- c(length(ytilde), 1)
    SigDict <- Sighalf %*% dict
    xtilde <- dense_sp_mult(SigDict, Xunfold, p)
    fit <- bigGlm(xtilde, ytilde, lower.limits = 0)
    bhat <- coef(fit)
    if (!returnall) return(bhat)
    else return(
      list(bhat = bhat, resids = ytilde - xtilde %*% bhat[-1] - bhat[1])
    )
  }
  negll <- function(rho) {
    Sinv <- solve(exp(-Dmat / rho))
    tr <- sum(t(Sinv) * Scatter)
    logdet <- determinant(Sinv)$modulus * ncol(resids)
    ll <- tr - logdet
    attributes(ll) <- NULL
    return(ll)
  }
   
  for (iter in 1:maxit) {
    print(iter)
    bhat <- big_lm(rhoinit)
    resids <- y - drop(X %*% bhat[-1]) - bhat[1]
    Scatter <- tcrossprod(resids)
    rho <- optimise(negll, interval = c(rhomin, rhomax), tol = tol)$minimum
    if (abs(rho - rhoinit) < eps) break
    rhoinit <- rho
  }
  
  
  null_model <- big_lm(rhomin, TRUE)
  sig_null <- sqrt(mean(null_model$resids^2))
  nulldev <- mean(dnorm(null_model$resids, 0, sig_null, log = TRUE))
  mle <- big_lm(rho, TRUE)
  sig_mle <- sqrt(mean(mle$resids^2))
  mledev <- mean(dnorm(mle$resids, 0, sig_mle, log = TRUE))
  return(null_model = null_model$bhat,
         mle = mle$bhat,
         nulldev = nulldev,
         mledev = mledev,
         rhohat = rho)
}

# dnegll <- function(rho) { # almost certainly wrong
#   Sinv <- solve(exp(-Dmat / rho))
#   tr <- sum(diag(Scatter %*% Sinv %*% Dmat)) / rho^2
#   logdet <- sum(diag(Dmat)) * ncol(resids) / rho^2
#   return(- tr - logdet)
# }