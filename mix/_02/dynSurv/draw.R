#' ####   mixture of GLMMs //
#' Functions for taking 'draws' from \Omega ~ N(\hat{Omega}, I(\Omega)^{-1}) and f(\b|...; \Omega).
#' ####

# sourceCpp('helper.cpp')

# Draw from \Omega ~ N(\hat{Omega}, I(\Omega)^{-1}) -----------------------
Omega.draw <- function(fit){
  Omega.mean <- setNames(c(
    vech(fit$coeffs$D),
    c(fit$coeffs$beta),
    fit$coeffs$var.e,
    fit$coeffs$gamma,
    fit$coeffs$eta
  ), names(fit$SE))
  
  Omega.vcov <- fit$vcov
  
  draw <- MASS::mvrnorm(1, Omega.mean, solve(Omega.vcov))
  
  Dinds <- grepl('^D\\[', names(draw))
  binds <- grepl('^G\\d?\\_|^B\\_', names(draw))
  sinds <- grepl('var.e', names(draw))
  ginds <- grepl('gamma', names(draw))
  einds <- grepl('^cont|^bin', names(draw))
  tinds <- grepl('^theta', names(draw))
  
  D <- vech2mat(draw[Dinds], ncol(fit$RE))
  beta <- draw[binds]
  var.e <- draw[sinds]
  gamma <- draw[ginds]
  eta <- draw[einds]
  
  return(list(
    D = D, beta = beta, var.e = var.e, gamma = gamma, eta = eta
  ))
}

# Draw b ------------------------------------------------------------------
b.draw <- function(b0, X, Y, Z, beta, var.e, D, Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta, Sigma0){
  
  b.hat <- ucminf::ucminf(b0,
                          joint_density, joint_density_ddb,
                          X, Z, beta, var.e, D, Y[, 1], Y[, 2], Y[, 3],
                          Delta, K, Fi, l0i, KK, Fu, haz, rep(gamma, each = 2), eta)$par
  
  Sigmai <- solve(joint_density_sdb(b.hat, X, Z, beta, var.e, D, 
                                    Y[, 1], Y[, 2], Y[, 3], 
                                    Delta, K, Fi, l0i,
                                    KK, Fu, haz, rep(gamma, each = 2), eta, 1e-3))
  
  # Check Sigmai is pos semi-def
  if(det(Sigmai) <= 0 || any(eigen(Sigmai)$value < 0)){
    Sigmai <- as.matrix(Matrix::nearPD(Sigmai)$mat)
  }
  
  # mvtnorm::rmvnorm(1, b.hat, Sigmai)
  list(
    b = MASS::mvrnorm(1, b.hat, Sigmai),
    S = Sigmai
  )
  
}

# Survival function -------------------------------------------------------
S <- function(b, Omega, surv){
  # b: Draw from f(\b|..., \Omega).
  # Omega: Draw from \Omega ~ N(\hat{Omega}, I(\Omega)^{-1}).
  # surv: Object created by prep.surv; dictates whether this is t (denom), or u(numer).
  
  l0u <- surv$l0u.t
  KK <- surv$KK.t
  eta <- Omega$eta
  gamma <- rep(Omega$gamma, each = 2)
  Fu <- do.call(cbind, replicate(3, surv$Fu.t, simplify = F))
  b <- b
  
  exp(-l0u %*% exp(KK %*% eta + Fu %*% (gamma * b)))
}

