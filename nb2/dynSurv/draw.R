#' ####   neg-binom //
#' Functions for taking 'draws' from \Omega ~ N(\hat{Omega}, I(\Omega)^{-1}) and f(\b|...; \Omega).
#' ####

# sourceCpp('helper.cpp')

# Draw from \Omega ~ N(\hat{Omega}, I(\Omega)^{-1}) -----------------------
Omega.draw <- function(fit){
  Omega.mean <- setNames(c(
    vech(fit$coeffs$D),
    c(fit$coeffs$beta),
    fit$coeffs$theta,
    fit$coeffs$gamma,
    fit$coeffs$eta
  ), names(fit$SE))
  
  Omega.vcov <- fit$vcov
  
  draw <- MASS::mvrnorm(1, Omega.mean, solve(Omega.vcov))
  
  Dinds <- grepl('^D\\[', names(draw))
  binds <- grepl('^beta\\_', names(draw))
  tinds <- grepl('^theta', names(draw))
  ginds <- grepl('gamma', names(draw))
  einds <- grepl('^cont|^bin', names(draw))
  
  D <- vech2mat(draw[Dinds], ncol(fit$RE))
  beta <- draw[binds]
  theta <- draw[tinds]
  gamma <- draw[ginds]
  eta <- draw[einds]
  
  return(list(
    D = D, beta = beta, theta = theta, gamma = gamma, eta = eta
  ))
}

# Draw b ------------------------------------------------------------------
b.draw <- function(b0, X, Y, Z, beta, theta, D, Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta){
  b.hat <- ucminf::ucminf(b0,
                          joint2_density, joint_density_ddb,
                          X, Y, Z, beta, theta, D, Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta)$par
  Sigmai <- solve(joint_density_sdb(b.hat, X, Y, Z, beta, theta, D, Delta, 
                                    K, Fi, l0i,
                                    KK, Fu, haz, gamma, eta, 1e-3))
  
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
  gamma <- Omega$gamma
  Fu <- surv$Fu.t
  b <- b
  
  exp(-l0u %*% exp(KK %*% eta + Fu %*% (gamma * b)))
}

