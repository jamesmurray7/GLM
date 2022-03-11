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
    fit$coeffs$eta,
    fit$coeffs$theta
  ), names(fit$SE))
  
  Omega.vcov <- fit$vcov
  
  draw <- MASS::mvrnorm(1, Omega.mean, solve(Omega.vcov))
  
  Dinds <- grepl('^D\\[', names(draw))
  binds <- grepl('^G\\_|^B\\_|^P\\_|^NB\\_', names(draw))
  sinds <- grepl('var.e', names(draw))
  ginds <- grepl('gamma', names(draw))
  einds <- grepl('^cont|^bin', names(draw))
  tinds <- grepl('^theta', names(draw))
  
  D <- vech2mat(draw[Dinds], ncol(fit$RE))
  beta <- draw[binds]
  var.e <- draw[sinds]
  gamma <- draw[ginds]
  eta <- draw[einds]
  theta <- draw[tinds]
  
  return(list(
    D = D, beta = beta, var.e = var.e, gamma = gamma, eta = eta, theta = theta
  ))
}

# Draw b ------------------------------------------------------------------
b.draw <- function(b0, X, Y, Z, beta, var.e, theta, D, Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta, Sigma0){
  
  # work out if count was fit by NB
  if(length(theta) == 0){
    nb <- F; theta <- 0
  }else{
    nb <- T; theta <- theta
  }
  
  # Create residual variance matrix, V
  V <- diag(x = var.e, nrow = nrow(Y), ncol = nrow(Y))
  b.hat <- ucminf::ucminf(b0,
                          joint_density, joint_density_ddb,
                          X, Z, beta, V, D, Y[, 1], Y[, 2], Y[, 3],
                          nb, theta, Delta, K, Fi, l0i, KK, Fu, haz, rep(gamma, each = 2), eta)$par
  
  # pp <- prod(pnorm(b.hat, b0, sqrt(diag(Sigma0)), lower.tail = F))
  # pp2 <- mvtnorm::pmvnorm(lower = -Inf, upper = b.hat, mean = b0, sigma = Sigma0, algorithm = mvtnorm::Miwa())
  # print(b0)
  # print(b.hat)
  # print(pp)
  # print(pp2)
  
  
  Sigmai <- solve(joint_density_sdb(b.hat, X, Z, beta, V, D, 
                                    Y[, 1], Y[, 2], Y[, 3], nb, theta, 
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

# Metropolis-Hastings scheme, this from joineRML and geared towards a N(m,S) RV
# Seems to be broken for non-Gaussian stuffs?
b.mh <- function(Omega.draw, Sigmai.prop, b.curr, pt){
  accept <- 0
  O <- Omega.draw
  b.prop <- MASS::mvrnorm(n = 1, mu = b.curr, Sigma = Sigmai.prop)
  log.a1 <- (-1 * joint_density(b.prop, pt$long$Xt, pt$long$Yt, pt$long$Zt,
                                O$beta, O$D, pt$surv$Delta, pt$surv$K,
                                pt$surv$Fi, pt$surv$l0i, pt$surv$KK.t, pt$surv$Fu.t,
                                pt$surv$l0u.t, O$gamma, O$eta)) - (-1 * joint_density(b.curr, pt$long$Xt, pt$long$Yt, pt$long$Zt,
                                                                                      O$beta, O$D, pt$surv$Delta, pt$surv$K,
                                                                                      pt$surv$Fi, pt$surv$l0i, pt$surv$KK.t, pt$surv$Fu.t,
                                                                                      pt$surv$l0u.t, O$gamma, O$eta))
  
  dens.curr <- mvtnorm::dmvnorm(x = b.curr, sigma = Sigmai.prop, log = T)
  dens.prop <- mvtnorm::dmvnorm(x = b.prop, sigma = Sigmai.prop, log = T)
  
  log.a2 <- dens.curr - dens.prop
  print(log.a1)
  print(dens.curr)
  print(dens.prop)
  a <- min(exp(log.a1 - log.a2), 1)
  randu <- runif(1)
  if (randu <= a) {
    b.curr <- b.prop
    accept <- 1
  }
  out <- list(b = b.curr, accept = accept)
  return(out)
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

