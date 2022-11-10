#' Working out the observed empirical information matrix.
vcov <- function(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  shape <- c(Omega$shape)
  gamma <- c(Omega$gamma)
  zeta <- c(Omega$zeta)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- dmats$Z
  X <- dmats$X
  Y <- dmats$Y
  m <- sapply(Y, length)
  
  #' Survival //
  S <- sv$S
  SS <- sv$SS
  Fi <- sv$Fi
  Fu <- sv$Fu
  Delta <- surv$Delta
  
  # Scores ------------------------------------------------------------------
  #' The RE covariance matrix, D
  Dinv <- solve(D)
  vech.indices <- which(lower.tri(D, diag = T), arr.ind = T)
  dimnames(vech.indices) <- NULL
  delta.D <- lapply(1:nrow(vech.indices), function(d){
    out <- matrix(0, nrow(D), ncol(D))
    ind <- vech.indices[d, 2:1]
    out[ind[1], ind[2]] <- out[ind[2], ind[1]] <- 1 # dD/dvech(d)_i
    out
  })
  
  lhs <- sapply(delta.D, function(d) {
    -0.5 * sum(diag(Dinv %*% d))
  })
  
  sDi <- function(i) {
    mapply(function(b) {
      out <- 0.5 * tcrossprod(b) %*% (Dinv %*% delta.D[[i]] %*% Dinv)   
      lhs[i] + sum(diag(out))
    },
    b = b,
    SIMPLIFY = T)
  }
  
  sD <- sapply(1:nrow(vech.indices), sDi)
  sD <- lapply(1:nrow(sD), function(x) sD[x, ]) # Cast to list
  
  #' Score for the shape parameter and score for the fixed effects, \beta 
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)
  
  beta.update <- mapply(function(b, X, Y, Z){
    long_derivs(b = b, X = X, Y = Y, Z = Z, beta = beta, shape = shape, design = X)
  }, b = b, X = X, Y = Y, Z = Z, SIMPLIFY = F)
  Sb <- lapply(beta.update, el, 1)
  
  shape.update <- mapply(function(b, X, Y, Z, tau){
    pracma::grad(E_shape.b, shape, X = X, Y = Y, Z = Z, tau = tau, beta = beta, b = b, w = w, v = v)
  }, b = b, X = X, Y = Y, Z = Z, tau = tau, SIMPLIFY = F)
  Ss <- lapply(shape.update, el, 1)
  
  #' Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Ss, Sgz){
    c(sD, c(Sb), c(Ss), Sgz)
  }, sD = sD, Sb = Sb, Ss = Ss, Sgz = Sgz)

  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  
  I
}
