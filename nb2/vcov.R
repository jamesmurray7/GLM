#' ###
#' Calculating vcov for EM fit using observed empirical information matrix
#' (negbin fit)
#' ###

vcov <- function(Omega, data.mat, b, Sigmai, l0u, gh.nodes, n){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  gamma <- c(Omega$gamma)
  zeta <- c(Omega$zeta)
  alpha <- c(Omega$alpha)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- data.mat$Z
  X <- data.mat$X
  Y <- data.mat$Y
  W <- data.mat$W
  m <- data.mat$m
  
  #' Survival //
  S <- data.mat$S
  SS <- data.mat$SS
  Fi <- data.mat$Fi
  Fu <- data.mat$Fu
  Delta <- data.mat$Delta
  
  gh <- statmod::gauss.quad.prob(gh.nodes, 'normal')
  w <- gh$w; v <- gh$n
  
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
      out <- 0.5 * tcrossprod(b) %*% (Dinv %*% delta.D[[i]] %*% Dinv)   # (Sigmai + tcrossprod(b))?
      lhs[i] + sum(diag(out))
    },
    b = b,
    SIMPLIFY = T)
  }
  
  sD <- sapply(1:nrow(vech.indices), sDi)
  sD <- lapply(1:nrow(sD), function(x) sD[x, ]) # Cast to list
  
  #' beta
  theta <- mapply(function(W) exp(W %*% alpha), W = W, SIMPLIFY = F)
  Sb <- mapply(function(X, Y, Z, b, theta){
    crossprod(X, Score_eta(b, X, Y, Z, beta, theta))
  }, X = X, Y = Y, Z = Z, b = b, theta = theta, SIMPLIFY = F)
  
  #' (gamma, eta)
  Sge <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammaeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi[1:2], l0u, Delta, w, v, 1e-4)
  }, b = b, Sigma = Sigmai, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  #' alpha
  Sa <- mapply(function(b, X, Y, Z, W, Sigma){
    Salpha(alpha, beta, X, Y, Z, W, b, Sigma, w, v, 1e-4)
  }, b = b, X = X, Y = Y, Z = Z, W = W, Sigma = Sigmai, SIMPLIFY = F)
  
  
  # Collate and form information --------------------------------------------
  
  S <- mapply(function(sD, Sb, Sa, Sge){
    c(sD, c(Sb), Sa, Sge)
  }, sD = sD, Sb = Sb, Sa = Sa, Sge = Sge)
  
  
  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  I
}
  