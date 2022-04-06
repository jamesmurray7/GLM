#' ###
#' Calculating vcov for EM fit using observed empirical information matrix
#' ###

vcov <- function(Omega, data.mat, b, Sigmai, l0u, gh.nodes, n, beta.quad){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  gamma <- c(Omega$gamma)
  eta <- c(Omega$eta)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- data.mat$Z
  X <- data.mat$X
  Y <- data.mat$Y
  m <- data.mat$m
  
  #' Survival //
  K <- data.mat$K
  KK <- data.mat$KK
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
  if(beta.quad){
    Sb <- mapply(function(X, Y, Z, b, S){
      Sbeta_quad(beta, X, Y, Z, b, S, w, v, 1e-4)
    }, X = X, Y = Y, Z = Z, b = b, S = Sigmai, SIMPLIFY = F)
  }else{
    Sb <- mapply(function(X, Y, Z, b, V){
      Sbeta(beta, X, Y, Z, b)
    }, X = X, Y = Y, Z = Z, b = b, SIMPLIFY = F)
  }
  #' (gamma, eta)
  Sge <- mapply(function(b, S, K, KK, Fu, Fi, l0u, Delta){
    Sgammaeta(c(gamma, eta), b, S, K, KK, Fu, Fi[1:2], l0u, Delta, w, v, 1e-4)
  }, b = b, S = Sigmai, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  
  # Collate and form information --------------------------------------------
  
  S <- mapply(function(sD, Sb, Sge){
    c(sD, c(Sb), Sge)
  }, sD = sD, Sb = Sb, Sge = Sge)
  
  
  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  I
}
  