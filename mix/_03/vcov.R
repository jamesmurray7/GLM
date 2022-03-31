#' ####
#' Calculating the observed emprical information matrix (-ve Hessian)
#' ####

vcov <- function(Omega, data.mat, V, b, bsplit, bmat, Sigmai, SigmaiSplit, l0u, gh.nodes, nb, n, quad){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  var.e <- Omega$var.e
  gamma <- c(Omega$gamma)
  eta <- c(Omega$eta)
  if(nb) theta <- Omega$theta else theta <- 0.0
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- data.mat$Z; Zk <- data.mat$Zk
  X <- data.mat$X; Xk <- data.mat$Xk
  Y <- data.mat$Y; Yk <- data.mat$Yk
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
  if(quad){
    Sb <- mapply(function(X, Y, Z, b, V, S){
      Sbeta_quad(beta, X, Y[,1], Y[,2], Y[,3], Z, b, V, S[[2]], w, v, nb, theta, eps = 1e-4)
    }, X = X, Y = Y, Z = Z, b = b, V = V, S = SigmaiSplit, SIMPLIFY = F)
  }else{
    Sb <- mapply(function(X, Y, Z, b, V){
      Sbeta(beta, X, Y[,1], Y[,2], Y[,3], Z, b, V, nb, theta)
    }, X = X, Y = Y, Z = Z, b = b, V = V, SIMPLIFY = F)
  }
  
  #' Residual variance for the Gaussian response
  tau.long <- mapply(function(Z, S){
    sqrt(diag(Z %*% S[[1]] %*% t(Z))) # just the first block is associated w/ Gaussian response
  }, Z = Z, S = SigmaiSplit, SIMPLIFY = F)
  
  Ss <- mapply(function(X, Y, Z, b, m, tau){
    rhs <- 0
    for(l in 1:gh.nodes) rhs <- rhs + w[l] * crossprod(Y[,1] - X %*% beta[1:4] - Z %*% b[[1]] - v[l] * tau)
    -m/(2*var.e) + 1/(2 * var.e^2) * rhs
  }, X = X, Y = Y, Z = Z, b = bsplit, m = m, tau = tau.long, SIMPLIFY = F)
  
  #' Dispersion parameter, theta
  if(nb){
    St <- mapply(function(b, X, Y, Z, S){
      Stheta(theta, beta[9:12], X, Y[, 3], Z, b[[3]], S[[3]], w, v, 1e-4)
    }, b = bsplit, X = X, Y = Y, Z = Z, S = SigmaiSplit)
  }else{
    St <- NULL
  }
  
  #' (gamma, eta)
  Sge <- mapply(function(b, S, K, KK, Fu, Fi, l0u, Delta){
    Sgammaeta(c(gamma, eta), b, S, K, KK, Fu, Fi[1:2], l0u, Delta, w, v, 1e-4)
  }, b = bmat, S = SigmaiSplit, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)


  # Collate and form information --------------------------------------------
  
  S <- mapply(function(sD, Sb, Ss, Sge){
    c(sD, c(Sb), Ss, Sge)
  }, sD = sD, Sb = Sb, Ss = Ss, Sge = Sge)

  if(nb) S <- rbind(S, St)
  
  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  I
}
