#' ####
#' Calculating the observed emprical information matrix (-ve Hessian)
#' ####

vcov <- function(Omega, data.mat, V, b, bsplit, bmat, Sigmai, SigmaiSplit, l0u, gh.nodes, n, quad){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  var.e <- Omega$var.e
  gamma <- c(Omega$gamma)
  gamma_r <- rep(gamma, c(2, 1, 2))
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- data.mat$Z; X <- data.mat$X; Y <- data.mat$Y; 
  m <- data.mat$m
  
  #' Survival //
  Fi <- data.mat$Fi
  Fu <- data.mat$Fu
  Fu.list <- data.mat$Fu.list
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
 
  #' \beta
  if(quad){
    Sb <- mapply(function(X, Y, Z, b, V, S){
      Sbeta_quad(beta, X, Y[,1], Y[,2], Y[,3], Z, b, V, S[[2]], w, v, .Machine$double.eps^(1/3))
    }, X = X, Y = Y, Z = Z, b = bsplit, V = V, S = SigmaiSplit, SIMPLIFY = F)
  }else{
    Sb <- mapply(function(X, Y, Z, b, V){
      Sbeta(beta, X, Y[,1], Y[,2], Y[,3], Z, b, V)
    }, X = X, Y = Y, Z = Z, b = bsplit, V = V, SIMPLIFY = F)
  }
  
  #' Residual variance for the Gaussian response
  tau.long <- mapply(function(Z, S){
    sqrt(diag(Z[["gc"]] %*% S[[1]] %*% t(Z[["gc"]]))) # just the first block is associated w/ Gaussian response
  }, Z = Z, S = SigmaiSplit, SIMPLIFY = F)
  
  Ss <- mapply(function(X, Y, Z, b, m, tau){
    rhs <- 0
    for(l in 1:gh.nodes) rhs <- rhs + w[l] * crossprod(Y[,1] - X %*% beta[1:4] - Z[["gc"]] %*% b[[1]] - v[l] * tau)
    -m/(2*var.e) + 1/(2 * var.e^2) * rhs
  }, X = X, Y = Y, Z = Z, b = bsplit, m = m, tau = tau.long, SIMPLIFY = F)
  
  #' \gamma
  Sg <- mapply(function(b, S, Fu, Fu.list, Fi, l0u, Delta){
    Sgamma(gamma, b, S, Fu[, 1:2, drop = F], Fu.list, Fi[1:2], l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = bmat, S = SigmaiSplit, Fu = Fu, Fu.list = Fu.list, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)


  # Collate and form information --------------------------------------------
  
  S <- mapply(function(sD, Sb, Ss, Sg){
    c(sD, c(Sb), Ss, Sg)
  }, sD = sD, Sb = Sb, Ss = Ss, Sg = Sg)
  
  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  I
}
