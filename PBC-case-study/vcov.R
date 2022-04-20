
vcov <- function(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n, family){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  if(!is.null(Omega$sigma)) sigma <- Omega$sigma else sigma <- 0
  gamma <- c(Omega$gamma)
  zeta <- c(Omega$zeta)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- dmats$Z
  X <- dmats$X
  Y <- dmats$Y
  m <- dmats$m
  
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
  
  #' The fixed effects, \beta 
  Sb <- mapply(function(X, Y, Z, b){
    Sbeta(beta, X, Y, Z, b, sigma, family)
  }, X = X, Y = Y, Z = Z, b = b, SIMPLIFY = F)
  
  #' The dispersion parameter, \sigma
  if(family == 'gaussian'){
    tau <- mapply(function(Sigma, Z) sqrt(diag(tcrossprod(Z %*% Sigma, Z))), Z = Z, Sigma = Sigma, SIMPLIFY = F)
    Ss <- list()
    for(i in 1:n){
      rhs <- 0
      for(l in 1:length(w)){
        rhs <- rhs + w[l] * crossprod(Y[[i]] - X[[i]] %*% beta - Z[[i]] %*% b - v[l] * tau[[i]])
      } 
      Ss[[i]] <- -m[i]/(2*sigma) + 1/(2 * sigma^2) * rhs
    }
  }else if(family == 'negative.binomial'){
    Ss <- mapply(function(X, Y, Z, b, Sigma){
      Stheta(sigma, beta, X, Y, Z, b, Sigma, w, v, .Machine$double.eps^(1/3))
    }, X = X, Y = Y, Z = Z, b = b.hat, Sigma = Sigma, SIMPLIFY = F)
  }else{
    Ss <- NULL
  }
  
  #' Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = b.hat, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Ss, Sgz){
    c(sD, c(Sb), Ss, Sgz)
  }, sD = sD, Sb = Sb, Ss = Ss, Sgz = Sgz)
  
  
  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  I
}
