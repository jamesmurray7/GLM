#' Working out the observed empirical information matrix.
vcov <- function(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n, include.all,
                 inds.met, phi.update.quad, beta.update.quad){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  phi <- c(Omega$phi)
  gamma <- c(Omega$gamma)
  zeta <- c(Omega$zeta)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- dmats$Z
  W <- dmats$W
  X <- dmats$X
  Y <- dmats$Y
  lY <- lapply(Y, lfactorial)
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
  
  #' Score for the dispersion parameter(s), \delta and score for the fixed effects, \beta 
  phivec <- mapply(function(W) W %*% phi, W = W, SIMPLIFY = F)
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)
  
  if(beta.update.quad){
    Sb <- mapply(function(b, X, Z, W, Y, lY, tau, summax){
      Sbeta2(beta, b, X, Z, W, Y, lY, delta, tau, w, v, summax)
    }, b = b, X = X, Z = Z, W = W, Y = Y, lY = lY, tau = tau, summax = summax, SIMPLIFY = F)
  }else{
    Sb <- mapply(function(b, X, Y, Z, phivec){
      long_derivs(b = b, X = X, Y = Y, Z = Z, beta = beta, phi = phivec, design = X)$grad
    }, b = b, X = X, Y = Y, Z = Z, phivec = phivec, SIMPLIFY = F)
  }
  
  # Information of subject-specific dispersion parameter
  Sp <- lapply(1:n, function(i){ 
    if(include.all){
      return(phi_update(b[[i]], X[[i]], Y[[i]], Z[[i]], W[[i]], beta, phi,
                        w=w, v=v, tau[[i]])$Score)
    }else if(i %in% inds.met){
      return(phi_update(b[[i]], X[[i]], Y[[i]], Z[[i]], W[[i]], beta, phi,
                        w, v, tau[[i]])$Score)
    }else{
      return(rep(0, dmats$w)) # No contribution.
    }
  })
  
  #' Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Sp, Sgz){
    c(sD, c(Sb), c(Sp), Sgz)
  }, sD = sD, Sb = Sb, Sp = Sp, Sgz = Sgz)

  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  
  I
}
