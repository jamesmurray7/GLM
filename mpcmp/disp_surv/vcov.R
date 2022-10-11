#' #################################################################
#' vcov.R
#' -------
#' Calculate the observed empirical information matrix for 
#' the psuedo-MLEs of a joint-model fit by approximate EM method.
#' This also calculates the variance of subject-specific \delta_i
#'           by (approximation of) the Hessian matrix.
#' #################################################################


vcov <- function(Omega, delta, dmats, surv, sv, Sigma, b, l0u, w, v, n, summax, inds.met, 
                 delta.update.quad, beta.update.quad){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  gamma <- c(Omega$gamma.surv)
  gamma.disp <- c(Omega$gamma.disp)
  zeta <- c(Omega$zeta)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- dmats$Z
  X <- dmats$X
  Y <- dmats$Y
  W <- dmats$W
  lY <- lapply(Y, lfactorial)
  m <- sapply(Y, length)
  
  #' Survival //
  S <- sv$S
  SS <- sv$SS
  Fi <- sv$Fi
  Fu <- sv$Fu
  Delta <- sv$Delta
  
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
  
  #' Calculate things we need for scores on \beta and \delta
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b, SIMPLIFY = F)
  nus <- mapply(function(W, delta) exp(W %*% delta), W = W, delta = delta, SIMPLIFY = F)
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)
  
  # Score for fixed effects (beta)
  if(beta.update.quad){
    Sb <- mapply(function(b, X, Z, W, Y, lY, delta, tau, summax){
      Sbeta2(beta, b, X, Z, W, Y, lY, delta, tau, w, v, summax)
    }, b = b, X = X, Z = Z, W = W, Y = Y, lY = lY, delta = delta, tau = tau, summax = summax, SIMPLIFY = F)
  }else{
    Sb <- mapply(function(b, X, Z, W, Y, lY, delta, summax){
      Sbeta_noquad(beta, b, X, Z, W, Y, lY, delta, summax)
    }, b = b, X = X, Z = Z, W = W, Y = Y, lY = lY, delta = delta, summax = summax, SIMPLIFY = F)
  }

  #' Survival parameters (\gamma, \zeta)
  inds.to.remove <- which(apply(do.call(rbind, delta), 1, function(i) any(abs(i) == max.val)))
  
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta, WFu, WFi, delta){
    Sgammazeta(c(gamma, gamma.disp, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, 
               WFu, WFi, delta, w, v, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, l0u = sv$l0u, Delta = sv$Delta,
  WFu = sv$WFu, WFi = sv$WFi, delta = delta, SIMPLIFY = F)
  
  Sgz2 <- lapply(1:n, function(i){
    if(i %in% inds.to.remove){
      S <- Sgz[[i]]
      S[2] <- 0           # Second element is gamma.disp
      return(S)
    }else{
      return(Sgz[[i]])
    }
  })
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Sgz){
    c(sD, c(Sb), c(Sgz))
  }, sD = sD, Sb = Sb, Sgz = Sgz2)

  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ Observed empirical information matrix (Mclachlan and Krishnan, 2008).
  
  # Information of subject-specific dispersion parameter
  Idelta <- lapply(1:n, function(i){
    if(i %in% inds.met){
      # Appraise the CMP part (which dominates in practice...)
      lls <- delta.update(delta[[i]], dmats$X[[i]], dmats$Z[[i]], dmats$W[[i]], dmats$Y[[i]],
                          b[[i]], beta, summax[[i]], w, v, tau[[i]], delta.update.quad)$Hessian
      # Contribution of presence of delta in survival log-likelihood...
      svs <- delta_update(delta[[i]], sv$WFu[[i]], sv$WFi[[i]], sv$Delta[[i]],
                          sv$SS[[i]], sv$Fu[[i]], b[[i]], gamma, gamma.disp, zeta,
                          sv$l0u[[i]], Sigma[[i]], w, v)$Hessian
      return(solve(-1*(-lls + svs)))
    }else{
      return(matrix(0, nr = dmats$w, nc = dmats$w))
    }
  })
  
  list(
    I = I,
    Idelta = Idelta
  )
}

