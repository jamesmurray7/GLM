
vcov <- function(Omega, delta, dmats, surv, sv, Sigma, b, l0u, w, v, n, summax, inds.met, delta.update.quad){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  gamma <- c(Omega$gamma)
  zeta <- c(Omega$zeta)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- dmats$Z
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
  
  #' Calculate things we need for scores on \beta and \delta
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b, SIMPLIFY = F)
  nus <- mapply(function(delta, Y) rep(exp(delta), length(Y)), delta = delta, Y = Y, SIMPLIFY = F)
  
  #' lambda, logZ and V calculations ----
  lambdas <- mapply(function(mu, nu, summax){
    lambda_appx(mu, nu, summax)
  }, mu = mus, nu = nus, summax = summax, SIMPLIFY = F)
  
  logZs <- mapply(function(mu, nu, lambda, summax){
    logZ_c(log(lambda), nu, summax)
  }, mu = mus, nu = nus, lambda = lambdas, summax = summax, SIMPLIFY = F)
  
  Vs <- mapply(function(mu, nu, lambda, logZ, summax){
    calc_V_vec(mu, lambda, nu, logZ, summax)
  }, mu = mus, nu = nus, lambda = lambdas, logZ = logZs, summax = summax, SIMPLIFY = F)
  
  #' Score for the dispersion parameter(s), \delta 
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)
  
  #' Score for the fixed effects, \beta 
  # Sb <- mapply(Sbeta, X, Y, mus, nus, lambdas, Vs, SIMPLIFY = F)
  # Sb <- mapply(function(b, X, Y, Z, lY, delta, tau, mu , nu, lam, V, summax){
  #   a <- Sbeta(X, Y, mu, nu, lam, V)
  #   if(any(a > 1e3) | any(is.nan(a))){
  #     return(Sbeta_cdiff(beta, b, X, Z, Y, lY, delta, tau, w, v, summax, .Machine$double.eps^(1/3)))
  #   }else
  #     return(a)
  # }, b = b, X = X, Z = Z, Y = Y, lY = lY, delta = delta, tau = tau, 
  # mu = mus, nu = nus, lam = lambdas, V = Vs, summax = summax, SIMPLIFY = F)
  
  Sb <- mapply(function(b, X, Z, Y, lY, delta, tau, summax){
    Sbeta2(beta, b, X, Z, Y, lY, delta, tau, w, v, summax)
  }, b = b, X = X, Z = Z, Y = Y, lY = lY, delta = delta, tau = tau, summax = summax, SIMPLIFY = F)

  #' Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Sgz){
    c(sD, c(Sb), Sgz)
  }, sD = sD, Sb = Sb, Sgz = Sgz)

  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  
  
  # Variances of subject-specific delta with/out quadrature (NB: THIS)
  Idelta <- lapply(1:n, function(i){
    if(i %in% inds.met){
      a <- delta.update(delta[[i]], X[[i]], Z[[i]], Y[[i]], b[[i]], beta, 
                        summax[[i]], w, v, tau[[i]], delta.update.quad)
      out <- -1 * a$Hessian
    }else{
      out <- 0
    }
    out
  })
  
  list(
    I = I,
    Idelta = unlist(Idelta)
  )
}

