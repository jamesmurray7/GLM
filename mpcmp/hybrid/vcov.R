
vcov <- function(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, num, summax,
                 all.mus, nu.vec, lambda.mat, logZ.mat, V.mat){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  delta <- c(Omega$delta)
  gamma <- c(Omega$gamma)
  zeta <- c(Omega$zeta)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- dmats$Z
  G <- dmats$G
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
  
  #' NEW \mus, \nus, and calculate \tau
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b, SIMPLIFY = F)
  nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)
  
  #' Indices for lookup given new \b.hat
  m <- mapply(function(a) get_indices(a, all.mus), a = mus, SIMPLIFY = F)
  n <- mapply(function(a) get_indices(a, round(nu.vec, 3)), a = nus, SIMPLIFY = F)
  
  #' \lambdas, Vs for \beta update. (With hardcode for NA/out-of-range values).
  lambdas <- mapply(function(mu, nu, m, n){
    out <- numeric(length(mu))
    if(any(is.na(m))){
      nas <- is.na(m)
      out[nas] <- lambda_appx(mu[nas], nu[nas], summax)
      out[!nas] <- mat_lookup(m[!nas], n[!nas], lambda.mat)
    }else{
      out <- mat_lookup(m, n, lambda.mat)
    }
    out
  }, mu = mus, nu = nus, m = m, n = n, SIMPLIFY = F)
  
  Vs <- mapply(function(mu, nu, m, n){
    out <- numeric(length(mu))
    if(any(is.na(m))){
      nas <- is.na(m)
      lams <- lambda_appx(mu[nas], nu[nas], summax)
      logZs <- logZ_c(log(lams), nu[nas], summax)
      out[nas] <- calc_V_vec(mu[nas], lams, nu[nas], logZs, summax)
      out[!nas] <- mat_lookup(m[!nas], n[!nas], V.mat)
    }else{
      out <- mat_lookup(m, n, V.mat)
    }
    out
  }, mu = mus, nu = nus, m = m, n = n, SIMPLIFY = F)
  
  #' Score for the fixed effects, \beta 
  Sb <- mapply(Sbeta, X, Y, mus, nus, lambdas, Vs, SIMPLIFY = F)
  
  #' Score for the dispersion parameter(s), \delta 
  Sd <- mapply(function(G, b, X, Z, Y, lY, tau){
    Sdelta_cdiff(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps = .Machine$double.eps^(1/3))
  }, G = G, b = b, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)

  #' Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Sd, Sgz){
    c(sD, c(Sb), c(Sd), Sgz)
  }, sD = sD, Sb = Sb, Sd = Sd, Sgz = Sgz)

  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:num, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/num
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  I
}
