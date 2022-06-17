
vcov <- function(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n, N, summax){
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
  nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)
  mus2 <- lapply(mus, mu_fix, N)
  nus2 <- lapply(nus, mu_fix, N)
  
  #' Grid lookups ----
  #' Grid lookups ----
  lambdas <- mapply(function(mu, nu){
    m <- (mu*(N/10)) - 1; n <- (nu*(N/10)) - 1
    mat_lookup(m, n, lambda.mat)
  }, mu = mus2, nu = nus2, SIMPLIFY = F)
  
  Vs <- mapply(function(mu, nu){
    m <- (mu*(N/10)) - 1; n <- (nu*(N/10)) - 1
    mat_lookup(m, n, V.mat)
  }, mu = mus2, nu = nus2, SIMPLIFY = F)
  
  logZ <- mapply(function(mu, nu){
    m <- (mu*(N/10)) - 1; n <- (nu*(N/10)) - 1
    mat_lookup(m, n, logZ.mat)
  }, mu = mus2, nu = nus2, SIMPLIFY = F)
  ABC <- mapply(calc.ABC, mus2, nus2, lambdas, logZ, 100, SIMPLIFY=F)
  
  #' Score for the fixed effects, \beta 
  Sb <- mapply(Sbeta, X, Y, mus2, nus2, lambdas, Vs, SIMPLIFY = F)
  #' Score for the dispersion parameter(s), \delta 
  Sd <- mapply(function(ABC, Y, mu, V, nu, G){
    crossprod(((ABC$A * (Y - mu) / V - lgamma(Y + 1) + ABC$B) * nu), G)
  }, ABC = ABC, Y = Y, mu = mus2, V = Vs, nu = nus2, G = G, SIMPLIFY = F)

  #' Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Sd, Sgz){
    c(sD, c(Sb), c(Sd), Sgz)
  }, sD = sD, Sb = Sb, Sd = Sd, Sgz = Sgz)

  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  I
}
