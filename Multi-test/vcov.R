#' Mclachlan and Krishnan (2008) approximation of the information matrix.
vcov <- function(Omega, dmats, surv, sv, 
                 Sigma, SigmaSplit, b, bsplit, 
                 l0u, w, v, n, family, K, q, beta.inds, b.inds){
  #' Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  sigma <- Omega$sigma
  gamma <- c(Omega$gamma)
  zeta <- c(Omega$zeta)
  
  #' Extract data objects ----
  #' Longitudinal //
  Z <- dmats$Z
  X <- dmats$X
  Y <- dmats$Y
  m <- lapply(Y, function(y) sapply(y, length))
  
  #' Survival //
  S <- sv$S
  SS <- sv$SS
  Fi <- sv$Fi
  Fu <- sv$Fu
  Delta <- surv$Delta
  
  beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1) 
  
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
    Sbeta(beta, X, Y, Z, b, sigma, family, beta.inds2, K)
  }, X = X, Y = Y, Z = Z, b = bsplit, SIMPLIFY = F)
  
  Ss <- vector('list', K)
  #' The dispersion parameter, \sigma
  if(any(unlist(family) == 'gaussian')){
    gauss.inds <- which(unlist(family) == 'gaussian')
    for(j in seq_along(gauss.inds)){
      temp <- numeric(n)
      gj <- gauss.inds[j]
      tau <- mapply(function(S, Z) sqrt(diag(tcrossprod(Z[[gj]] %*% S[[gj]],
                                                        Z[[gj]]))), Z = Z, S = SigmaSplit, SIMPLIFY = F)
      for(i in 1:n){
        rhs <- 0
        for(l in 1:length(w)){
          rhs <- rhs + w[l] * crossprod(Y[[i]][[gj]] - X[[i]][[gj]] %*% beta[beta.inds[[gj]]] - Z[[i]][[gj]] %*% bsplit[[i]][[gj]] - v[l] * tau[[i]])
        }
        temp[i] <- -m[[i]][gj]/(2 * unlist(sigma)[gj]) + 1/(2 * unlist(sigma)[gj]^2) * rhs
      }
      Ss[[gj]] <- temp
    }
  }
  #' Negative binomial case
  if(any(unlist(family) == 'negative.binomial')){
    nb.inds <- which(unlist(family) == 'negative.binomial')
    for(j in seq_along(nb.inds)){
      nj <- nb.inds[j]
      Ss[[nj]] <- mapply(function(X, Y, Z, b, Sigma){
        Stheta(sigma[[nj]], beta[beta.inds[[nj]]], X[[nj]], Y[[nj]], Z[[nj]], b[[nj]], Sigma[[nj]], w, v, .Machine$double.eps^(1/3))
      }, X = X, Y = Y, Z = Z, b = bsplit, Sigma = Sigmasplit, SIMPLIFY = T)
    }
  }

  #' Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, q, .Machine$double.eps^(1/3))
  }, b = b, Sigma = SigmaSplit, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, q, .Machine$double.eps^(1/3))
  }, b = b, Sigma = SigmaSplit, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  Hgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Hgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, q, .Machine$double.eps^(1/3))
  }, b = b, Sigma = SigmaSplit, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  HH <<- Reduce('+', Hgz) 
  print(solve(-HH))
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Sgz){
    c(sD, c(Sb), Sgz)
  }, sD = sD, Sb = Sb, Sgz = Sgz)

  if(any(unlist(family) %in% c('gaussian', 'negative.binomial'))) S <- rbind(S, do.call(rbind, Ss))
  
  SS <- rowSums(S) # sum S
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  # ^ observed empirical information matrix (Mclachlan and Krishnan, 2008).
  I
}

  