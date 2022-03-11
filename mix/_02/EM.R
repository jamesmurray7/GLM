#' ####    MIXTURE/
#' Approximate EM for joint model of negative-binomial (gamma-poisson deviate) x survival response.
#' ####
library(dplyr)
library(survival)
library(glmmTMB)
library(Rcpp)
library(RcppArmadillo)
source('simData.R')
source('inits.R')
source('survFns.R')
source('vcov.R')
sourceCpp('mix.cpp')
vech <- function(x) x[lower.tri(x, diag = T)]

EMupdate <- function(b, Y, X, Z, V,
                     D, beta, var.e, gamma, eta, 
                     Delta, l0i, l0u, Fi, Fu, K, KK, survdata, sv, w, v){
  n <- length(b)

  #' ### ------
  #' E step
  #' ### ------
  
  b.hat <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, haz){
    ucminf::ucminf(b, joint_density, joint_density_ddb,
                   X = X, Z = Z, beta = beta, vars = var.e, D = D,
                   Y_1 = Y[, 1], Y_2 = Y[, 2], Y_3 = Y[, 3],
                   Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu,
                   haz = haz, gamma = rep(gamma, each = 2), eta = eta)$par
  }, b = b, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
  KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)
  # split into consistuent parts in two ways
  bmat <- lapply(b.hat, matrix, nc = 2, byr = T)
  bsplit <- lapply(b.hat, function(y) lapply(split(seq(6), rep(seq(3), each = 2)), function(x) y[x]))
  
  Sigmai <- mapply(function(b, X, Z, Y, Delta, K, Fi, l0i, KK, Fu, l0u){
    solve(joint_density_sdb(b = b, X = X, Z = Z, beta = beta, vars = var.e, D = D,
                            Y_1 = Y[,1], Y_2 = Y[,2], Y_3 = Y[,3],
                            Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu, 
                            haz = l0u, gamma = rep(gamma, each = 2), eta = eta, eps = 1e-3))
  }, b = b.hat, X = X, Z = Z, Y = Y, Delta = Delta, K = K, Fi = Fi,
  l0i= l0i, KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)
  # Split out into constituent block-diag pieces.
  SigmaiSplit <- lapply(Sigmai, function(y) lapply(split(seq(6), rep(seq(3), each = 2)), function(x) y[x,x]))
  
  # Update to D
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, S = Sigmai, b = b.hat, SIMPLIFY = F)
  
  # Update to beta
  Sb <- mapply(function(X, Y, Z, b){
    Sbeta(beta, X, Y[,1], Y[,2], Y[,3], Z, b, var.e)
  }, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)
  
  Hb <- mapply(function(X, Y, Z, b){
    Hbeta(beta, X, Y[,1], Y[,2], Y[,3], Z, b, var.e, 1e-4)
  }, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)
  
  # Update to var.e (Gaussian submodels)
  tauL <- mapply(function(S, Z){
    out <- vector('list', 2)
    for(i in c(1,3)) out[[i]] <- sqrt(diag(tcrossprod(Z %*% S[[i]], Z)))
    out
  }, S = SigmaiSplit, Z = Z, SIMPLIFY = F)
  
  var.e.update <- mapply(function(b, Y, X, Z, tau){
    out <- numeric(3)
    for(l in 1:length(w)){
      out[1] <- out[1] + w[l] * crossprod(Y[,1] - X %*% beta[1:4]  - Z %*% b[1:2] - tau[[1]] * v[l])
      out[2] <- out[2] + w[l] * crossprod(Y[,3] - X %*% beta[9:12] - Z %*% b[5:6] - tau[[3]] * v[l])
    }
    out[3] <- length(Y[,1])
    out
  }, b = b.hat, Y = Y, X = X, Z = Z, tau = tauL)

  # gamma, eta
  Sge <- mapply(function(b, S, K, KK, Fu, Fi, l0u, Delta){
    Sgammaeta(c(gamma, eta), b, S, K, KK, Fu, Fi[1:2], l0u, Delta, w, v, 1e-4)
  }, b = bmat, S = SigmaiSplit, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)
  
  Hge <- mapply(function(b, S, K, KK, Fu, Fi, l0u, Delta){
    Hgammaeta(c(gamma, eta), b, S, K, KK, Fu, Fi[1:2], l0u, Delta, w, v, 1e-4)
  }, b = bmat, S = SigmaiSplit, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  #' ### ------
  #' M step
  #' ### ------
  
  # D
  D.new <- Reduce('+', Drhs)/n
  # beta
  beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
  # var.e
  var.e.update <- rowSums(var.e.update)
  var.e.new <- c(var.e.update[1]/var.e.update[3], var.e.update[2]/var.e.update[3])
  # (gamma, eta)
  gammaeta.new <- c(gamma, eta) - solve(Reduce('+', Hge), rowSums(Sge))
  # baseline hazard
  lambda <- lambdaUpdate(sv$surv.times, sv$ft, gamma, eta, K, SigmaiSplit, bsplit, w, v)
  # Baseline hazard objects
  l0.new <- sv$nev/rowSums(lambda)
  l0u.new <- lapply(l0u, function(x){
    ll <- length(x); l0.new[1:ll]
  })
  l0i.new <- c()
  l0i.new[which(unlist(Delta) == 0)] <- 0 
  l0i.new[which(unlist(Delta) == 1)] <- l0.new[match(survdata[which(unlist(Delta)==1), 'survtime'], sv$ft)]
  l0i.new <- as.list(l0i.new)
  
  return(list(
    D.new = D.new, beta.new = beta.new, var.e.new = var.e.new, 
    gamma.new = gammaeta.new[1:3], eta.new = gammaeta.new[4:5], 
    b.hat = b.hat, l0.new = l0.new, l0i.new = l0i.new, l0u.new = l0u.new
  ))
}

EM <- function(data, ph, survdata, gh = 9, tol = 0.01, verbose = F, post.process = T){
  start <- proc.time()[3]
  n <- nrow(survdata)
  # Get data matrices
  m <- Y <- X <- Z <- K <- list()
  for(i in 1:n){
    i.dat <- data[data$id == i, ]
    m[[i]] <- nrow(i.dat)
    Y[[i]] <- cbind(i.dat$Y.1, i.dat$Y.2, i.dat$Y.3)
    X[[i]] <- model.matrix(~time+cont+bin, i.dat)
    Z[[i]] <- model.matrix(~time, i.dat)
    K[[i]] <- cbind(unique(i.dat$cont), unique(i.dat$bin))
  }
  
  # initial conditions
  inits.long <- Longit.inits(data)
  b <- Ranefs(inits.long)
  beta <- inits.long$beta.init
  var.e <- inits.long$var.e.init
  D <- inits.long$D.init
  inits.surv <- TimeVarCox(data, b)
  b <- lapply(1:n, function(i) b[i, ])
  rm(inits.long) # large object
  
  # Populate off-block-diagonal elements (i.e. between responses) of D?
  # I think comment this out and see if it speeds up convergence / reduces number iterations?
  b.inds <- split(seq(1:6), rep(1:3, each = 2))
  temp <- cov(do.call(rbind, b))
  DD <- D
  for(i in 1:3){
    for(j in 1:3){
      if(i == j){
        break
      }else{
        bi <- b.inds[[i]]; bj <- b.inds[[j]]
        DD[bi, bj] <- temp[bi, bj]
        DD[bj, bi] <- t(temp[bi, bj])
      }
    }
  }

  if(any(eigen(DD)$values < 0) || det(DD) <= 0) DD <- as.matrix(Matrix::nearPD(DD, keepDiag = T)$mat)
  D <- DD
  
  
  # Survival objects
  sv <- surv.mod(ph, data, inits.surv$l0.init)
  Delta <- as.list(sv$Di)
  l0i <- as.list(sv$l0i)
  l0u <- sv$l0u
  Fi <- lapply(1:n, function(i) do.call(c, replicate(3, sv$Fi[i, ], simplify = F)))
  Fu <- sv$Fu
  KK <- sapply(1:n, function(i){
    x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
    if('numeric'%in%class(x)) x <- t(as.matrix(x))
    x
  })
  gamma <- inits.surv$inits[3:5]
  eta <- inits.surv$inits[1:2]
  
  # Quadrature //
  aa <- statmod::gauss.quad.prob(gh, 'normal')
  w <- aa$w; v <- aa$n
  
  # Collect parameters
  vD <- vech(D);
  names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
  params <- c(vD, beta, var.e, gamma, eta)
  
  # Collect data objects
  data.mat <- list(Y = Y, X = X, Z = Z, m = m,                          # Longit.
                   K = K, KK = KK, Fi = Fi, Fu = Fu, Deltai = Delta)   # Survival

  EMstart <- proc.time()[3]
  diff <- 100; iter <- 0
  while(diff > tol){
    update <- EMupdate(b, Y, X, Z, V, D, beta, var.e, gamma, eta, Delta, l0i, l0u, Fi, Fu, K, KK, survdata, sv, w, v)
    params.new <- c(vech(update$D.new), update$beta.new, update$var.e.new, update$gamma.new, update$eta.new)
    names(params.new) <- names(params)
    if(verbose) print(sapply(params.new, round, 4))
    diff <- max(
      abs(params.new-params)/(abs(params)+1e-3)
    )
    message('\nIteration ', iter + 1)
    message('Maximum relative difference: ', round(diff, 5))
    b <- update$b.hat
    D <- update$D.new
    beta <- update$beta.new
    var.e <- update$var.e.new
    gamma <- update$gamma.new
    eta <- update$eta.new
    l0 <- update$l0.new; l0u <- update$l0u.new; l0i <- update$l0i.new
    params <- params.new
    iter <- iter + 1
  }
  EMend <- proc.time()[3]
  coeffs <- list(D = D, beta = beta, var.e = var.e, gamma = gamma, eta = eta)
  out <- list(coeffs = coeffs,
              RE = do.call(rbind, b),
           #   inits = inits.long,
              EMtime = EMend-EMstart,
              iter = iter,
              totaltime = proc.time()[3] - start)
  out$hazard <- cbind(ft = sv$ft, haz = l0)
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time <- proc.time()[3]
    #' hat{b} at MLEs
    b <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, haz){
      ucminf::ucminf(b, joint_density, joint_density_ddb,
                     X = X, Z = Z, beta = beta, V = V, D = D,
                     Y_1 = Y[, 1], Y_2 = Y[, 2], Y_3 = Y[, 3], nb, theta,
                     Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu,
                     haz = haz, gamma = rep(gamma, each = 2), eta = eta)$par
    }, b = b, X = X, Y = Y, Z = Z, V = V, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
    KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)
    bmat <- lapply(b, matrix, nc = 2, byr = T)
    bsplit <- lapply(b, function(y) lapply(split(seq(6), rep(seq(3), each = 2)), function(x) y[x]))
    # And its covariance matrix
    Sigmai <- mapply(function(b, X, Z, V, Y, Delta, K, Fi, l0i, KK, Fu, l0u){
      solve(joint_density_sdb(b = b, X = X, Z = Z, beta = beta, V = V, D = D,
                              Y_1 = Y[,1], Y_2 = Y[,2], Y_3 = Y[,3], nb, theta,
                              Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu, 
                              haz = l0u, gamma = rep(gamma, each = 2), eta = eta, eps = 1e-3))
    }, b = b, X = X, Z = Z, V = V, Y = Y, Delta = Delta, K = K, Fi = Fi,
    l0i= l0i, KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)
    # Split out into constituent block-diag pieces.
    SigmaiSplit <- lapply(Sigmai, function(y) lapply(split(seq(6), rep(seq(3), each = 2)), function(x) y[x,x]))
    
    # Calculate Standard Errors
    I <- structure(vcov(coeffs, data.mat, V, b, bsplit, bmat, Sigmai, SigmaiSplit, l0u, gh, nb, n),
                   dimnames = list(names(params), names(params)))
    out$SE <- sqrt(diag(solve(I)))
    out$vcov <- I
    out$RE <- do.call(rbind, b)
    out$postprocess.time <- round(proc.time()[3]-start.time, 2)
  }
  out
}