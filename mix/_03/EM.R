#' ####    MIXTURE/
#' Approximate EM for mixture of joint models
#' This directory an implementation of Rustand et al. Simulation (ArXiv 220306256).
#' -- Also trying some new things wrt data set-up now b unequal.
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
                     D, beta, var.e, gamma,
                     Delta, l0i, l0u, Fi, Fu, Fu.list, survdata, sv, w, v, inds, quad){
  n <- length(b)
  gamma_r <- rep(gamma, c(2, 1, 2))
  
  #' ### ------
  #' E step
  #' ### ------
  
  b.hat <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, haz){
    ucminf::ucminf(b, joint_density, joint_density_ddb,
                   X = X, Z = Z, beta = beta, V = V, D = D,
                   Y_1 = Y[, 1], Y_2 = Y[, 2], Y_3 = Y[, 3], 
                   Delta = Delta, Fi = Fi, l0i = l0i, Fu = Fu,
                   haz = haz, gamma = gamma_r)$par
  }, b = b, X = X, Y = Y, Z = Z, V = V, Delta = Delta, Fi = Fi, l0i = l0i, Fu = Fu, haz = l0u, SIMPLIFY = F)
  bmat <- lapply(b.hat, function(x){
    matrix(c(x[1:3], 0, x[4:5]), nc = 2, by = T)
  })
  bsplit <- lapply(b.hat, function(y) lapply(inds, function(x) y[x + 1]))
  Sigmai <- mapply(function(b, X, Z, V, Y, Delta, Fi, l0i, Fu, l0u){
    solve(joint_density_sdb(b = b, X = X, Z = Z, beta = beta, V = V, D = D,
                            Y_1 = Y[,1], Y_2 = Y[,2], Y_3 = Y[,3],
                            Delta = Delta, Fi = Fi, l0i = l0i, Fu = Fu, 
                            haz = l0u, gamma = gamma_r, eps = 1e-3))
  }, b = b.hat, X = X, Z = Z, V = V, Y = Y, Delta = Delta, Fi = Fi,
  l0i= l0i, Fu = Fu, l0u = l0u, SIMPLIFY = F)
  # Split out into constituent block-diag pieces.
  SigmaiSplit <- lapply(Sigmai, function(y) lapply(inds, function(x) as.matrix(y[x + 1, x + 1])))
  
  # Update to D
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, S = Sigmai, b = b.hat, SIMPLIFY = F)
  
  #' \beta update
  if(quad){
    Sb <- mapply(function(X, Y, Z, b, V, S){
      Sbeta_quad(beta, X, Y[, 1], Y[, 2], Y[, 3], Z, b, V, S[[2]], w, v, .Machine$double.eps^(1/3))
    }, X = X, Y = Y, Z = Z, b = bsplit, V = V, S = SigmaiSplit)
    Hb <- mapply(function(X, Y, Z, b, V, S){
      Hbeta_quad(beta, X, Y[, 1], Y[, 2], Y[, 3], Z, b, V, S[[2]], w, v, .Machine$double.eps^(1/3))
    }, X = X, Y = Y, Z = Z, b = bsplit, V = V, S = SigmaiSplit, SIMPLIFY = F)
  }else{
    # Without quadrature on binomial leg.
    Sb <- mapply(function(X, Y, Z, b, V){
      Sbeta(beta, X, Y[, 1], Y[, 2], Y[, 3], Z, b, V)
    }, X = X, Y = Y, Z = Z, b = bsplit, V = V)
    Hb <- mapply(function(X, Y, Z, b, V){
      Hbeta(beta, X, Y[, 1], Y[, 2], Y[, 3], Z, b, V, .Machine$double.eps^(1/3))
    }, X = X, Y = Y, Z = Z, b = bsplit, V = V, SIMPLIFY = F)
  }
  
  # \sigma^2 (Gaussian submodel)
  tauL <- mapply(function(S, Z){
    unname(sqrt(diag(tcrossprod(Z[["gc"]] %*% S[[1]], Z[["gc"]]))))
  }, S = SigmaiSplit, Z = Z, SIMPLIFY = F)
  
  var.e.update <- mapply(function(b, Y, X, Z, tau){
    out <- numeric(2)
    for(l in 1:length(w)) out[1] <- out[1] + w[l] * crossprod(Y[,1] - X %*% beta[1:4] - Z[["gc"]] %*% b[[1]] - tau * v[l])
    out[2] <- length(Y[,1])
    out
  }, b = bsplit, Y = Y, X = X, Z = Z, tau = tauL)

  #' \gamma 
  Sg <- mapply(function(b, S, Fu, Fu.list, Fi, l0u, Delta){
    Sgamma(gamma, b, S, Fu[, 1:2, drop = F], Fu.list, Fi[1:2], l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = bmat, S = SigmaiSplit, Fu = Fu, Fu.list = Fu.list, Fi = Fi, l0u = l0u, Delta = Delta)
  
  Hg <- mapply(function(b, S, Fu, Fu.list, Fi, l0u, Delta){
    Hgamma(gamma, b, S, Fu[, 1:2, drop = F], Fu.list, Fi[1:2], l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = bmat, S = SigmaiSplit, Fu = Fu, Fu.list = Fu.list, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  #' ### ------
  #' M step
  #' ### ------
  
  # D
  D.new <- Reduce('+', Drhs)/n
  # beta
  beta.new <- beta - solve(Reduce('+', Hb), rowSums(Sb))
  # var.e
  var.e.update <- rowSums(var.e.update)
  var.e.new <- var.e.update[1]/var.e.update[2]
  # gamma
  gamma.new <- c(gamma) - solve(Reduce('+', Hg), rowSums(Sg))
  # baseline hazard
  ft <- cbind(1, sv$ft, 1, 1, sv$ft)
  lambda <- lambdaUpdate(sv$surv.times, ft, gamma, SigmaiSplit, bsplit, inds, w, v)
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
    D.new = D.new, beta.new = beta.new, var.e.new = var.e.new, gamma.new = gamma.new, 
    b.hat = b.hat, l0.new = l0.new, l0i.new = l0i.new, l0u.new = l0u.new
  ))
}

EM <- function(data, ph, survdata, gh = 3, tol = 0.01, verbose = F, post.process = T, quad = F){
  start <- proc.time()[3]
  n <- nrow(survdata)
  # Get data matrices
  m <- Y <- X <- Z <- list()
  for(i in 1:n){
    i.dat <- data[data$id == i, ]
    m[[i]] <- nrow(i.dat)
    Y[[i]] <- cbind(i.dat$Y.1, i.dat$Y.2, i.dat$Y.3)
    X[[i]] <- model.matrix(~time+cont+bin, i.dat)
    Z[[i]] <- setNames(
      list(model.matrix(~time, i.dat), model.matrix(~1, i.dat)),
      c('gc', 'b')
    )
  }
  # indices for extracting sub-matrices/columns.
  inds <- lapply(split(seq(5), rep(seq(3), c(2,1,2))), function(x) x-1)
  # initial conditions
  inits.long <- Longit.inits(data)
  b <- Ranefs(inits.long)
  beta <- inits.long$beta.init
  var.e <- inits.long$var.e.init
  V <- lapply(m, function(x) diag(x = var.e, nr = x, nc = x))
  D <- inits.long$D.init
  
  # pre-populate D step?
  # pre-populate D step?
  D[lower.tri(D, F)] <- cov(b)[lower.tri(cov(b), F)]
  D[upper.tri(D, F)] <- t(D)[upper.tri(D, F)]
  
  inits.surv <- TimeVarCox(data, b)
  b <- lapply(1:n, function(i) b[i, ])
  rm(inits.long) # large object
  
  # Survival objects
  sv <- surv.mod(ph, data, inits.surv$l0.init)
  Delta <- as.list(sv$Di)
  l0i <- as.list(sv$l0i)
  l0u <- sv$l0u
  Fi <- lapply(1:n, function(i) do.call(c, replicate(3, sv$Fi[i, ], simplify = F)))
  Fi <- lapply(Fi, function(x) x[-4]) # removing slope from binary part.
  Fu <- sv$Fu
  Fu <- lapply(Fu, function(x) cbind(x, 1, x))
  Fu.list <- lapply(Fu, function(y) lapply(inds, function(x) y[, x + 1, drop = F]))

  gamma <- inits.surv$inits
  
  # Quadrature //
  aa <- statmod::gauss.quad.prob(gh, 'normal')
  w <- aa$w; v <- aa$n
  
  # Collect parameters
  vD <- vech(D);
  names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
  params <- c(vD, beta, var.e, gamma)
  
  # Collect data objects
  data.mat <- list(Y = Y, X = X, Z = Z, m = m,                            # Longit.
                   Fi = Fi, Fu = Fu, Fu.list = Fu.list, Deltai = Delta)   # Survival
  
  EMstart <- proc.time()[3]
  diff <- 100; iter <- 0
  while(diff > tol){
    update <- EMupdate(b, Y, X, Z, V, D, beta, var.e, gamma, Delta, l0i, l0u, Fi, Fu, Fu.list, survdata, sv, w, v, inds, quad)
    params.new <- c(vech(update$D.new), update$beta.new, update$var.e.new, update$gamma.new)
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
    var.e <- update$var.e.new; V <- lapply(m, function(x) diag(x = var.e, nr = x, nc = x))
    gamma <- update$gamma.new
    l0 <- update$l0.new; l0u <- update$l0u.new; l0i <- update$l0i.new
    params <- params.new
    iter <- iter + 1
  }
  EMend <- proc.time()[3]
  gamma_r <- rep(gamma, c(2, 1, 2))
  coeffs <- list(D = D, beta = beta, var.e = var.e, gamma = gamma)
  out <- list(coeffs = coeffs,
              RE = do.call(rbind, b),
              EMtime = EMend-EMstart,
              iter = iter,
              totaltime = proc.time()[3] - start)
  out$hazard <- cbind(ft = sv$ft, haz = l0)
  out$quad <- quad
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time <- proc.time()[3]
    #' \hat{b} at MLEs
    b <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, haz){
      ucminf::ucminf(b, joint_density, joint_density_ddb,
                     X = X, Z = Z, beta = beta, V = V, D = D,
                     Y_1 = Y[, 1], Y_2 = Y[, 2], Y_3 = Y[, 3], 
                     Delta = Delta, Fi = Fi, l0i = l0i, Fu = Fu,
                     haz = haz, gamma = gamma_r)$par
    }, b = b, X = X, Y = Y, Z = Z, V = V, Delta = Delta, Fi = Fi, l0i = l0i, Fu = Fu, haz = l0u, SIMPLIFY = F)
    bmat <- lapply(b, function(x){
      matrix(c(x[1:3], 0, x[4:5]), nc = 2, by = T)
    })
    bsplit <- lapply(b, function(y) lapply(inds, function(x) y[x + 1]))
    # And its covariance matrix.
    Sigmai <- mapply(function(b, X, Z, V, Y, Delta, Fi, l0i, Fu, l0u){
      solve(joint_density_sdb(b = b, X = X, Z = Z, beta = beta, V = V, D = D,
                              Y_1 = Y[,1], Y_2 = Y[,2], Y_3 = Y[,3],
                              Delta = Delta, Fi = Fi, l0i = l0i, Fu = Fu, 
                              haz = l0u, gamma = gamma_r, eps = 1e-3))
    }, b = b, X = X, Z = Z, V = V, Y = Y, Delta = Delta, Fi = Fi,
    l0i= l0i, Fu = Fu, l0u = l0u, SIMPLIFY = F)
    # Split out into constituent block-diag pieces.
    SigmaiSplit <- lapply(Sigmai, function(y) lapply(inds, function(x) as.matrix(y[x + 1, x + 1])))
    
    # Calculate Standard Errors
    I <- structure(vcov(coeffs, data.mat, V, b, bsplit, bmat, Sigmai, SigmaiSplit, l0u, gh, n, quad),
                   dimnames = list(names(params), names(params)))
    out$SE <- sqrt(diag(solve(I)))
    out$vcov <- I
    out$RE <- do.call(rbind, b)
    out$postprocess.time <- round(proc.time()[3]-start.time, 2)
  }
  out
}