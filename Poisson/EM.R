#' ####    Poisson/
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
sourceCpp('po.cpp')
vech <- function(x) x[lower.tri(x, diag = T)]

EMupdate <- function(b, Y, X, Z, 
                     D, beta, gamma, eta,
                     Delta, l0i, l0u, Fi, Fu, K, KK, survdata, sv, w, v, beta.quad = F){
  n <- length(b)
  
  #' ### ------
  #' E step
  #' ### ------
  
  b.hat <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, haz){
    ucminf::ucminf(b, joint_density, joint_density_ddb,
                   X = X, Y = Y, Z = Z, beta = beta, D = D,
                   Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu,
                   haz = haz, gamma = gamma, eta = eta)$par
  }, b = b, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
  KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)
  
  Sigmai <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, haz){
    solve(joint_density_sdb(b, X, Y, Z, beta, D, Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta, eps = 1e-4))
  }, b = b.hat, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
  KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)
  
  # Update to D
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, S = Sigmai, b = b.hat, SIMPLIFY = F)

  # Score and Hessian for \beta
  if(beta.quad){
    Sb <- mapply(function(b, X, Y, Z, S){
      pracma::grad(pois_ll_quad, beta, X = X, Y = Y, Z = Z, b = b, S = S, w = w$w2, v = v$v2)
    }, b = b.hat, X = X, Y = Y, Z = Z, S = Sigmai, SIMPLIFY = F)
    
    Hb <- mapply(function(b, X, Y, Z, S){
      pracma::hessian(pois_ll_quad, beta, X = X, Y = Y, Z = Z, b = b, S = S, w = w$w2, v = v$v2)
    }, b = b.hat, X = X, Y = Y, Z = Z, S = Sigmai, SIMPLIFY = F)
  }  else{
    Sb <- mapply(function(X, Y, Z, b){
      Sbeta(beta, X, Y, Z, b)
    }, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)
    
    Hb <- mapply(function(X, Y, Z, b){
      Hbeta(beta, X, Y, Z, b, 1e-4)
    }, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)
  }
  
  # Score and Hessian for (gamma, eta)
  Sge <- mapply(function(b, Delta, Fi, K, KK, Fu, l0u, S){
    Sgammaeta(c(gamma, eta), b, S, K, KK, Fu, Fi, l0u, Delta, w$w, v$v, 1e-4)
  }, b = b.hat, Delta = Delta, Fi = Fi, K = K, KK = KK, Fu = Fu, l0u = l0u, S = Sigmai)
  
  Hge <- mapply(function(b, Delta, Fi, K, KK, Fu, l0u, S){
    Hgammaeta(c(gamma, eta), b, S, K, KK, Fu, Fi, l0u, Delta, w$w, v$v, 1e-4)
  }, b = b.hat, Delta = Delta, Fi = Fi, K = K, KK = KK, Fu = Fu, l0u = l0u, S = Sigmai, SIMPLIFY = F)
  
  #' ### ------
  #' M step
  #' ### ------
  
  D.new <- Reduce('+', Drhs)/n
  beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
  gammaeta.new <- c(gamma, eta) - solve(Reduce('+', Hge), rowSums(Sge))
  lambda <- lambdaUpdate(sv$surv.times, sv$ft, gamma, eta, K, Sigmai, b.hat, w$w, v$v)
  # Baseline hazard
  l0.new <- sv$nev/rowSums(lambda)
  l0u.new <- lapply(l0u, function(x){
    ll <- length(x); l0.new[1:ll]
  })
  l0i.new <- c()
  l0i.new[which(unlist(Delta) == 0)] <- 0 
  l0i.new[which(unlist(Delta) == 1)] <- l0.new[match(survdata[which(unlist(Delta)==1), 'survtime'], sv$ft)]
  l0i.new <- as.list(l0i.new)
  
  return(list(
    D.new = D.new, beta.new = beta.new, 
    gamma.new = gammaeta.new[1], eta.new = gammaeta.new[2:3], 
    b.hat = b.hat, l0.new = l0.new, l0i.new = l0i.new, l0u.new = l0u.new
  ))
}

EM <- function(data, ph, survdata, gh = 9, tol = 0.01, post.process = T, verbose = F, beta.quad = F){
  start <- proc.time()[3]
  inits.long <- Longit.inits(data)
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- Ranefs(inits.long); n <- nrow(b)
  inits.surv <- TimeVarCox(data, b)
  gamma <- inits.surv$inits[3]
  eta <- inits.surv$inits[1:2]
  b <- lapply(1:n, function(i) b[i, ])
  
  # Get data matrices
  K <- Y <- X <- Z <- list()
  for(i in 1:n){
    i.dat <- data[data$id == i, ]
    Y[[i]] <- i.dat$Y
    X[[i]] <- model.matrix(~time + cont + bin, i.dat)
    Z[[i]] <- model.matrix(~time, i.dat)
    K[[i]] <- unname(cbind(unique(i.dat$cont), unique(i.dat$bin)))
  }   
  m <- lapply(Y, length)
  
  # Survival objects
  sv <- surv.mod(ph, data, inits.surv$l0.init)
  Fi <- sv$Fi; Fu <- sv$Fu;
  Delta <- as.list(sv$Di)
  l0i <- as.list(sv$l0i)
  l0u <- sv$l0u
  Fi <- lapply(1:n, function(i) Fi[i, ])
  KK <- sapply(1:n, function(i){
    x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
    if('numeric' %in% class(x)) x <- t(as.matrix(x))
    x
  })
  
  # Quadrature //
  aa <- statmod::gauss.quad.prob(gh, 'normal') # For survival part
  w <- aa$w; v <- aa$n
  bb <- statmod::gauss.quad.prob(gh, 'normal', sigma = sqrt(0.5))
  w2 <- bb$w; v2 <- bb$n
  w <- list(w=w,w2=w2); v <- list(v=v,v2=v2)
  
  
  # Collect parameters
  vD <- vech(D);
  names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
  params <- c(vD, beta, gamma, eta)
  
  # Collect data objects
  data.mat <- list(Y = Y, X = X, Z = Z, m = m, K = K, KK = KK, Fi = Fi, Fu = Fu, Delta = Delta)
  
  # Start EM
  diff <- 100; iter <- 0
  EMstart <- proc.time()[3]
  while(diff > tol){
    update <- EMupdate(b, Y, X, Z, D, beta, gamma, eta, Delta, l0i, l0u, Fi, Fu, K, KK, survdata, sv, w, v, beta.quad)
    params.new <- c(vech(update$D.new), update$beta.new, update$gamma.new, update$eta.new)
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
    gamma <- update$gamma.new
    eta <- update$eta.new
    l0 <- update$l0.new; l0u <- update$l0u.new; l0i <- update$l0i.new
    params <- params.new
    iter <- iter + 1
  }
  EMend <- proc.time()[3]
  coeffs <- list(D = D, beta = beta, gamma = gamma, eta = eta)
  out <- list(coeffs = coeffs,
              RE = do.call(rbind, b),
              inits = inits.long,
              EMtime = EMend-EMstart,
              iter = iter,
              totaltime = proc.time()[3] - start)
  out$hazard <- cbind(ft = sv$ft, haz = l0)
  out$beta.quad <- beta.quad
  
  # Post Processing
  if(post.process){
    message('\nCalculating SEs...')
    start.time <- proc.time()[3]
    #' b and Sigmai at MLEs
    b <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, haz){
      ucminf::ucminf(b, joint_density, joint_density_ddb,
                     X = X, Y = Y, Z = Z, beta = beta, D = D,
                     Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu,
                     haz = haz, gamma = gamma, eta = eta)$par
    }, b = b, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
    KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)
    Sigmai <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, haz){
      solve(joint_density_sdb(b, X, Y, Z, beta, D, Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta, eps = 1e-4))
    }, b = b, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
    KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)
   
    I <- structure(vcov(coeffs, data.mat, b, Sigmai, l0u, gh, n),
                   dimnames = list(names(params), names(params)))
    out$vcov <- I
    out$SE <- sqrt(diag(solve(I)))
    out$postprocess.time <- round(proc.time()[3]-start.time, 2)
    out$RE <- do.call(rbind, b)
  }
  out
}