#' ####    negbinom/
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
sourceCpp('nb.cpp')
vech <- function(x) x[lower.tri(x, diag = T)]

EMupdate <- function(b, Y, X, Z, W,
                     D, beta, alpha, gamma, zeta,
                     Delta, l0i, l0u, Fi, Fu, S, SS, survdata, sv, w, v){
  n <- length(b)
  
  #' ### ------
  #' E step
  #' ### ------
  
  b.hat <- mapply(function(b, X, Y, Z, W, Delta, S, Fi, l0i, SS, Fu, haz){
    optim(b, joint_density, joint_density_ddb,
                   X = X, Y = Y, Z = Z, W = W, beta = beta, alpha = alpha, D = D,
                   Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu,
                   haz = haz, gamma = gamma, zeta = zeta, method = 'BFGS')$par
  }, b = b, X = X, Y = Y, Z = Z, W = W,Delta = Delta, S = S, Fi = Fi, l0i = l0i,
  SS = SS, Fu = Fu, haz = l0u, SIMPLIFY = F)
  
  Sigmai <- mapply(function(b, X, Y, Z, W, Delta, S, Fi, l0i, SS, Fu, haz){
    solve(joint_density_sdb(b, X, Y, Z, W, beta, alpha, D, Delta, S, Fi, l0i, SS, Fu, haz, gamma, zeta, eps = 1e-4))
  }, b = b.hat, X = X, Y = Y, Z = Z, W = W,Delta = Delta, S = S, Fi = Fi, l0i = l0i,
  SS = SS, Fu = Fu, haz = l0u, SIMPLIFY = F)
  
  # Update to D
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, S = Sigmai, b = b.hat, SIMPLIFY = F)
  
  # Update for \beta
  theta <- mapply(function(W) exp(W %*% alpha), W = W, SIMPLIFY = F)
  Sb <- Sbeta(beta, X, Y, Z, b.hat, theta)
  Hb <- Hbeta(beta, X, Y, Z, b.hat, theta, 1e-4)
  
  # Score and Hessian for \theta
  Sa <- mapply(function(X, Y, Z, W, b, Sigma){
    Salpha(alpha, beta, X, Y, Z, W, b, Sigma, w, v, 1e-4)
  }, X = X, Y = Y, Z = Z, W = W,b = b.hat, Sigma = Sigmai, SIMPLIFY = F)
    
  Ha <- mapply(function(X, Y, Z, W, b, Sigma){
    Halpha(alpha, beta, X, Y, Z, W, b, Sigma, w, v, 1e-3, 1e-5)
  }, X = X, Y = Y, Z = Z, W = W, b = b.hat, Sigma = Sigmai, SIMPLIFY = F)
  
  # Score and Hessian for (gamma, eta)
  Sge <- mapply(function(b, Delta, Fi, S, SS, Fu, l0u, Sigma){
    Sgammaeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 1e-4)
  }, b = b.hat, Delta = Delta, Fi = Fi, S = S, SS = SS, Fu = Fu, l0u = l0u, Sigma = Sigmai)
  
  Hge <- mapply(function(b, Delta, Fi, S, SS, Fu, l0u, Sigma){
    Hgammaeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 1e-4)
  }, b = b.hat, Delta = Delta, Fi = Fi, S = S, SS = SS, Fu = Fu, l0u = l0u, Sigma = Sigmai, SIMPLIFY = F)
  
  #' ### ------
  #' M step
  #' ### ------
  
  D.new <- Reduce('+', Drhs)/n
  beta.new <- beta - solve(Hb, Sb)
  alpha.new <- alpha - solve(Reduce('+', Ha), Reduce('+', Sa))
  gammaeta.new <- c(gamma, zeta) - solve(Reduce('+', Hge), rowSums(Sge))
  lambda <- lambdaUpdate(sv$surv.times, sv$ft, gamma, zeta, S, Sigmai, b.hat, w, v)
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
    D.new = D.new, beta.new = beta.new, alpha.new = alpha.new,
    gamma.new = gammaeta.new[1], zeta.new = gammaeta.new[2:3], 
    b.hat = b.hat, l0.new = l0.new, l0i.new = l0i.new, l0u.new = l0u.new
  ))
}

EM <- function(data, ph, survdata, gh = 3, dispformula = ~time,
               tol = 0.01, post.process = T, verbose = F, beta.quad = F){
  start <- proc.time()[3]
  inits.long <- Longit.inits(data, dispformula)
  beta <- inits.long$beta.init
  alpha <- inits.long$alpha.init # dispersion parameter (pre-exp)
  D <- inits.long$D.init
  b <- Ranefs(inits.long); n <- nrow(b)
  inits.surv <- TimeVarCox(data, b)
  gamma <- inits.surv$inits[3]
  zeta <- inits.surv$inits[1:2]
  b <- lapply(1:n, function(i) b[i, ])
  
  # Get data matrices
  S <- Y <- X <- Z <- W <- list()
  for(i in 1:n){
    i.dat <- data[data$id == i, ]
    Y[[i]] <- i.dat$Y
    X[[i]] <- model.matrix(~time + cont + bin, i.dat)
    Z[[i]] <- model.matrix(~time, i.dat)
    W[[i]] <- model.matrix(dispformula, i.dat)
    S[[i]] <- unname(cbind(unique(i.dat$cont), unique(i.dat$bin)))
  }   
  m <- lapply(Y, length)
  # Survival objects
  sv <- surv.mod(ph, data, inits.surv$l0.init)
  Fi <- sv$Fi; Fu <- sv$Fu;
  Delta <- as.list(sv$Di)
  l0i <- as.list(sv$l0i)
  l0u <- sv$l0u
  Fi <- lapply(1:n, function(i) Fi[i, ])
  SS <- sapply(1:n, function(i){
    x <- apply(S[[i]], 2, rep, nrow(Fu[[i]]))
    if('numeric' %in% class(x)) x <- t(as.matrix(x))
    x
  })
  
  # Quadrature //
  aa <- statmod::gauss.quad.prob(gh, 'normal')
  w <- aa$w; v <- aa$n
  # Collect parameters
  vD <- vech(D);
  names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
  params <- c(vD, beta, alpha, gamma, zeta)
  # Collect data objects
  data.mat <- list(Y = Y, X = X, Z = Z, W = W, m = m, S = S, SS = SS, Fi = Fi, Fu = Fu, Delta = Delta)
  # EM
  diff <- 100; iter <- 0
  EMstart <- proc.time()[3]
  while(diff > tol){
    update <- EMupdate(b, Y, X, Z, W, D, beta, alpha, gamma, zeta, Delta, l0i, l0u, Fi, Fu, S, SS, survdata, sv, w, v)
    params.new <- c(vech(update$D.new), update$beta.new, update$alpha.new, update$gamma.new, update$zeta.new)
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
    alpha <- update$alpha.new
    gamma <- update$gamma.new
    zeta <- update$zeta.new
    l0 <- update$l0.new; l0u <- update$l0u.new; l0i <- update$l0i.new
    params <- params.new
    iter <- iter + 1
  }
  EMend <- proc.time()[3]
  coeffs <- list(D = D, beta = beta, alpha = alpha, gamma = gamma, zeta = zeta)
  out <- list(coeffs = coeffs,
              RE = do.call(rbind, b),
              inits = inits.long,
              EMtime = EMend-EMstart,
              iter = iter,
              totaltime = proc.time()[3] - start)
  out$hazard <- cbind(ft = sv$ft, haz = l0)
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time <- proc.time()[3]
    #' b and Sigmai at MLEs
    b <- mapply(function(b, X, Y, Z, W, Delta, S, Fi, l0i, SS, Fu, haz){
      optim(b, joint_density, joint_density_ddb,
            X = X, Y = Y, Z = Z, W = W, beta = beta, alpha = alpha, D = D,
            Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu,
            haz = haz, gamma = gamma, zeta = zeta, method = 'BFGS')$par
    }, b = b, X = X, Y = Y, Z = Z, W = W,Delta = Delta, S = S, Fi = Fi, l0i = l0i,
    SS = SS, Fu = Fu, haz = l0u, SIMPLIFY = F)
    
    Sigmai <- mapply(function(b, X, Y, Z, W, Delta, S, Fi, l0i, SS, Fu, haz){
      solve(joint_density_sdb(b, X, Y, Z, W, beta, alpha, D, Delta, S, Fi, l0i, SS, Fu, haz, gamma, zeta, eps = 1e-4))
    }, b = b, X = X, Y = Y, Z = Z, W = W,Delta = Delta, S = S, Fi = Fi, l0i = l0i,
    SS = SS, Fu = Fu, haz = l0u, SIMPLIFY = F)
    
    I <- structure(vcov(coeffs, data.mat, b, Sigmai, l0u, gh, n),
                   dimnames = list(names(params), names(params)))
    out$SE <- sqrt(diag(solve(I)))
    out$vcov <- I
    out$RE <- do.call(rbind, b)
    out$postprocess.time <- round(proc.time()[3]-start.time, 2)
  }
  
  out
}
