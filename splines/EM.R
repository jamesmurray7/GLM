#' ####    quad/
#' Approximate EM for joint model of (three) quadratic Gaussian responses
#' ####
library(dplyr)
library(survival)
library(Rcpp)
library(RcppArmadillo)
source('simData.R')
source('inits.R')
source('survFns.R')
source('MVLME.R')
sourceCpp('splines.cpp')
sourceCpp('mvlme.cpp')
vech <- function(x) x[lower.tri(x, diag = T)]

EMupdate <- function(b, Y, X, XtX, Z, V, m,        # for beta update and minimisation (X and Z are BLOCK matrices)
                     Ymat, X_single, Z_single,     # for var.e update                 (X, Z design matrices on single outcomes)
                     D, beta, var.e, gamma, eta,   # parameters
                     Delta, l0i, l0u, Fi, Fu, K, KK, 
                     survdata, sv, w, v, nK, beta.inds, degree, basis){
  n <- length(b)
  gh <- length(w)
  #' ### ------
  #' E step
  #' ### ------
  b.hat <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, l0u){
    ucminf::ucminf(b, joint_density, joint_density_db,
                   X = X, Y = Y, Z = Z, beta = beta, V = V, D = D,
                   K = K, KK = KK, Fi = Fi, Fu = Fu, l0u = l0u, l0i = l0i,
                   Delta = Delta, gamma = rep(gamma, each = (degree + 1)), eta = eta,
                   control = list(xtol = 1e-4, grtol = 1e-5), hessian = 0)$par
  }, b = b, X = X, Y = Y, Z = Z, V = V, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
  KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)
  # Split out into response-oriented objects.
  bmat <- lapply(b.hat, matrix, nc = (degree + 1), byr = T)
  bsplit <- lapply(b.hat, function(y) lapply(split(seq(nK * (degree + 1)), rep(seq(nK), each = (degree + 1))), function(x) y[x]))
  
  # Posterior variance
  Sigmai <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, l0u){
    solve(joint_density_sdb(b=b, X = X, Y = Y, Z = Z, beta = beta, V = V, D = D,
                            K = K, KK = KK, Fi = Fi, Fu = Fu, l0u = l0u, l0i = l0i,
                            Delta = Delta, gamma = rep(gamma, each = (degree + 1)), eta = eta))
  }, b = b.hat, X = X, Y = Y, Z = Z, V = V, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
  KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)
  
  # Split into nK constituent sub-matrices
  S <- lapply(Sigmai, function(y) lapply(split(seq(nK * (degree + 1)), rep(seq(nK), each = (degree + 1))), function(x) y[x, x]))
  
  #' Step to update D -----
  Drhs <- mapply(function(S, b) S + tcrossprod(b), S  = Sigmai, b = b.hat, SIMPLIFY = F)
  
  #' Step to update beta ------
  tau.long <- mapply(function(S, Z) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigmai, Z = Z, SIMPLIFY = F)
  
  beta.rhs <- mapply(function(b, X, Y, Z, tau){
    mu <- Z %*% b
    rhs <- numeric(length(tau))
    for(l in 1:gh) rhs <- rhs + w[l] * tau * v[l]
    crossprod(X, Y - mu - rhs)
  }, b = b.hat, X = X, Y = Y, Z = Z, tau = tau.long, SIMPLIFY = F)
  
  #' Step to update sigma^2 ------
  tau.longK <- mapply(function(S, Z){
    out <- vector('list', nK)
    for(k in 1:nK) out[[k]] <- sqrt(diag(tcrossprod(Z %*% S[[k]], Z)))   # exploiting set-up of Z being the same
    out
  }, S = S, Z = Z_single, SIMPLIFY = F)
  
  mu.longK <- mapply(function(X, Z, b){
    out <- list()
    for(k in 1:nK) out[[k]] <- X %*% beta[beta.inds[[k]]] + Z %*% b[[k]]
    out
  }, X = X_single, Z = Z_single, b = bsplit , SIMPLIFY = F)
  
  Esigma <- mapply(function(Y, mu, tau){
    temp <- matrix(NA, nr = gh, nc = nK)
    for(k in 1:nK){
      for(l in 1:gh){
        temp[l, k] <- w[l] * crossprod(Y[, k] - mu[[k]] - tau[[k]] * v[l])
      }
    }
    colSums(temp)
  }, Y = Ymat, mu = mu.longK, tau = tau.longK)
  
  #' Step to update (gamma, eta) ------
  Sge <- mapply(function(bmat, S, K, KK, Fu, Fi, l0u, Delta){
    Sgammaeta(c(gamma, eta), bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi[1:(degree + 1)], haz = l0u,
              Delta = Delta, w = w, v = v, eps = 1e-4)
  }, bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)
  
  Hge <- mapply(function(bmat, S, K, KK, Fu, Fi, l0u, Delta){
    Hgammaeta(c(gamma, eta), bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi[1:(degree + 1)], haz = l0u,
              Delta = Delta, w = w, v = v, eps = 1e-4)
  }, bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  #' ###
  #' M-step
  #' ###
  
  D.new <- Reduce('+', Drhs) / n
  beta.new <- c(solve(Reduce('+', XtX)) %*% Reduce('+', beta.rhs))
  var.e.new <- rowSums(Esigma) / colSums(do.call(rbind, m))
  gammaeta.new <- c(gamma, eta) - solve(Reduce('+', Hge), rowSums(Sge))
  
  # baseline hazard
  lambda <- lambdaUpdate(sv$surv.times, sv$ft, basis, gamma, eta, K, S, bsplit, w, v)
  # Baseline hazard objects
  l0 <- sv$nev/rowSums(lambda)
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
    gamma.new = gammaeta.new[1:nK], eta.new = gammaeta.new[(nK + 1):length(gammaeta.new)], 
    b.hat = b.hat, l0 = l0.new, l0i.new = l0i.new, l0u.new = l0u.new
  ))
}

EM <- function(data, ph, survdata, gh = 3, tol = 0.01, nK = 3, degree = 3, MVLME = T, tol.mvlme = 1e-2, verbose = F){
  start <- proc.time()[3]
  n <- nrow(survdata)
  basis <- getsurvbasis(data, degree)
  # Get data matrices
  Z <- getZfromsurvbasis(basis, data)
  m <- Y <- Ymat <- X <- K <- list()
  for(i in 1:n){
    i.dat <- data[data$id == i, ]
    m[[i]] <- rep(nrow(i.dat), 3)
    Y[[i]] <- c(i.dat$Y.1, i.dat$Y.2, i.dat$Y.3)
    Ymat[[i]] <- cbind(i.dat$Y.1, i.dat$Y.2, i.dat$Y.3)
    colnames(Ymat[[i]]) <- paste0('Y.', 1:nK)
    X[[i]] <- cbind(Z[[i]], i.dat$cont, i.dat$bin)
    colnames(X[[i]])[(degree + 2):ncol(X[[i]])] <- c('cont', 'bin')
    K[[i]] <- unname(cbind(unique(i.dat$cont), unique(i.dat$bin)))
  }
  # Block matrices
  Xblock <- lapply(X, function(x) as.matrix(Matrix::bdiag(replicate(nK, x, simplify = F))))
  Zblock <- lapply(Z, function(x) as.matrix(Matrix::bdiag(replicate(nK, x, simplify = F))))
  XtX <- lapply(Xblock, crossprod)
  
  # initial conditions
  message('Getting initial conditions')
  df.basis <- as.data.frame(cbind(id = do.call(c, lapply(1:n, function(i) rep(i, m[[i]][1]))), do.call(rbind, Ymat), do.call(rbind, X)))
  inits.long <- Longit.inits(nK, df.basis, degree)
  # TO DO: Add methods for doing this 'auto' way and add option to EM call.
  inits.surv <- TimeVarCox2(data, Ranefs(inits.long), basis, nK = nK)
  
  # Get initial conditions from additional MVLME step, or just from univariate LME fits?
  if(MVLME){
    mvlme.inits <- mvlme(data, Y, Xblock, Zblock, Ymat, X, Z, m, inits.long, nK, degree, tol.mvlme = tol.mvlme)
    b <- do.call(rbind, lapply(mvlme.inits$b, c))
    beta <- c(mvlme.inits$beta); names(beta) <- names(inits.long$beta.init)
    var.e <- mvlme.inits$var.e; names(var.e) <- names(inits.long$var.e.init)
    V <- lapply(m, function(iii) {
      diag(x = rep(var.e, iii), ncol = sum(iii))
    })
    D <- mvlme.inits$D
    b <- lapply(1:n, function(i) b[i, ])
  }else{
    b <- Ranefs(inits.long)
    beta <- inits.long$beta.init
    var.e <- inits.long$var.e.init
    V <- lapply(m, function(iii) {
      diag(x = rep(var.e, iii), ncol = sum(iii))
    })
    D <- inits.long$D.init
    b <- lapply(1:n, function(i) b[i, ])
  }
  
  sv <- surv.mod(ph, data, inits.surv$l0.init, degree = degree)
  Delta <- as.list(sv$Di)
  l0i <- as.list(sv$l0i)
  l0u <- sv$l0u
  Fi <- lapply(1:n, function(i) do.call(c, replicate(nK, sv$Fi[i, ], simplify = F)))
  Fu <- sv$Fu
  KK <- sapply(1:n, function(i){
    x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
    if('numeric'%in%class(x)) x <- t(as.matrix(x))
    x
  })
  gamma <- inits.surv$inits[3:length(inits.surv$inits)]
  eta <- inits.surv$inits[1:2]
  # basis on just failure times, for update to \lambda
  basis.ft <- cbind(1, sv$basisdf[match(sv$ft, sv$basisdf[, 1]), -1])
  
  # Quadrature //
  aa <- statmod::gauss.quad.prob(gh, 'normal')
  w <- aa$w; v <- aa$n
  
  # Initialise parameter vector
  vD <- vech(D);
  names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
  params <- c(vD, beta, var.e, gamma, eta)
  
  # indices for beta and b
  beta.inds <- lapply(1:nK, function(k){
    which(grepl(paste0('^beta', k, '_'), names(beta)))
  }) 
  
  message('Done!\nStarting EM...\n')
  EMstart <- proc.time()[3]
  diff <- 100; iter <- 0
  while(diff > tol){
    update <- EMupdate(b, Y, Xblock, XtX, Zblock, V, m, 
                       Ymat, X, Z, 
                       D, beta, var.e, gamma, eta, 
                       Delta, l0i, l0u, Fi, Fu, K, KK, survdata, sv, w, v, nK, beta.inds, degree, basis.ft)
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
    var.e <- update$var.e.new; V <- lapply(m, function(iii) diag(x = rep(var.e, iii), ncol = sum(iii)))
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
              EMtime = EMend-EMstart,
              iter = iter,
              totaltime = proc.time()[3] - start)
  out
}
