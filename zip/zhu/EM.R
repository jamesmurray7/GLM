#' #############
#' EM.R // Approximate EM algorithm for univariate ZIP longitudinal process and 
#'         a survival model
#'     --> Mimicks simulation 2 in Zhu et al. (2018).
#' Assumes working directory is set to ~/.../GLMM/zip/zhu
#' #############
library(survival)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)

sourceCpp('zip.cpp')
source('./_Functions.R')
source('../survFnsInt.R')
source('./inits.R')
vech <- function(x) x[lower.tri(x,diag=T)]

EMupdate <- function(b, Y, X, Z, Xz, Zz, 
                     beta, D, alpha, indzi,
                     gamma, K, eta, l0u, KK, Fu, l0i, Delta, w, v, survdata, sv){
  n <- length(Y)
  # E-step
  b.hat <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta){
    ucminf::ucminf(b, joint_density, NULL, # joint_density_ddb, 
                   Y, X, Z, Xz, Zz, beta, alpha, D, indzi, gamma, K, eta, l0u, KK, Fu, l0i, Delta)$par
  }, b = b, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
  KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), SIMPLIFY = F)

  Sigmai <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta){
    solve(joint_density_sdb(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, gamma, K, eta, l0u, KK, Fu, l0i, Delta, eps = 1e-4))
  }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
  KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), SIMPLIFY = F)
  
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, b = b.hat, S = Sigmai, SIMPLIFY = F)
  
  ba <- mapply(function(b, Y, X, Z, Xz, Zz){
    out <- list()
    out[[1]] <- pracma::grad(b_logdensity2, c(beta, alpha), b = b, Y=Y,X=X,Z=Z,Xz=Xz,Zz=Zz,
                             beta_length = length(beta), alpha_length = length(alpha), D=D, indzi=2)
    out[[2]] <- pracma::hessian(b_logdensity2, c(beta, alpha), b = b, Y=Y,X=X,Z=Z,Xz=Xz,Zz=Zz,
                                beta_length = length(beta), alpha_length = length(alpha), D=D, indzi=2)
    out
  }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, SIMPLIFY = F)
  
  tau <- mapply(function(Fu, S){
    diag(tcrossprod(Fu %*% S, Fu))
  }, Fu = Fu, S = Sigmai, SIMPLIFY = F)
  
  Sge <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta, tau){
    Sgammaeta(c(gamma,eta), b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, K, 
              l0u, KK, Fu, l0i, Delta, tau, w, v)
  },  b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
  KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), tau = tau, SIMPLIFY = F)
  
  Hge <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta, tau){
    Hgammaeta(c(gamma,eta), b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, K, 
              l0u, KK, Fu, l0i, Delta, tau, w, v, 1e-4)
  },  b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
  KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), tau = tau, SIMPLIFY = F)
  
  # M-step
  # D covariance matrix
  D.new <- (Reduce('+', Sigmai) + Reduce('+', lapply(b.hat, tcrossprod)))/250
  # ZIP and Poisson processes coefficients
  Sba <- colSums(do.call(rbind, lapply(ba, '[[', 1)))
  Hba <- Reduce('+', lapply(ba, '[[', 2))
  
  beta.alpha.new <- c(beta,alpha) - solve(Hba,Sba)
  beta.new <- beta.alpha.new[1:3]
  alpha.new <- beta.alpha.new[4:6]

  # Survival parameters (gamma, eta)
  Sge <- rowSums(do.call(cbind, Sge))
  Imat <- Reduce('+', Hge)
  gammaeta.new <- c(gamma, eta) - solve(Imat, Sge)
  gamma.new <- gammaeta.new[1:2]
  eta.new <- gammaeta.new[3]
  
  # Update for baseline hazard
  lambda <- lambdaUpdate(sv$surv.times, sv$ft, gamma, eta, K, Sigmai, b.hat, id=n, w, v)
  l0.new <- sv$nev/rowSums(lambda)
  l0u.new <- lapply(l0u, function(x){
    ll <- length(x); l0.new[1:ll]
  })
  l0i.new <- c()
  l0i.new[which(Delta == 0)] <- 0 
  l0i.new[which(Delta == 1)] <- l0.new[match(survdata[which(Delta==1), 'survtime'], sv$ft)]
  
  return(list(
    D.new = D.new, beta.new = beta.new, alpha.new = alpha.new, b.hat = b.hat,
    gamma.new = gamma.new, eta.new = eta.new,
    l0.new = l0.new, l0u.new = l0u.new, l0i.new = l0i.new
  ))
}

EM <- function(data, ph, survdata, indzi = 2, gh = 9, tol = 1e-2, verbose = F){
  start <- proc.time()[3]
  inits.long <- Longit.inits(data)
  beta <- inits.long$beta.init
  alpha <- inits.long$alpha.init
  D <- inits.long$D.init
  b <- Ranefs(inits.long); n <- nrow(b)
  inits.surv <- TimeVarCox(data, b)
  b <- lapply(1:n, function(i) b[i, ])
  # Get data matrices...
  X <- Y <- Z <- Xz <- Zz <- K <- list()
  for(i in 1:n){
    i.dat <- data[data$id == i, ]
    X[[i]] <- model.matrix(~time + bin, i.dat)
    Xz[[i]] <- model.matrix(~time + bin, i.dat)
    Z[[i]] <- Zz[[i]] <- model.matrix(~1, i.dat)
    Y[[i]] <- i.dat$Y
    K[[i]] <- as.matrix(unname(unique(i.dat[,'bin'])))
    row.names(K[[i]]) <- NULL
  }      
  # Survival objects
  sv <- surv.mod(ph, data, inits.surv$l0.init)
  Fu <- sv$Fu
  Delta <- sv$Di
  l0i <- sv$l0i
  l0u <- sv$l0u
  KK <- sapply(1:n, function(i){
    x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
    if('numeric' %in% class(x)) x <- t(as.matrix(x))
    x
  })
  gamma <- inits.surv$inits[2:3]
  eta <- inits.surv$inits[1]
  
  aa <- statmod::gauss.quad.prob(gh, 'normal')
  w <- aa$w; v <- aa$n
  
  EMstart <- proc.time()[3]
  params <- c(vech(D), beta, alpha, gamma, eta)
  diff <- 100; iter <- 0
  while(diff > tol){
    update <- EMupdate(b, Y, X, Z, Xz, Zz, beta, D, alpha, 2, gamma, K, eta, l0u, KK, Fu, l0i, Delta, w, v, survdata, sv)
    params.new <- c(vech(update$D.new), update$beta.new, update$alpha.new, update$gamma.new, update$eta.new)
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
    eta <- update$eta.new
    l0 <- update$l0.new; l0u <- update$l0u.new; l0i <- update$l0i.new
    params <- c(vech(D), beta, alpha, gamma, eta)
    iter <- iter + 1
  }
  EMend <- proc.time()[3]
  coeffs <- list(D = D, beta = beta, alpha = alpha, gamma = gamma, eta = eta)
  out <- list(coeffs = coeffs,
              RE = do.call(rbind, b),
              inits = inits.long,
              EMtime = EMend-EMstart,
              iter = iter,
              totaltime = proc.time()[3] - start)
  out
}
