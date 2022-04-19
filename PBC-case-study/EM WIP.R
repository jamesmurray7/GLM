library(survival)
library(glmmTMB)
library(Rcpp)
library(RcppArmadillo)
source('_Functions.R')
source('simData.R')
source('inits.R')
sourceCpp('funs.cpp')
vech <- function(x) x[lower.tri(x, T)]


EMupdate <- function(Omega, family, X, Y, Z, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, sv, w, v, n, m){
  
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta

  #' Posterior mode via ucminf ----
  b.hat <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    ucminf::ucminf(b, joint_density, joint_density_ddb,
                   Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
                   Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma = gamma, zeta = zeta)$par
  }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
     Fu = Fu, l0u = l0u, SIMPLIFY = F)
  
  #' And its variance ----
  Sigma <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    solve(
      joint_density_sdb(b = b, Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
                        Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma = gamma, 
                        zeta = zeta, eps = 1e-3))
  }, b = b.hat, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
  Fu = Fu, l0u = l0u, SIMPLIFY = F)
  
  #' #########
  #' E-step ##
  #' #########
  
  # D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  # \beta
  Sb <- mapply(function(X, Y, Z, b){
    Sbeta(beta, X, Y, Z, b, sigma, family)
  }, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)
  Hb <- mapply(function(X, Y, Z, b){
    Hbeta(beta, X, Y, Z, b, sigma, family, .Machine$double.eps^(1/4))
  }, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)
  
  # Dispersion ('\sigma')
  if(family == 'gaussian'){
    tau <- mapply(function(Sigma, Z) sqrt(diag(tcrossprod(Z %*% Sigma, Z))), Z = Z, Sigma = Sigma, SIMPLIFY = F)
    sigma.update <- mapply(function(X, Y, Z, b, tau){
      vare_update(X, Y, Z, b, beta, tau, w, v)
    }, X = X, Y = Y, Z = Z, b = b.hat, tau = tau, SIMPLIFY = T)
  }else if(family == 'negative.binomial'){
    sigma.update <- mapply(function(X, Y, Z, b, Sigma){
      out <- vector('list', 2)
      out[[1]] <- Stheta(sigma, beta, X, Y, Z, b, Sigma, w, v, .Machine$double.eps^(1/3))
      out[[2]] <- Htheta(sigma, beta, X, Y, Z, b, Sigma, w, v, .Machine$double.eps^(1/4))
      out
    }, X = X, Y = Y, Z = Z, b = b.hat, Sigma = Sigma, SIMPLIFY = F)
  }else{
    sigma.update <- NULL
  }
  
  # Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, .Machine$double.eps^(1/3))
  }, b = b.hat, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)
  
  Hgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Hgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, .Machine$double.eps^(1/4))
  }, b = b.hat, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  #' #########
  #' M-step ##
  #' #########
  
  # D
  D.new <- Reduce('+', D.update)/n
  
  # beta
  beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
  
  # Dispersion
  if(family == 'gaussian'){
    sigma.new <- sum(sigma.update)/sum(m)
  }else if(family == 'negative.binomial'){
    sigma.new <- sigma - sum(unlist(lapply(sigma.update, el, 1)))/sum(unlist(lapply(sigma.update, el, 2)))
  }else{
    sigma.new <- 0
  }
  
  # Survival parameters (gamma, zeta)
  gammazeta.new <- c(gamma, zeta) - solve(Reduce('+', Hgz), rowSums(Sgz))
  
  # The baseline hazard and related objects
  lambda.update <- lambdaUpdate(sv$surv.times, sv$ft.mat, gamma, zeta, S, Sigma, b.hat, w, v)
  l0.new <- sv$nev/rowSums(lambda.update)
  l0u.new <- lapply(l0u, function(ll){
    l0.new[1:length(ll)]
  })
  l0i.new <- l0.new[match(sv$Ti, sv$ft)] 
  l0i.new[is.na(l0i.new)] <- 0
  
  # Return
  list(
    D = D.new, beta = beta.new, sigma = sigma.new,       # <Y>
    gamma = gammazeta.new[1], zeta = gammazeta.new[-1],  # Survival
    l0 = l0.new, l0u = l0u.new, l0i = as.list(l0i.new),  # ---""---
    b = b.hat                                            #   REs.
  )
  
}



EM <- function(long.formula, surv.formula, data, family, post.process = T, control = list()){
  
  start.time <- proc.time()[3]
  
  #' Initial parsing ----
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- nrow(surv$survdata)
  
  #' Initial conditons ----
  if(!is.null(control$dispformula)) dispformula <- control$dispformula else dispformula <- NULL
  inits.long <- Longit.inits(long.formula, data, family, dispformula = dispformula)
  inits.surv <- TimeVarCox(data, inits.long$b, surv$ph)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  sigma <- inits.long$sigma.init # dispersion / resid. variance / 0 otherwise.
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  # Survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma', names(inits.surv$inits))]
  
  #' Longitudinal and survival data objects ----
  sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
  dmats <- createDataMatrices(data, formulas)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  m <- sapply(Y, length)
  Fi <- sv$Fi; Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u; Delta <- surv$Delta # survival
  l0 <- sv$l0
  S <- surv$S
  SS <- lapply(1:n, function(i){
    out <- apply(S[[i]], 2, rep, nrow(Fu[[i]]))
    if("numeric"%in%class(out)) out <- t(out)
    out
  })
  
  #' Parameter vector ----
  Omega <- list(D=D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, sigma, gamma, zeta)

  #' Gauss-Hermite Quadrature ----
  if(!is.null(control$gh.nodes)) gh <- control$gh.nodes else gh <- 3
  if(!is.null(control$gh.sigma)) .sigma <- control$gh.sigma else .sigma <- 1
  
  GH <- statmod::gauss.quad.prob(gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  #' Assign family to joint density and parameter updates ----
  if("function"%in%class(family)) family <- family()$family
  if(family == 'nbinom2') family <- 'negative.binomial' ## Catch-all.
  
  #' Begin EM ----
  diff <- 100; iter <- 0;
  if(!is.null(control$tol)) tol <- control$tol else tol <- 1e-2
  if(!is.null(control$maxit)) maxit <- control$maxit else maxit <- 100
  if(!is.null(control$conv)) conv <- control$conv else conv <- "relative"
  if(!conv%in%c('absolute', 'relative')) stop('Only "absolute" and "relative" difference convergence criteria implemented.')
  if(!is.null(control$verbose)) verbose <- control$verbose else verbose <- F
  
  EMstart <- proc.time()[3]
  while(diff > tol && iter < maxit){
    update <- EMupdate(Omega, family, X, Y, Z, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, sv, w, v, n, m)
    params.new <- setNames(c(vech(update$D), update$beta, update$sigma, update$gamma, update$zeta),
                           names(params))
    diff <- difference(params, params.new, conv)
    if(verbose){
      print(sapply(params.new, round, 4))
      message("Iteration ", iter + 1, ' maximum ', conv, ' difference: ', round(diff, 4))
    }

    # Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; sigma <- update$sigma
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    iter <- iter + 1
    Omega <- list(D=D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
    params <- params.new
    
  }
  EMend <- proc.time()[3]
  coeffs <- list(D = D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
  out <- list(coeffs = coeffs,
              RE = do.call(rbind, b),
              iter = iter,
              EMtime = EMend-EMstart,
              totaltime = proc.time()[3] - start.time)
  out$hazard <- cbind(ft = sv$ft, haz = l0)
  out$family <- family
  if(post.process){
    message('Post processing not yet implemented!')
  }
  
  out
}
