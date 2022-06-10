#' #######
#' EM.R
#' --
#' This largely based on my PBC-case-study/EM, inits, _Functions.R files
#' ######

library(survival)
library(glmmTMB)
library(Rcpp)
library(RcppArmadillo)
source('_Functions.R')
source('simData.R')
source('inits.R')
source('vcov.R')
sourceCpp('cmp.cpp')
vech <- function(x) x[lower.tri(x, T)]

EMupdate <- function(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                     sv, w, v, n, m, summax, debug){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; delta <- Omega$delta; gamma <- Omega$gamma; zeta <- Omega$zeta
  
  #' Find b.hat and Sigma
  b.hat <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
    optim(b, joint_density, joint_density_ddb,
          X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
          S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
          gamma = gamma, zeta = zeta, summax = summax, method = 'BFGS')$par
  }, b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
  l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
    solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
                            S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                            gamma = gamma, zeta = zeta, summax = summax, eps = 1e-3))
  }, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
  l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  if(debug){DEBUG.b.hat <<- b.hat; DEBUG.Sigma <<- Sigma}
  #' #########
  #' E-step ##
  #' #########
  
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b.hat, SIMPLIFY = F)
  nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)
  lambdas <- mapply(function(mu, nu) lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax), mu = mus, nu = nus, SIMPLIFY = F)
  V <- mapply(V_mu_lambda, mus, lambdas, nus, summax = summax, SIMPLIFY = F)
  ABC <- mapply(calc.ABC, mus, nus, lambdas, summax, SIMPLIFY = F)
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), Z = Z, S = Sigma, SIMPLIFY = F)
  # ABC <- mapply(function(mu, nu, lambda, tau){
  #   suppressWarnings(calc2.ABC(mu, nu, lambda, summax, tau, w, v))
  # }, mu = mus, nu = nus, lambda = lambdas, tau = tau, SIMPLIFY = F)
  # D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  # \beta
  Sb <- mapply(Sbeta, X, Y, mus, nus, lambdas, V, SIMPLIFY = F)
  Hb <- mapply(getW1, X, mus, nus, lambdas, V, SIMPLIFY = F)
  # Sb <- mapply(function(X, Y, Z, G, b){
  #   Sbeta(beta, X, Y, Z, G, b, delta, summax = summax)
  # }, X = X, Y = Y, Z = Z, G = G, b = b.hat, SIMPLIFY = F)
  # Hb <- mapply(function(X, Z, G, b){
  #   getW1(X, Z, G, b, beta, delta, summax)
  # }, X = X, Z = Z, G = G, b = b.hat, SIMPLIFY = F)
  # Hb <- mapply(function(X, Y, Z, G, b){
  #   GLMMadaptive:::fd_vec(beta, Sbeta, X, Y, Z, G, b, delta, summax = summax)
  # }, X = X, Y = Y, Z = Z, G = G, b = b.hat, SIMPLIFY = F)
  
  # \delta
  # Sd <- mapply(function(ABC, Y, mu, V, nu, G){
  #   crossprod(((ABC$A * (Y - mu) / V - lgamma(Y + 1) + ABC$B) * nu), G)
  # }, ABC = ABC, Y = Y, mu = mus, V = V, nu = nus, G = G)
  Sd <- mapply(function(ABC, Y, mu, V, nu, G, tau){
    lhs <- numeric(length(mu))
    for(l in 1:length(w)) lhs <- lhs + w[l] * ABC$A * (Y - mu * exp(tau * v[l])) / V
    crossprod(((lhs - lgamma(Y + 1) + ABC$B) * nu), G)
  }, ABC = ABC, Y = Y, mu = mus, V = V, nu = nus, G = G, tau= tau, SIMPLIFY = F)
  Hd <- mapply(getW2, ABC, V, nus, G, SIMPLIFY = F)
  # tau <- mapply(function(Z, S) sqrt(diag(tcrossprod(Z %*% S, Z))), Z = Z, S = Sigma, SIMPLIFY = F)
  # Sd <- mapply(function(X, Y, lY, Z, b, G){
  #   Sdelta(delta, X, Y, lY, Z, b, G, beta, summax)
  # }, X = X, Y = Y, lY = lY, Z = Z, b = b.hat, G = G, SIMPLIFY = F)
  # Hd <- mapply(function(X, Y, lY, Z, b, G){
  #   GLMMadaptive:::fd_vec(delta, Sdelta, X, Y, lY, Z, b, G, beta, summax)
  # }, X = X, Y = Y, lY = lY, Z = Z, b = b.hat, G = G, SIMPLIFY = F)
  
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
  
  # \beta and \delta
  beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
  delta.new <- delta - solve(Reduce('+', Hd), Reduce('+', Sd))

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
  if(debug){plot(l0.new~l0,pch=20);abline(0,1,col='red')}
  e <- proc.time()[3]
  # Return
  list(
    D = D.new, beta = beta.new, delta = delta.new,       # <Y>
    gamma = gammazeta.new[1], zeta = gammazeta.new[-1],  # Survival
    l0 = l0.new, l0u = l0u.new, l0i = as.list(l0i.new),  # ---""---
    b = b.hat,                                           #   REs.
    t = round(e-s,3)
  )
  
}

EM <- function(long.formula, disp.formula, surv.formula, data, summax = 100, post.process = T, control = list()){
  start.time <- proc.time()[3]
  
  #' Parsing formula objects ----
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- surv$n
  
  #' Initial conditions ----
  inits.long <- Longit.inits(long.formula, disp.formula, data)
  inits.surv <- TimeVarCox(data, inits.long$b, surv$ph, formulas)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  delta <- inits.long$delta.init
  # Survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma', names(inits.surv$inits))]
  
  #' Data objects ----
  sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
  dmats <- createDataMatrices(data, formulas, disp.formula)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  lY <- lapply(Y, lfactorial)
  G <- dmats$G                             # Dispersion data matrix
  m <- sapply(Y, length)
  # survival
  Fi <- sv$Fi; Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u; Delta <- surv$Delta 
  l0 <- sv$l0
  S <- sv$S; SS <- sv$SS
  
  #' Parameter vector and list ----
  Omega <- list(D=D, beta = beta, delta = delta, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, delta, gamma, zeta)
  
  #' Unpack control args ----
  if(!is.null(control$tol)) tol <- control$tol else tol <- 1e-2
  if(!is.null(control$verbose)) verbose <- control$verbose else verbose <- F
  if(!is.null(control$debug)) debug <- control$debug else debug <- F
  if(!is.null(control$gh.nodes)) .gh <- control$gh.nodes else .gh <- 3
  if(!is.null(control$gh.sigma)) .sigma <- control$gh.sigma else .sigma <- 1
  if(!is.null(control$maxit)) maxit <- control$maxit else maxit <- 200
  if(!is.null(control$conv)) conv <- control$conv else conv <- "relative"
  if(!conv%in%c('absolute', 'relative')) stop('Only "absolute" and "relative" difference convergence criteria implemented.')
  if(!is.null(control$summax.override)) summax.override <- control$summax.override else summax.override <- F
  
  
  if(debug){DEBUG.inits <<- params; DEBUG.dmats <<- dmats; DEBUG.sv <<- sv; DEBUG.surv <<- surv}
  #' Gauss-hermite quadrature ----
  GH <- statmod::gauss.quad.prob(.gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  #' Overwrite summax if too low ----
  summax.old <- summax
  summax.new <- max(ceiling(c(max(data[, formulas$response]) + 20 * sqrt(var(data[, formulas$response])))), 100)
  if(!summax.override){
    if(summax.new>summax.old){
      if(verbose){
        cat(sprintf('Original summax %d too low, running with new value --> %d.\n', summax.old, summax.new))
      }
    summax <- summax.new
    }
  }else{
    summax <- summax.old
  }
  
  #' Begin EM ----
  diff <- 100; iter <- 0
  EMstart <- proc.time()[3]
  step.times <- c()
  while(diff > tol && iter < maxit){
    update <- EMupdate(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                       sv, w, v, n, m, summax, debug)
    if(debug) DEBUG.update <<- update
    params.new <- c(vech(update$D), update$beta, update$delta, update$gamma, update$zeta)
    names(params.new) <- names(params)
    # Convergence criteria + print (if wanted).
    diff <- difference(params, params.new, conv)
    if(verbose){
      print(sapply(params.new, round, 4))
      message("Iteration ", iter + 1, ' maximum ', conv, ' difference: ', round(diff, 4))
    }
    
    #' Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; delta <- update$delta
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    iter <- iter + 1
    Omega <- list(D=D, beta = beta, delta = delta, gamma = gamma, zeta = zeta)
    params <- params.new
    step.times[iter] <- update$t # Store timings, as this could be interesting(?)
  }
  if(debug) DEBUG.Omega <<- Omega
  EMend <- proc.time()[3]
  out <- list(coeffs = Omega,
              EMtime = round(EMend - EMstart, 3),
              iter = iter,
              comp.time = round(proc.time()[3]-start.time, 3))
  modelInfo <- list(
    summax = summax, 
    forms = formulas
  )
  out$modelInfo <- modelInfo
  out$hazard <- cbind(ft = sv$ft, haz = l0)
  out$stepmat <- cbind(iter = 1:iter, time = step.times)
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time <- proc.time()[3]
    #' Calculating \b and \Sigma at MLEs
    b.hat <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
      optim(b, joint_density, joint_density_ddb,
            X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
            S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
            gamma = gamma, zeta = zeta, summax = summax, method = 'BFGS')$par
    }, b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
    l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
    
    Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
      solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
                              S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                              gamma = gamma, zeta = zeta, summax = summax, eps = 1e-3))
    }, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
    l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
    # The Information matrix
    I <- structure(vcov(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n, summax),
                   dimnames = list(names(params), names(params)))
    I.inv <- tryCatch(solve(I), error = function(e) e)
    if(inherits(I.inv, 'error')) I.inv <- structure(MASS::ginv(I),
                                                    dimnames = dimnames(I))
    out$SE <- sqrt(diag(I.inv))
    out$vcov <- I
    out$RE <- do.call(rbind, b)
    out$postprocess.time <- round(proc.time()[3]-start.time, 2)
    
    #' Calculate the log-likelihood.
    # Timing done separately as EMtime + postprocess.time is for EM + SEs.
    out$logLik <- log.lik(Omega, dmats, b, surv, sv, l0u, l0i, summax)
    out$lltime <- round(proc.time()[3] - start.time, 2) - out$postprocess.time
  }
  
  out
}