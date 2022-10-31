#' #######
#' EM.R
#' ----
#' GP-1 implementation.
#' ######

library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
source('_Functions.R')
source('vcov.R')
source('simData.R')
sourceCpp('funs.cpp')
source('inits.R')
vech <- function(x) x[lower.tri(x, T)]

E_shape.b <- function(shape, X, Y, Z, tau, beta, b, w, v){
  shape <- shape
  eta <- X %*% beta + Z %*% b
  out <- numeric(length(w))
  for(l in 1:length(w)){
    mu.l <- exp(eta + tau * v[l])
    out[l] <- w[l] * ll_Gamma2(Y, shape, mu.l)
  }
  sum(out)
}

EMupdate <- function(Omega, X, Y, Z, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                     sv, w, v, n, m, debug){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; shape <- Omega$shape; gamma <- Omega$gamma; zeta <- Omega$zeta
  # qb <- mvtnorm::qmvnorm(.975, mean = 0, sigma = D, tail = 'lower.tail')$quantile
  
  b.hat <- mapply(function(b, X, Y, Z, S, SS, Fi, Fu, l0i, l0u, Delta){
    optim(b, joint_density, joint_density_ddb,
          X = X, Y = Y, Z = Z, beta = beta, shape = shape, D = D, S = S, SS = SS, Fi = Fi, 
          Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta, gamma = gamma, zeta = zeta,
          method = 'BFGS', hessian = T)
  }, b = b, X = X, Y = Y, Z = Z, S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, l0u = l0u, Delta = Delta,
  SIMPLIFY = F)
  Sigma <- lapply(b.hat, function(x) solve(x$hessian))
  b.hat <- lapply(b.hat, function(x) x$par)
  
  if(debug){DEBUG.b.hat <<- b.hat; DEBUG.Sigma <<- Sigma}
  check <- sapply(Sigma, det)
  if(any(check < 0 | is.nan(check))){
    ind <- which(check < 0 | is.nan(check))
    message('Some non pos-def or NaN Sigma for ids: ', paste(ind, collapse = ', '))
    for(i in seq_along(ind)){
      Sigma[[ind[i]]] <- as.matrix(Matrix::nearPD(Sigma[[ind[i]]], corr = T)$mat)
    }
  }
  
  #' #########
  #' E-step ##
  #' #########
  
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)
  
  # D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  #' \beta update...
  beta.update <- mapply(function(b, X, Y, Z){
    long_derivs(b = b, X = X, Y = Y, Z = Z, beta = beta, shape = shape, design = X)
  }, b = b.hat, X = X, Y = Y, Z = Z, SIMPLIFY = F)
  Sb <- lapply(beta.update, el, 1)
  Hb <- lapply(beta.update, el, 2)
  
  # Shape 
  shape.update <- mapply(function(b, X, Y, Z, tau){
    out <- setNames(vector('list', 2), c('Score', 'H'))
    out$Score <- pracma::grad(E_shape.b, shape, X = X, Y = Y, Z = Z, tau = tau, beta = beta, b = b, w = w, v = v)
    out$H <- pracma::hessian(E_shape.b, shape, X = X, Y = Y, Z = Z, tau = tau, beta = beta, b = b, w = w, v = v)
    out
  }, b = b.hat, X = X, Y = Y, Z = Z, tau = tau, SIMPLIFY = F)
  
  #' Survival parameters (\gamma, \zeta)
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
  (D.new <- Reduce('+', D.update)/n)
  
  # \beta and \phi
  (beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb)))
  # Shape
  (shape.new <- shape - solve(Reduce('+', lapply(shape.update, el, 2)),
                             Reduce('+', lapply(shape.update, el, 1))))

  # Survival parameters (gamma, zeta)
  (gammazeta.new <- c(gamma, zeta) - solve(Reduce('+', Hgz), rowSums(Sgz)))
  
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
    D = D.new, beta = beta.new, shape = shape.new,       # <Y>
    gamma = gammazeta.new[1], zeta = gammazeta.new[-1],  # Survival
    l0 = l0.new, l0u = l0u.new, l0i = as.list(l0i.new),  # ---""---
    b = b.hat,                                           #   REs.
    t = round(e-s,3)
  )
  
}

EM <- function(long.formula, surv.formula, data, post.process = T, control = list()){
  start.time <- proc.time()[3]
  
  #' Parsing formula objects ----
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- surv$n
  
  #' Initial conditions --> Longitudinal... ----
  inits.long <- Longit.inits(long.formula, data)
  
  #' Data objects ----
  dmats <- createDataMatrices(data, formulas)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  m <- sapply(Y, length)
  
  #' Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  #' If genpois then extract here, else set to length of dmats$w
  shape <- inits.long$shape.init
  
  #' Unpack control args ----
  if(!is.null(control$tol)) tol <- control$tol else tol <- 1e-2
  if(!is.null(control$verbose)) verbose <- control$verbose else verbose <- F
  if(!is.null(control$debug)) debug <- control$debug else debug <- F
  if(!is.null(control$gh.nodes)) .gh <- control$gh.nodes else .gh <- 3
  if(!is.null(control$gh.sigma)) .sigma <- control$gh.sigma else .sigma <- 1
  if(!is.null(control$maxit)) maxit <- control$maxit else maxit <- 200
  if(!is.null(control$conv)) conv <- control$conv else conv <- "relative"
  if(!conv%in%c('absolute', 'relative')) stop('Only "absolute" and "relative" difference convergence criteria implemented.')
  
  #' Initial conditions --> survival... ----
  inits.surv <- TimeVarCox(data, do.call(rbind, b), surv$ph, formulas)
  # Survival data objects ...
  sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
  Fi <- sv$Fi; Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u; Delta <- surv$Delta 
  l0 <- sv$l0
  S <- sv$S; SS <- sv$SS
  # Extract survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma', names(inits.surv$inits))]
   
  #' Parameter vector and list ----
  Omega <- list(D=D, beta = beta, shape = shape, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, shape, gamma, zeta)
  
  #' Gauss-hermite quadrature ----
  GH <- statmod::gauss.quad.prob(.gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  if(debug){DEBUG.inits <<- params; DEBUG.dmats <<- dmats; DEBUG.sv <<- sv; DEBUG.surv <<- surv}
  
  startup.time <- round(proc.time()[3] - start.time, 3)
  #' Begin EM ----
  diff <- 100; iter <- 0
  EMstart <- proc.time()[3]
  step.times <- c()
  if(verbose) message('Starting EM algorithm.')
  while(diff > tol && iter < maxit){
    update <- EMupdate(Omega, X, Y, Z, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                       sv, w, v, n, m, debug)
    if(debug) DEBUG.update <<- update
    params.new <- c(vech(update$D), update$beta, update$shape, update$gamma, update$zeta)
    names(params.new) <- names(params)
    
    # Convergence criteria + print (if wanted).
    diffs <- difference(params, params.new, conv)
    diff <- max(diffs)
    if(verbose){
      print(sapply(params.new, round, 4))
      message("Iteration ", iter + 1, ' maximum ', conv, ' difference: ', round(diff, 4), ' for parameter ', names(params.new)[which.max(diffs)])
      message("Largest RE relative difference: ", round(max(difference(do.call(cbind, update$b), do.call(cbind, b), conv)), 4))
    }
      
    #' Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; shape <- update$shape
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    iter <- iter + 1
    Omega <- list(D=D, beta = beta, shape = shape, gamma = gamma, zeta = zeta)
    params <- params.new
    step.times[iter] <- update$t # Store timings, as this could be interesting(?)
  }
  if(debug) DEBUG.Omega <<- Omega
  EMend <- proc.time()[3]
  EMtime <- round(EMend - EMstart, 3)
  out <- list(coeffs = Omega,
              iter = iter)
  modelInfo <- list(
    forms = c(formulas),
    names = list(beta = dmats$np, rand = dmats$nq),
    n = n, mi = m
  )
  out$modelInfo <- modelInfo
  out$hazard <- cbind(ft = sv$ft, haz = l0, nev = sv$nev)
  out$stepmat <- cbind(iter = 1:iter, time = step.times)
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time.p <- proc.time()[3]
    #' Calculating \b and \Sigma at MLEs
    b.hat <- mapply(function(b, X, Y, Z, S, SS, Fi, Fu, l0i, l0u, Delta){
      optim(b, joint_density, joint_density_ddb,
            X = X, Y = Y, Z = Z, beta = beta, shape = shape, D = D, S = S, SS = SS, Fi = Fi, 
            Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta, gamma = gamma, zeta = zeta,
            method = 'BFGS', hessian = T)
    }, b = b, X = X, Y = Y, Z = Z, S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, l0u = l0u, Delta = Delta,
    SIMPLIFY = F)
    Sigma <- lapply(b.hat, function(x) solve(x$hessian))
    b <- lapply(b.hat, function(x) x$par)
    # The Information matrix
    I <- structure(vcov(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n),
                   dimnames = list(names(params), names(params)))
    I.inv <- tryCatch(solve(I), error = function(e) e)
    if(inherits(I.inv, 'error')) I.inv <- structure(MASS::ginv(I),
                                                    dimnames = dimnames(I))
    out$SE <- sqrt(diag(I.inv))
    out$vcov <- I
    out$RE <- do.call(rbind, b)
    postprocess.time <- round(proc.time()[3]-start.time.p, 2)
    out$logLik <- log.lik(Omega, dmats, b, surv, sv, l0u, l0i)
  }
  comp.time <- round(proc.time()[3]-start.time, 3)
  
  out$elapsed.time <- c(`startup time` = unname(startup.time),
                        `EM time` = unname(EMtime),
                        `post processing` = if(post.process) unname(postprocess.time) else NULL,
                        `Total computation time` = unname(comp.time))
  
  out
}
