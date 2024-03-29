#' #######
#' EM.R
#' ----
#' Doesn't use a grid for quicker lookups for \lambda(\mu,\nu); logZ(\mu,\nu) and V(\mu, \nu) values.
#' ######

library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
source('_Functions.R')
source('vcov.R')
source('simData.R')
sourceCpp('nogrids.cpp')
source('inits.R')
vech <- function(x) x[lower.tri(x, T)]

EMupdate <- function(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                     sv, w, v, n, m, summax, debug, optimiser.arguments){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; delta <- Omega$delta; gamma <- Omega$gamma; zeta <- Omega$zeta
  
  #' Find b.hat and Sigma
  b.update <- b.minimise(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta, 
                         Omega, summax, method = optimiser.arguments$optimiser, 
                         obj = 'joint_density', Hessian = optimiser.arguments$Hessian, Hessian.eps = optimiser.arguments$eps)
  b.hat <- b.update$b.hat
  Sigma <- b.update$Sigma
  
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
  
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b.hat, SIMPLIFY = F)
  nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)

  #' lambda, logZ and V lookups ----
  lambdas <- mapply(function(mu, nu){
    lambda_appx(mu, nu, summax)
  }, mu = mus, nu = nus, SIMPLIFY = F)
  
  logZs <- mapply(function(mu, nu, lambda){
    logZ_c(log(lambda), nu, summax)
  }, mu = mus, nu = nus, lambda = lambdas, SIMPLIFY = F)
  
  Vs <- mapply(function(mu, nu, lambda, logZ){
    calc_V_vec(mu, lambda, nu, logZ, summax)
  }, mu = mus, nu = nus, lambda = lambdas, logZ = logZs, SIMPLIFY = F)
  
  # D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  # \beta
  Sb <- mapply(Sbeta, X, Y, mus, nus, lambdas, Vs, SIMPLIFY = F)
  Hb <- mapply(getW1, X, mus, nus, lambdas, Vs, SIMPLIFY = F)
  
  # \delta
  Sd <- mapply(function(G, b, X, Z, Y, lY, tau){
    Sdelta_cdiff(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps=.Machine$double.eps^(1/3))
  }, G = G, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)

  Hd <- mapply(function(G, b, X, Z, Y, lY, tau){
    Hdelta(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps=.Machine$double.eps^(1/4))
  }, G = G, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)
  
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
  (D.new <- Reduce('+', D.update)/n)
  
  # \beta and \delta
  (beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb)))
  (delta.new <- delta - solve(Reduce('+', Hd), c(Reduce('+', Sd))))

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
    D = D.new, beta = beta.new, delta = delta.new,       # <Y>
    gamma = gammazeta.new[1], zeta = gammazeta.new[-1],  # Survival
    l0 = l0.new, l0u = l0u.new, l0i = as.list(l0i.new),  # ---""---
    b = b.hat, mus = mus,                                #   REs.
    t = round(e-s,3)
  )
  
}

EM <- function(long.formula, disp.formula, surv.formula, data, summax = 100, post.process = T, 
               control = list(), disp.control = list(), delta.init = NULL, optim.control = list()){
  start.time <- proc.time()[3]
  
  #' Parsing formula objects ----
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- surv$n
  
  #' Initial conditions --> Longitudinal... ----
  inits.long <- Longit.inits(long.formula, disp.formula, data)
  
  #' Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  
  #' Data objects ----
  dmats <- createDataMatrices(data, formulas, disp.formula)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  lY <- lapply(Y, lfactorial)
  G <- dmats$G                             # Dispersion data matrix
  m <- sapply(Y, length)
  
  #' Unpack control args ----
  if(!is.null(control$tol)) tol <- control$tol else tol <- 1e-2
  if(!is.null(control$verbose)) verbose <- control$verbose else verbose <- F
  if(!is.null(control$debug)) debug <- control$debug else debug <- F
  if(!is.null(control$gh.nodes)) .gh <- control$gh.nodes else .gh <- 3
  if(!is.null(control$gh.sigma)) .sigma <- control$gh.sigma else .sigma <- 1
  if(!is.null(control$maxit)) maxit <- control$maxit else maxit <- 200
  if(!is.null(control$conv)) conv <- control$conv else conv <- "relative"
  if(!conv%in%c('absolute', 'relative')) stop('Only "absolute" and "relative" difference convergence criteria implemented.')
  if(!is.null(control$auto.summax)) auto.summax <- control$auto.summax else auto.summax <- T
  #' Control arguments specific to dispersion estimates ----
  if(!is.null(disp.control$delta.method)) delta.method <- disp.control$delta.method else delta.method <- 'optim'
  if(!delta.method %in% c('bobyqa', 'optim')) stop('delta.method must be either "optim" or "bobyqa".\n')
  if(!is.null(disp.control$min.profile.length)) min.profile.length <- disp.control$min.profile.length else min.profile.length <- 1
  if(!is.null(disp.control$what)) what <- disp.control$what else what <- 'median'
  if(!what %in% c('mean', 'median')) stop("what must be either 'mean' or 'median'.")
  if(!is.null(disp.control$cut)) cut <- disp.control$cut else cut <- T
  if(!is.null(disp.control$re.maximise)) re.maximise <- disp.control$re.maximise else re.maximise <- T
  if(!is.null(disp.control$max.val)) max.val <- disp.control$max.val else max.val <- 2
  
  if(auto.summax){
    summax.old <- summax
    summax <- max(sapply(Y, max)) * 2
    if(verbose) cat(sprintf("Automatically setting summax to %d\n", summax))
  }
  
  # Initial conditon for delta
  if(is.null(delta.init)){
    delta.inits.raw <- get.delta.inits(dmats, beta, b, delta.method, summax, verbose, min.profile.length, max.val)  
    # Return the user-specified estimate (Defaulting to cut + median)
    if(cut)
      initdelta <- if(what == 'mean') delta.inits.raw$mean.cut.estimate else delta.inits.raw$median.cut.estimate
    else
      initdelta <- if(what == 'mean') delta.inits.raw$mean.estimate else delta.inits.raw$median.estimate
    delta <- setNames(initdelta, names(inits.long$delta.init))
    if(verbose) cat(sprintf('Initial condition for delta: %.3f.\n', delta))
  }else{
    delta <- setNames(delta.init, names(inits.long$delta.init))
  }
  
  #' Set default optimiser arguments if _not_ specified.                      By default we...
  if(is.null(optim.control$optimiser)) optim.control$optimiser <- 'optim'             # Use optim by default (BFGS).
  if(is.null(optim.control$Hessian))   optim.control$Hessian <- 'obj'                 # Appraise 2nd deriv on objective function (not grad).
  if(is.null(optim.control$eps))       optim.control$eps <- .Machine$double.eps^(1/4) # Usual for fd.
  
  #' Re-maximise Marginal wrt b for CMP rather than Poisson (obtained thus far) ---------
  if(re.maximise){
    s <- proc.time()[3]
    if(verbose) cat('Re-maximising in light of new delta estimate')
    b <- b.minimise(b, X, Y, lY, Z, G, S = NULL, SS = NULL, Fi = NULL, Fu = NULL, l0i = NULL, l0u = NULL, Delta = NULL, 
                    Omega = list(beta = beta, delta = delta, D = D), summax, 
                    method = optim.control$optimiser, obj = 'marginal', optim.control$Hessian, optim.control$eps)$b.hat
    remaximisation.time <- round(proc.time()[3] - s, 3)
    if(verbose) cat(sprintf(", this took %.2f seconds.\n", remaximisation.time))
  }
  
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
  Omega <- list(D=D, beta = beta, delta = delta, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, delta, gamma, zeta)
  
  
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
    update <- EMupdate(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                       sv, w, v, n, m, summax, debug, optim.control)
    if(debug) DEBUG.update <<- update
    params.new <- c(vech(update$D), update$beta, update$delta, update$gamma, update$zeta)
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
  EMtime <- round(EMend - EMstart, 3)
  out <- list(coeffs = Omega,
              iter = iter)
  modelInfo <- list(
    forms = formulas,
    summax = summax
  )
  if(auto.summax) modelInfo$summax.type <- 'automatic' else modelInfo$summax.type <- 'manual'
  if(is.null(delta.init)){
    modelInfo$delta.init <- delta.inits.raw # Return ALL information.
  }else{
    modelInfo$delta.init <- delta.init
  }
  modelInfo$optimHess <- c(optimiser = optim.control$optimiser, Hessian = optim.control$Hessian, epsilon = optim.control$eps)
  out$modelInfo <- modelInfo
  out$hazard <- cbind(ft = sv$ft, haz = l0, nev = sv$nev)
  out$stepmat <- cbind(iter = 1:iter, time = step.times)
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time.p <- proc.time()[3]
    #' Calculating \b and \Sigma at MLEs
    b.update <- b.minimise(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta, 
                           Omega, summax, method = optim.control$optimiser, obj = 'joint_density', 
                           Hessian = optim.control$Hessian, optim.control$eps)
    b <- b.update$b.hat
    Sigma <- b.update$Sigma
  
    # The Information matrix
    I <- structure(vcov(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n, summax),
                   dimnames = list(names(params), names(params)))
    I.inv <- tryCatch(solve(I), error = function(e) e)
    if(inherits(I.inv, 'error')) I.inv <- structure(MASS::ginv(I),
                                                    dimnames = dimnames(I))
    out$SE <- sqrt(diag(I.inv))
    out$vcov <- I
    out$RE <- do.call(rbind, b)
    postprocess.time <- round(proc.time()[3]-start.time.p, 2)
  }
  comp.time <- round(proc.time()[3]-start.time, 3)
  
  out$elapsed.time <- c(`delta optimisation` = if(is.null(delta.init)) unname(delta.inits.raw$time) else NULL,
                        `re-maximisation` = if(re.maximise) unname(remaximisation.time) else NULL,
                        `startup time` = unname(startup.time),
                        `EM time` = unname(EMtime),
                        `post processing` = if(post.process) unname(postprocess.time) else NULL,
                        `Total computation time` = unname(comp.time))
  
  out
}
