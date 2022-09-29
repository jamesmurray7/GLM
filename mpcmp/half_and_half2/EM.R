#' #######
#' EM.R                     half_and_half2/
#' ----
#' Doesn't use a grid for quicker lookups for \lambda(\mu,\nu); logZ(\mu,\nu) and V(\mu, \nu) values.
#' Treats dispersion as a **subject specific scalar**.
#' Treats truncation amounts sub-specific, too.
#' ######

library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
source('_Functions.R')
source('vcov.R')
source('simData.R')
sourceCpp('halfhalf.cpp')
source('inits.R')
vech <- function(x) x[lower.tri(x, T)]

EMupdate <- function(Omega, X, Y, lY, Z, delta, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                     sv, w, v, n, m, summax, debug, optimiser.arguments, inds.met, delta.update.quad, beta.update.quad){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; gamma <- Omega$gamma; zeta <- Omega$zeta
  
  #' Find b.hat and Sigma
  b.update <- b.minimise(b, X, Y, lY, Z, delta, S, SS, Fi, Fu, l0i, l0u, Delta, 
                         Omega, summax, method = optimiser.arguments$optimiser, 
                         obj = 'joint_density', 
                         Hessian = optimiser.arguments$Hessian, Hessian.eps = optimiser.arguments$eps)
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
  nus <- mapply(function(delta, Y) rep(exp(delta), length(Y)), delta = delta, Y = Y, SIMPLIFY = F)
  if (beta.update.quad || delta.update.quad)
    tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)

  #' D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  #' \beta update...
  if(beta.update.quad){ #' Taken **with** quadrature
    Sb <- mapply(function(b, X, Z, Y, lY, delta, tau, summax){
      Sbeta2(beta, b, X, Z, Y, lY, delta, tau, w, v, summax)
    }, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, delta = delta, tau = tau, summax = summax, SIMPLIFY = F)
    
    Hb <- mapply(function(b, X, Z, Y, lY, delta, tau, summax){
      Hbeta2(beta, b, X, Z, Y, lY, delta, tau, w, v, summax)
    }, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, delta = delta, tau = tau, summax = summax, SIMPLIFY = F)
  }else{                #' Taken **without** quadrature
    Sb <- mapply(function(b, X, Z, Y, lY, delta, summax){
      Sbeta_noquad(beta, b, X, Z, Y, lY, delta, summax)
    }, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, delta = delta, summax = summax, SIMPLIFY = F)
    
    Hb <- mapply(function(b, X, Z, Y, lY, delta, summax){
      Hbeta_noquad(beta, b, X, Z, Y, lY, delta, summax)
    }, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, delta = delta, summax = summax, SIMPLIFY = F)
  }
  
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
  
  # \beta and \delta
  (beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb)))
  
  delta.new <- lapply(1:n, function(i){ # (All-in-one function!)
    if(i %in% inds.met)
      return(delta.update(delta[[i]], X[[i]], Z[[i]], Y[[i]],
                          b.hat[[i]], beta, summax[[i]], w, v, tau[[i]], delta.update.quad)$new)
    else 
      return(0)
  })

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
    b = b.hat,                                #   REs.
    t = round(e-s,3)
  ) # -> update
  
}

EM <- function(long.formula, surv.formula, data, post.process = T, 
               control = list(), disp.control = list(), optim.control = list(),
               delta.update.quad = T, beta.update.quad = T,
               summax.fn = NULL, min.summax = 20, initialise.delta = T){
  #' Defaults: 
  #'   _not_ updating \delta, and doing so every iteration in in instances where it `update.deltas` is `TRUE`.
  #'   Truncation amount `summax` is taken as max(2 * max(Y_i), 20).
  start.time <- proc.time()[3]
  
  #' Parsing formula objects ----
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- surv$n
  
  #' Initial conditions --> Longitudinal... ----
  inits.long <- Longit.inits(long.formula, data)
  
  #' Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  
  #' Data objects ----
  dmats <- createDataMatrices(data, formulas)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  lY <- lapply(Y, lfactorial)
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
  #' Control arguments specific to dispersion estimates ----
  if(!is.null(disp.control$delta.method)) delta.method <- disp.control$delta.method else delta.method <- 'bobyqa'
  if(!delta.method %in% c('bobyqa', 'optim')) stop('delta.method must be either "optim" or "bobyqa".\n')
  if(!is.null(disp.control$min.profile.length)) min.profile.length <- disp.control$min.profile.length else min.profile.length <- 3
  if(!is.null(disp.control$re.maximise)) re.maximise <- disp.control$re.maximise else re.maximise <- T
  if(!is.null(disp.control$max.val)) max.val <- disp.control$max.val else max.val <- 3
  if(!is.null(disp.control$truncated)) truncated <- disp.control$truncated else truncated <- F
  
  #' Set default optimiser arguments if _not_ specified.                      By default we...
  if(is.null(optim.control$optimiser)) optim.control$optimiser <- 'optim'             # Use optim by default (BFGS method).
  if(is.null(optim.control$Hessian))   optim.control$Hessian <- 'obj'                 # Appraise 2nd deriv on objective function (not `grad`).
  if(is.null(optim.control$eps))       optim.control$eps <- .Machine$double.eps^(1/4) # Usual for fd.
  
  #' Summax __for each subject__ (sometimes called `xi`).
  if(!is.null(summax.fn) & !'function'%in%class(summax.fn)) stop("Provide 'summax.fn' as a function.") 
  if(is.null(summax.fn)) summax.fn <- function(y) max(y) * 2
  summax <- as.list(pmax(min.summax, sapply(Y, summax.fn)))
  
  # Initial conditon for deltas
  if(initialise.delta){  #' Initialise \delta_i if requested (default)
    delta.inits.raw <- get.delta.inits(dmats, beta, b, delta.method, summax, verbose, min.profile.length, max.val)  
    # delta <- if(truncated) as.list(delta.inits.raw$truncated.estimates) else as.list(delta.inits.raw$subject.estimates)
    delta <- sapply(delta.inits.raw$subject.estimates,
                    function(d) if(abs(d) > max.val || (max.val - abs(d)) <= 1e-5) sign(d) else d)
    
    #' Re-maximise Marginal wrt b for CMP rather than Poisson (obtained thus far) ---------
    if(re.maximise){
      s <- proc.time()[3]
      if(verbose) cat('Re-maximising in light of new delta estimates')
      b <- b.minimise(b, X, Y, lY, Z, delta, S = NULL, SS = NULL, Fi = NULL, Fu = NULL, l0i = NULL, l0u = NULL, Delta = NULL, 
                      Omega = list(beta = beta, D = D), summax, 
                      method = optim.control$optimiser, obj = 'marginal', optim.control$Hessian, optim.control$eps)$b.hat
      remaximisation.time <- round(proc.time()[3] - s, 3)
      if(verbose) cat(sprintf(", this took %.2f seconds.\n", remaximisation.time))
    }
  }else{                 #' Otherwise set delta_i = 0 \forall i = 1, ..., n.
    delta <- as.list(rep(0, n))
  }
  
  #' Indices of those who meet `min.profile.length` inclusion criterion (therefore treated as true CMP)
  inds.met <- which(sapply(Y, function(y) length(unique(y))) > min.profile.length)
  inds.not <- which(sapply(Y, function(y) length(unique(y))) <= min.profile.length)
  criteria.met <- sprintf("%.2f%%", length(inds.met)/n*100)
  
  
  # plot(sapply(1:n, function(i) b[[i]] - inits.long$b[i,]))
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
  Omega <- list(D=D, beta = beta, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, gamma, zeta)
  
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
    s.em.i <- proc.time()[3]
    update <- EMupdate(Omega, X, Y, lY, Z, delta, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                       sv, w, v, n, m, summax, debug, optim.control, inds.met, delta.update.quad, beta.update.quad)
    if(debug) DEBUG.update <<- update
    params.new <- c(vech(update$D), update$beta, update$gamma, update$zeta)
    names(params.new) <- names(params)
    
    delta.old <- unlist(delta)
    delta <- unlist(update$delta)
    # Set back to max.val if necessary
    #if(truncated) 
    delta <- sapply(delta, function(d) if(abs(d) > max.val) return(sign(d) * max.val) else return(d))
    delta.diffs <- difference(delta.old, delta, 'relative')
    
    if(debug){
      plot(delta.old, delta, xlab = 'old', ylab = 'new', pch = 20, main = bquote('Iteration'~.(iter+1)))
      abline(0, 1)
      abline(h = c(-1, 1) * max.val, col = 'red'); abline(v = c(-1, 1) * max.val, col = 'red')
    }
    # Convergence criteria + print (if wanted).
    diffs <- difference(params, params.new, conv)
    diff <- max(diffs)
    if(verbose){
      print(sapply(params.new, round, 4))
      message("Iteration ", iter + 1, ' maximum ', conv, ' difference: ', round(diff, 4), ' for parameter ', names(params.new)[which.max(diffs)])
      message("Largest RE relative difference: ", round(max(difference(do.call(cbind, update$b), do.call(cbind, b), conv)), 4))
      message("Largest relative difference for subject-specific dispersion: ", round(max(delta.diffs), 4))
    }
      
    #' Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; delta <- as.list(delta)
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    iter <- iter + 1
    Omega <- list(D=D, beta = beta, gamma = gamma, zeta = zeta)
    params <- params.new
    step.times[iter] <- proc.time()[3] - s.em.i # Store timings, as this could be interesting?
  }
  if(debug) DEBUG.Omega <<- Omega
  EMend <- proc.time()[3]
  EMtime <- round(EMend - EMstart, 3)
  out <- list(coeffs = Omega,
              iter = iter)
  modelInfo <- list(
    forms = formulas,
    summax = unlist(summax),
    n = n, mi = m,
    inds.met = inds.met,
    max.val = max.val,
    delta.update.quad = delta.update.quad,
    beta.update.quad = beta.update.quad,
    summax.fn = summax.fn,
    min.summax = min.summax
  )
  if(initialise.delta) modelInfo$delta.init <- delta.inits.raw # Return ALL information.
  modelInfo$optimHess <- c(optimiser = optim.control$optimiser, Hessian = optim.control$Hessian, epsilon = optim.control$eps)
  out$modelInfo <- modelInfo
  out$hazard <- cbind(ft = sv$ft, haz = l0, nev = sv$nev)
  out$stepmat <- cbind(iter = 1:iter, time = step.times)
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time.p <- proc.time()[3]
    #' Calculating \b and \Sigma at MLEs
    b.update <- b.minimise(b, X, Y, lY, Z, delta, S, SS, Fi, Fu, l0i, l0u, Delta, 
                           Omega, summax, method = optim.control$optimiser, obj = 'joint_density', 
                           Hessian = optim.control$Hessian, optim.control$eps)
    b <- b.update$b.hat
    Sigma <- b.update$Sigma
  
    # The Information matrix
    vcv <- vcov(Omega, delta, dmats, surv, sv, Sigma, b, l0u, w, v, n, summax, inds.met, delta.update.quad,
                beta.update.quad)
    I <- structure(vcv$I, dimnames = list(names(params), names(params)))
    I.inv <- tryCatch(solve(I), error = function(e) e)
    if(inherits(I.inv, 'error')) I.inv <- structure(MASS::ginv(I),
                                                    dimnames = dimnames(I))
    out$SE <- sqrt(diag(I.inv))
    out$vcov <- I
    out$RE <- do.call(rbind, b)
    delta.df <- data.frame(
      id = inds.met,
      delta = unlist(delta)[inds.met],
      SE  = sqrt(1/vcv$Idelta)[inds.met],
      truncated = ifelse(abs(unlist(delta)[inds.met]) == max.val, '*', '')
    )
    out$delta.df <- delta.df
    postprocess.time <- round(proc.time()[3]-start.time.p, 2)
    
  }
  comp.time <- round(proc.time()[3]-start.time, 3)
  
  out$elapsed.time <- c(`delta optimisation` = if(initialise.delta) unname(delta.inits.raw$time) else NULL,
                        `re-maximisation` = if(re.maximise & initialise.delta) unname(remaximisation.time) else NULL,
                        `startup time` = unname(startup.time),
                        `EM time` = unname(EMtime),
                        `post processing` = if(post.process) unname(postprocess.time) else NULL,
                        `Total computation time` = unname(comp.time))
  
  out
}
