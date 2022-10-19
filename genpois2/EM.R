#' #######
#' EM.R
#' ----
#' TMB-esque implementation...
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

EMupdate <- function(Omega, X, Y, Z, W, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                     sv, w, v, n, m, debug, include.all, inds.met, delta.update.quad,
                     beta.update.quad, ww){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; delta <- Omega$delta; 
  gamma <- Omega$gamma; zeta <- Omega$zeta
  phi <- exp(delta)
  
  b.hat <- mapply(function(b, X, Y, Z, S, SS, Fi, Fu, l0i, l0u, Delta){
    optim(b, joint_density, joint_density_ddb,
          X = X, Y = Y, Z = Z, beta = beta, phi = phi, D = D, S = S, SS = SS, Fi = Fi, 
          Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta, gamma = gamma, zeta = zeta,
          method = 'BFGS', hessian = T)
  }, b = b, X = X, Y = Y, Z = Z,  S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, l0u = l0u, Delta = Delta,
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
  if(beta.update.quad){ #' Taken **with** quadrature
    # Sb <- mapply(function(b, X, Z, W, Y, lY, tau, summax){
    #   Sbeta2(beta, b, X, Z, W, Y, lY, delta, tau, w, v, summax)
    # }, b = b.hat, X = X, Z = Z, W = W, Y = Y, lY = lY, tau = tau, summax = summax, SIMPLIFY = F)
    # 
    # Hb <- mapply(function(b, X, Z, W, Y, lY, tau, summax){
    #   Hbeta2(beta, b, X, Z, W, Y, lY, delta, tau, w, v, summax)
    # }, b = b.hat, X = X, Z = Z, W = W, Y = Y, lY = lY, tau = tau, summax = summax, SIMPLIFY = F)
    stop('beta.update.quad = T not yet implemented!')
  }else{                #' Taken **without** quadrature
    beta.update <- mapply(function(b, X, Y, Z){
      long_derivs(b = b, X = X, Y = Y, Z = Z, beta = beta, phi = phi, design = X)
    }, b = b.hat, X = X, Y = Y, Z = Z, SIMPLIFY = F)
    Sb <- lapply(beta.update, el, 1)
    Hb <- lapply(beta.update, el, 2)
  }
  
  # \phi
  delta.updates <- lapply(1:n, function(i){ # (All-in-one function!)
    if(include.all){
      return(delta_update(b.hat[[i]], X[[i]], Y[[i]], Z[[i]], beta, delta,
                        w, v, tau[[i]]))
    }else if(i %in% inds.met){
      return(delta_update(b.hat[[i]], X[[i]], Y[[i]], Z[[i]], beta, delta,
                        w, v, tau[[i]]))
    }else{
      return(list(Score = rep(0, ww), Hessian = matrix(0, ww, ww))) # No contribution.
    }
  })
  
  
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
  (delta.new <- delta - solve(Reduce('+', lapply(delta.updates, el, 2)), 
                            c(Reduce('+', lapply(delta.updates, el, 1)))))
  
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
    b = b.hat,                                           #   REs.
    t = round(e-s,3)
  ) #-> update
  
}

EM <- function(long.formula, disp.formula = ~1, surv.formula, data, summax = NULL, post.process = T, 
               control = list(), optim.control = list(), genpois.inits = F){
  start.time <- proc.time()[3]
  
  #' Parsing formula objects ----
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- surv$n
  
  #' Initial conditions --> Longitudinal... ----
  inits.long <- Longit.inits(long.formula, disp.formula, data, genpois.inits = genpois.inits)
  
  #' Data objects ----
  dmats <- createDataMatrices(data, formulas, disp.formula)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  lY <- lapply(Y, lfactorial)
  W <- dmats$W                             # Dispersion data matrix
  m <- sapply(Y, length)
  
  #' Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  #' If genpois then extract here, else set to length of dmats$w
  if(genpois.inits){
    delta <- inits.long$delta
  }else{
    delta <- setNames(rep(0, dmats$w), dmats$nw)
  }
  
  #' Unpack control args ----
  if(!is.null(control$tol)) tol <- control$tol else tol <- 1e-2
  if(!is.null(control$verbose)) verbose <- control$verbose else verbose <- F
  if(!is.null(control$debug)) debug <- control$debug else debug <- F
  if(!is.null(control$gh.nodes)) .gh <- control$gh.nodes else .gh <- 3
  if(!is.null(control$gh.sigma)) .sigma <- control$gh.sigma else .sigma <- 1
  if(!is.null(control$maxit)) maxit <- control$maxit else maxit <- 200
  if(!is.null(control$conv)) conv <- control$conv else conv <- "relative"
  if(!conv%in%c('absolute', 'relative')) stop('Only "absolute" and "relative" difference convergence criteria implemented.')
  if(!is.null(control$min.profile.length)) min.profile.length <- control$min.profile.length else min.profile.length <- 3
  if(!is.null(control$include.all)) include.all <- control$include.all else include.all <- T
  if(!is.null(control$delta.update.quad)) delta.update.quad <- control$delta.update.quad else delta.update.quad <- T
  if(!is.null(control$beta.update.quad)) beta.update.quad <- control$beta.update.quad else beta.update.quad <- F
  
  #' Indices of those who meet `min.profile.length` inclusion criterion (therefore treated as true genpois)
  #' Only these subjects contribute to update for \delta.
  inds.met <- which(sapply(Y, function(y) length(unique(y))) > min.profile.length)
  inds.not <- which(sapply(Y, function(y) length(unique(y))) <= min.profile.length)
  (criteria.met <- sprintf("%.2f%%", length(inds.met)/n*100))
  
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
    update <- EMupdate(Omega, X, Y, Z, W, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                       sv, w, v, n, m, debug, include.all,
                       inds.met, delta.update.quad, beta.update.quad, dmats$w)
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
    forms = c(formulas, disp.formula = disp.formula),
    names = list(disp = dmats$nw, beta = dmats$np, rand = dmats$nq),
    n = n, mi = m,
    inds.met = inds.met,
    include.all = include.all,
    delta.update.quad = delta.update.quad,
    beta.update.quad = beta.update.quad
  )
  out$modelInfo <- modelInfo
  out$hazard <- cbind(ft = sv$ft, haz = l0, nev = sv$nev)
  out$stepmat <- cbind(iter = 1:iter, time = step.times)
  
  if(post.process){
    message('\nCalculating SEs...')
    phi <- exp(delta)
    start.time.p <- proc.time()[3]
    #' Calculating \b and \Sigma at MLEs
    b <- mapply(function(b, X, Y, Z, S, SS, Fi, Fu, l0i, l0u, Delta){
      optim(b, joint_density, joint_density_ddb,
            X = X, Y = Y, Z = Z, beta = beta, phi = phi, D = D, S = S, SS = SS, Fi = Fi, 
            Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta, gamma = gamma, zeta = zeta,
            method = 'BFGS', hessian = T)
    }, b = b, X = X, Y = Y, Z = Z,  S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, l0u = l0u, Delta = Delta,
    SIMPLIFY = F)
    Sigma <- lapply(b, function(x) solve(x$hessian))
    b <- lapply(b, function(x) x$par)
    # The Information matrix
    I <- structure(vcov(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n,
                        include.all, inds.met, delta.update.quad, beta.update.quad),
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
