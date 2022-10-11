#' #######
#' EM.R                     time_var_nu/
#' ----
#' Treats dispersion as a **subject specific vector** of form (intercept, time).
#'                     (i.e. subjects variation in relation to mean value can change).
#' Subject specific truncation amounts, too, denoted `xi`.
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

EMupdate <- function(Omega, dmats, delta, b, sv, w, v, n, m, summax, debug, 
                     optimiser.arguments, inds.met, delta.update.quad, beta.update.quad, max.val, iter){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; gamma <- Omega$gamma.surv; zeta <- Omega$zeta;
  gamma.disp <- Omega$gamma.disp
  
  #' Find b.hat and Sigma
  b.update <- b.minimise(b, dmats, sv, delta, 
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
  
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), 
                X = dmats$X, Z = dmats$Z, b = b.hat, SIMPLIFY = F)
  nus <- mapply(function(W, delta) exp(W %*% delta), W = dmats$W, delta = delta, SIMPLIFY = F)
  if (beta.update.quad || delta.update.quad)
    tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = dmats$Z, SIMPLIFY = F)

  #' D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  #' \beta update...
  if(beta.update.quad){ #' Taken **with** quadrature
    Sb <- mapply(function(b, X, Z, W, Y, lY, delta, tau, summax){
      Sbeta2(beta, b, X, Z, W, Y, lY, delta, tau, w, v, summax)
    }, b = b.hat, X = dmats$X, Z = dmats$Z, W = dmats$W, Y = dmats$Y, 
       lY = dmats$lY, delta = delta, tau = tau, summax = summax, SIMPLIFY = F)
    
    Hb <- mapply(function(b, X, Z, W, Y, lY, delta, tau, summax){
      Hbeta2(beta, b, X, Z, W, Y, lY, delta, tau, w, v, summax)
    }, b = b.hat, X = dmats$X, Z = dmats$Z, W = dmats$W, Y = dmats$Y, lY = dmats$lY, 
       delta = delta, tau = tau, summax = summax, SIMPLIFY = F)
  }else{                #' Taken **without** quadrature
    Sb <- mapply(function(b, X, Z, W, Y, lY, delta, summax){
      Sbeta_noquad(beta, b, X, Z, W, Y, lY, delta, summax)
    }, b = b.hat, X = dmats$X, Z = dmats$Z, W = dmats$W, Y = dmats$Y, lY = dmats$lY, 
       delta = delta, summax = summax, SIMPLIFY = F)
    
    Hb <- mapply(function(b, X, Z, W, Y, lY, delta, summax){
      Hbeta_noquad(beta, b, X, Z, W, Y, lY, delta, summax)
    }, b = b.hat, X = dmats$X, Z = dmats$Z, W = dmats$W, Y = dmats$Y, lY = dmats$lY, 
       delta = delta, summax = summax, SIMPLIFY = F)
  }
  
  #' Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta, WFu, WFi, delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, 
               WFu, WFi, delta, gamma.disp, w, v, .Machine$double.eps^(1/3))
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, l0u = sv$l0u, Delta = sv$Delta,
     WFu = sv$WFu, WFi = sv$WFi, delta = delta)
  
  Hgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta, WFu, WFi, delta){
    Hgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, 
               WFu, WFi, delta, gamma.disp, w, v, 1e-3)#.Machine$double.eps^(1/4))
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, l0u = sv$l0u, Delta = sv$Delta,
  WFu = sv$WFu, WFi = sv$WFi, delta = delta, SIMPLIFY = F)
  
  #' #########
  #' M-step ##
  #' ######### 
  
  # D
  (D.new <- Reduce('+', D.update)/n)
  
  # \beta and \delta
  (beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb)))
  
  delta.new <- lapply(1:n, function(i){ # (All-in-one function!)
    # print(i)
    if(i %in% inds.met){
      # Appraise the CMP part (which dominates in practice...)
      lls <- delta.update(delta[[i]], dmats$X[[i]], dmats$Z[[i]], dmats$W[[i]], dmats$Y[[i]],
                   b.hat[[i]], beta, summax[[i]], w, v, tau[[i]], delta.update.quad)
      # Contribution of presence of delta in survival log-likelihood...
      svs <- delta_update(delta[[i]], sv$WFu[[i]], sv$WFi[[i]], sv$Delta[[i]],
                          sv$SS[[i]], sv$Fu[[i]], b.hat[[i]], gamma, gamma.disp, zeta,
                          sv$l0u[[i]], Sigma[[i]], w, v)
      return(delta[[i]] - solve(-lls$H + svs$Hessian, c(lls$Score + c(svs$Score))))
    }else{
      return(delta[[i]]) # Return pre-assigned rep(0, w). 
    } 
  })

  # Survival parameters (gamma, zeta)
  # inds.to.remove <- which(apply(do.call(rbind, delta), 1, function(i) any(abs(i) == max.val)))
  # qw <- length(c(gamma, gamma.disp, zeta))
  
  gamma.disp.update <- lapply(1:n, function(i){
    if(i %in% inds.met && abs(delta[[i]]) < max.val){
      disp_gamma_update(delta[[i]], b.hat[[i]], Sigma[[i]], sv$S[[i]], sv$SS[[i]],
                        sv$Fu[[i]], sv$Fi[[i]], sv$l0u[[i]], sv$Delta[[i]], sv$WFu[[i]], sv$WFi[[i]],
                        gamma.surv, gamma.disp, zeta, w, v)
    }else{
      return(list(Score = 0, Hessian = 0)) # Those at max value do __not__ contribute to update g_2
    }
  })
  
  # (gamma.disp.new <- gamma.disp - sum(sapply(gamma.disp.update, el, 1))/sum(sapply(gamma.disp.update, el, 2)))
  (gammazeta.new <- c(gamma, zeta) - solve(Reduce('+', Hgz), rowSums(Sgz)))
  if(iter >= 1){
    (gamma.disp.new <- gamma.disp - sum(sapply(gamma.disp.update, el, 1)) / sum(sapply(gamma.disp.update, el, 2)))
  }else{
    (gamma.disp.new <- gamma.disp)
  }
    
  gamma.surv.new <- gammazeta.new[1]; 
  zeta.new <- gammazeta.new[-1]

  
  # The baseline hazard and related objects
  lambda.update <- lambdaUpdate(sv$surv.times, sv$ft.mat, delta, sv$WFu, gamma, gamma.disp, zeta, sv$S, Sigma, b.hat, w, v)
  l0.new <- sv$nev/rowSums(lambda.update)
  l0u.new <- lapply(sv$l0u, function(ll){
    l0.new[1:length(ll)]
  })
  l0i.new <- l0.new[match(sv$Ti, sv$ft)] 
  l0i.new[is.na(l0i.new)] <- 0
  if(debug){plot(l0.new~sv$l0,pch=20);abline(0,1,col='red')}
  e <- proc.time()[3]
  # Return
  list(
    # Longit.
    D = D.new, beta = beta.new, delta = delta.new, b = b.hat,                              
    # Survival
    gamma.surv = gamma.surv.new, gamma.disp = gamma.disp.new,
    zeta = zeta.new,
    l0 = l0.new, l0u = l0u.new, l0i = as.list(l0i.new),  
    # Timing
    t = round(e-s,3)
  )  -> update
  
}

EM <- function(long.formula, surv.formula, disp.formula,
               data, post.process = T, 
               control = list(), disp.control = list(), optim.control = list(),
               delta.update.quad = T, beta.update.quad = F,
               summax.fn = NULL, min.summax = 20, initialise.delta = F){
  #' Defaults:
  #'   Truncation amount `xi`_i is taken as max(2 * max(Y_i), 20).
  #'   delta is not initialised for each subject
  #'   delta is updated using quadrature
  #'   beta is __not__ updated using quadrature.
  start.time <- proc.time()[3]
  
  #' Parsing formula objects ----
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data, disp.formula)
  n <- surv$n
  
  #' Initial conditions --> Longitudinal... ----
  inits.long <- Longit.inits(long.formula, data)
  
  #' Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  
  #' Data objects ----
  dmats <- createDataMatrices(data, formulas, disp.formula)
  m <- sapply(dmats$Y, length)
  
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
  if(!is.null(disp.control$min.profile.length)) min.profile.length <- disp.control$min.profile.length else min.profile.length <- 3
  if(!is.null(disp.control$max.val)) max.val <- disp.control$max.val else max.val <- 3
  if(!is.null(disp.control$truncated)) truncated <- disp.control$truncated else truncated <- F
  
  #' Set default optimiser arguments if _not_ specified.                      By default we...
  if(is.null(optim.control$optimiser)) optim.control$optimiser <- 'optim'             # Use optim by default (BFGS method).
  if(is.null(optim.control$Hessian))   optim.control$Hessian <- 'obj'                 # Appraise 2nd deriv on objective function (not `grad`).
  if(is.null(optim.control$eps))       optim.control$eps <- .Machine$double.eps^(1/4) # Usual for fd.
  
  #' Summax __for each subject__ (sometimes called `xi`).
  if(!is.null(summax.fn) & !'function'%in%class(summax.fn)) stop("Provide 'summax.fn' as a function.") 
  if(is.null(summax.fn)) summax.fn <- function(y) max(y) * 2
  summax <- as.list(pmax(min.summax, sapply(dmats$Y, summax.fn)))
  
  #' Indices of those who meet `min.profile.length` inclusion criterion (therefore treated as 'true' CMP)
  inds.met <- which(sapply(dmats$Y, function(y) length(unique(y))) > min.profile.length)
  inds.not <- which(sapply(dmats$Y, function(y) length(unique(y))) <= min.profile.length)
  (criteria.met <- sprintf("%.2f%%", length(inds.met)/n*100))
  
  #' An initial condition for delta ---------------------------------------
  delta <- lapply(1:n, function(i) setNames(rep(0, dmats$w), dmats$nw))
  # delta <- lapply(1:n, function(i){ # Quickly find delta at initial conditions (not yet conditioned on survival part)
  #   if(i %in% inds.met){            # and eschew quadrature for now...
  #     return(delta.update(delta[[i]],
  #                         dmats$X[[i]], dmats$Z[[i]], dmats$W[[i]], dmats$Y[[i]],
  #                         b[[i]], beta, summax[[i]], NULL, NULL, NULL, F)$new)
  #   }else{
  #     return(delta[[i]])
  #   }
  # })
  # delta <- lapply(delta, function(x) ifelse(abs(x) > max.val, sign(x) * max.val, x))
  
  #' Initial conditions --> survival... ----
  inits.surv <- TimeVarCox(data, do.call(rbind, b), surv$ph, formulas, disp.formula, delta, max.val)
  # Survival data objects ...
  sv <- surv.mod(surv$ph, surv$survdata, formulas, disp.formula, surv, inits.surv$l0.init, dmats$w)
  sv$Delta <- surv$Delta
  # Extract survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma.surv <- inits.surv$inits[grepl('gamma1$', names(inits.surv$inits))]
  gamma.disp <- c('gamma2' = 0)#inits.surv$inits[grepl('gamma2$', names(inits.surv$inits))]
   
  #' Parameter vector and list ----
  Omega <- list(D=D, beta = beta, gamma.surv = gamma.surv, gamma.disp = gamma.disp, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, gamma.surv, gamma.disp, zeta)
  
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
    update <- EMupdate(Omega, dmats, delta, b, sv, w, v, n, m, summax, debug, 
                       optim.control, inds.met, delta.update.quad, beta.update.quad, max.val, iter)
    if(debug) DEBUG.update <<- update
    params.new <- c(vech(update$D), update$beta, update$gamma.surv, update$gamma.disp, update$zeta)
    names(params.new) <- names(params)
    
    delta.old <- do.call(rbind, delta)
    delta <- do.call(rbind, update$delta)
    # Set back to max.val if necessary
    #if(truncated) 
    delta <- apply(delta, 2, function(d) ifelse(abs(d) > max.val, sign(d) * max.val, d))
    delta.diffs <- difference(delta.old, delta, 'relative')
    
    if(debug){
      par(mfrow = c(1, dmats$w))
      for(pp in 1:dmats$w){
        plot(delta.old[,pp], delta[,pp], 
             xlab = paste0('old delta_', dmats$nw[pp]), 
             ylab = paste0('new delta_', dmats$nw[pp]), pch = 20, main = bquote('Iteration'~.(iter+1)))
        abline(0, 1)
        abline(h = c(-1, 1) * max.val, col = 'red'); abline(v = c(-1, 1) * max.val, col = 'red')
      }
    }
    # Convergence criteria + print (if wanted).
    diffs <- difference(params, params.new, conv)
    diff <- max(diffs)
    if(verbose){
      print(sapply(params.new, round, 4))
      message("Iteration ", iter + 1, ' maximum ', conv, ' difference: ', round(diff, 4), ' for parameter ', names(params.new)[which.max(diffs)])
      message("Largest RE relative difference: ", round(max(difference(do.call(cbind, update$b), do.call(cbind, b), conv)), 4))
      if(iter > 0) message("Largest relative difference for subject-specific dispersion: ", round(max(delta.diffs), 4))
    }
      
    #' Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; delta <- lapply(1:n, function(i) delta[i,])
    gamma <- update$gamma.surv; gamma.disp <- update$gamma.disp; zeta <- update$zeta
    sv$l0 <- update$l0; sv$l0u <- update$l0u; sv$l0i <- update$l0i
    iter <- iter + 1
    Omega <- list(D=D, beta = beta, gamma.surv = gamma, gamma.disp = gamma.disp, zeta = zeta)
    params <- params.new
    step.times[iter] <- proc.time()[3] - s.em.i # Store timings, as this could be interesting?
  }
  if(debug) DEBUG.Omega <<- Omega
  EMend <- proc.time()[3]
  EMtime <- round(EMend - EMstart, 3)
  out <- list(coeffs = Omega,
              iter = iter)
  modelInfo <- list(
    forms = c(formulas, disp.formula = disp.formula),
    names = list(disp = dmats$nw, beta = dmats$np, rand = dmats$nq),
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
  out$hazard <- cbind(ft = sv$ft, haz = sv$l0, nev = sv$nev)
  out$stepmat <- cbind(iter = 1:iter, time = step.times)
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time.p <- proc.time()[3]
    #' Calculating \b and \Sigma at MLEs
    b.update <- b.minimise(b, dmats, sv, delta, 
                           Omega, summax, method = optim.control$optimiser, 
                           obj = 'joint_density', 
                           Hessian = optim.control$Hessian, Hessian.eps = optim.control$eps)
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
    delta.SEs <- as.data.frame(structure(do.call(rbind, invisible(lapply(inds.met, function(i){
      Idel <- vcv$Idelta[[i]]
      c(delta[[i]], suppressWarnings(sqrt(diag(solve(Idel)))))
    }))), dimnames = list(as.character(inds.met), c(dmats$nw, paste0('SE_', dmats$nw)))))
    
    delta.dfs <- setNames(lapply(1:dmats$w, function(w){
      estimate <- delta.SEs[,w]
      SE <- delta.SEs[,(dmats$w + w)] # We know structure of delta.SEs: (delta), (SE delta)
      SE <- ifelse(abs(estimate) == max.val | is.nan(SE), 0, SE)
      lb <- estimate - qnorm(.975) * SE; ub <- estimate + qnorm(.975) * SE
      data.frame(id = inds.met, 'estimate' = estimate, 'SE' = SE, 'lb' = lb, 'ub' = ub)
    }), dmats$nw)
    
    out$delta.df <- delta.dfs
    postprocess.time <- round(proc.time()[3]-start.time.p, 2)
    
    out$logLik <- log.lik(Omega, dmats, b, delta, surv, sv, l0u, l0i, summax)
    
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
