#' #######
#' EM.R                 ((HYBRID))
#' ----
#' At each iteration, constructs grid, or adds to existing grid...
#' ######

# rm(list=ls())
library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
source('_Functions.R')
source('simData.R')
source('../inits.R')
source('vcov.R')
sourceCpp('test.cpp')
vech <- function(x) x[lower.tri(x, T)]

EMupdate <- function(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                     sv, w, v, num, m, summax, debug, all.mus, nu.vec, lambda.mat, logZ.mat, V.mat){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; delta <- Omega$delta; gamma <- Omega$gamma; zeta <- Omega$zeta
  
  #' Find b.hat and Sigma
  b.hat <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
    optim(b, joint_density, joint_density_ddb,
          X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
          S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
          gamma = gamma, zeta = zeta, 
          lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, 
          all_mus = all.mus, all_nus = round(nu.vec, 3),
          summax = summax, method = 'BFGS')$par
  }, b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
  l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)


  Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
    solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
                            S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                            gamma = gamma, zeta = zeta, lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, 
                            all_mus = all.mus, all_nus = round(nu.vec, 3), summax = summax, eps = .001))
  }, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
  l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  if(debug){DEBUG.b.hat <<- b.hat; DEBUG.Sigma <<- Sigma}
  
  check <- sapply(Sigma, det)
  if(any(check < 0)){
    message('Some non pos def Sigma...')
    ind <- which(check<0)
    for(i in seq_along(ind)){
      Sigma[[ind[i]]] <- matrix(0, length(b.hat[[ind[i]]]), length(b.hat[[ind[i]]]))
    }
  }
  #' #########
  #' E-step ##
  #' #########
  
  #' NEW \mus, \nus, and calculate \tau
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b.hat, SIMPLIFY = F)
  nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)
  
  #' Indices for lookup given new \b.hat
  m <- mapply(function(a) get_indices(a, all.mus), a = mus, SIMPLIFY = F)
  n <- mapply(function(a) get_indices(a, round(nu.vec, 3)), a = nus, SIMPLIFY = F)
  
  #' \lambdas, Vs for \beta update.
  lambdas <- mapply(function(m, n) mat_lookup(m, n, lambda.mat), m = m, n = n, SIMPLIFY = F)
  Vs <- mapply(function(m, n) mat_lookup(m, n, V.mat), m = m, n = n, SIMPLIFY = F)
  
  #' D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  #' \beta
  Sb <- mapply(Sbeta, X, Y, mus, nus, lambdas, Vs, SIMPLIFY = F)
  Hb <- mapply(getW1, X, mus, nus, lambdas, Vs, SIMPLIFY = F)
  
  #' \delta
  Sd <- mapply(function(G, b, X, Z, Y, lY, tau){
    Sdelta_cdiff(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps=.Machine$double.eps^(1/3))
  }, G = G, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)
  
  Hd <- mapply(function(G, b, X, Z, Y, lY, tau){
    Hdelta(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps=.Machine$double.eps^(1/4))
  }, G = G, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)
  
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
  (D.new <- Reduce('+', D.update)/num)
  
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

EM <- function(long.formula, disp.formula, surv.formula, data, summax = 100, post.process = T, control = list(),
               delta.init=NULL){
  start.time <- proc.time()[3]
  
  #' Parsing formula objects ----
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  num <- surv$n
  
  #' Initial conditions ----
  inits.long <- Longit.inits(long.formula, disp.formula, data)
  inits.surv <- TimeVarCox(data, inits.long$b, surv$ph, formulas)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:num, function(i) inits.long$b[i, ])
  delta <- inits.long$delta.init
  if(!is.null(delta.init)) delta <- c(delta.init)
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
  if(!is.null(control$auto.summax)) auto.summax <- control$auto.summax else auto.summax <- F
  if(!is.null(control$net)) net <- control$net else net <- 0.5
  
  if(auto.summax){
    summax.old <- summax
    summax <- as.integer(max(sapply(Y, max)) * 2)
    if(verbose) cat(sprintf("Setting summax to %d\n", summax))
  }
  
  
  if(debug){DEBUG.inits <<- params; DEBUG.dmats <<- dmats; DEBUG.sv <<- sv; DEBUG.surv <<- surv}
  #' Gauss-hermite quadrature ----
  GH <- statmod::gauss.quad.prob(.gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  #' Matrices ---
  # Define mus and nus one expects at initial conditions...
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b)
  nus <- mapply(function(G) exp(G %*% delta), G = G)
  mm <- do.call(c, mus); nn <- do.call(c, nus); 
  nu.vec <- as.vector(unique(nn))
  # Generate mu vector
  .min <- max(0, min(mm) - net); .max <- max(mm) + net
  all.mus <- generate_mus(.min, .max)
  
  if(verbose) start.grids <- proc.time()[3]
  lambda.mat <- structure(gen_lambda_mat(.min, .max, nu.vec, summax),
                          dimnames = list(as.character(all.mus),
                                          as.character(nu.vec))
  )
  
  logZ.mat <- structure(gen_logZ_mat(.min, .max, nu.vec, lambda.mat, summax),
                        dimnames = list(as.character(all.mus),
                                        as.character(nu.vec))
  )
  
  V.mat <- structure(gen_V_mat(.min, .max, nu.vec, lambda.mat, logZ.mat, summax),
                     dimnames = list(as.character(all.mus),
                                     as.character(nu.vec))
  )
  if(verbose){
    end.grids <- proc.time()[3]
    d <- dim(lambda.mat)
    cat(sprintf("First %d x %d lambda, logZ and V grids calculated in %.3f seconds.\n", d[1], d[2], round(end.grids - start.grids, 3)))
  }
  
  #' Begin EM ----
  diff <- 100; iter <- 0
  EMstart <- proc.time()[3]
  step.times <- c()
  while(diff > tol && iter < maxit){
    update <- EMupdate(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                       sv, w, v, num, m, summax, debug, all.mus, nu.vec,
                       lambda.mat, logZ.mat, V.mat)
    
    if(debug) DEBUG.update <<- update
    params.new <- c(vech(update$D), update$beta, update$delta, update$gamma, update$zeta)
    names(params.new) <- names(params)
    
    # Convergence criteria + print (if wanted).
    diff <- difference(params, params.new, conv)
    if(verbose){
      print(sapply(params.new, round, 4))
      message("Iteration ", iter + 1, ' maximum ', conv, ' difference: ', round(diff, 4))
      message("Largest RE relative difference: ", round(max(difference(do.call(cbind, update$b), do.call(cbind, b), conv)), 4))
    }
    
    # create new matrices
    mm <- do.call(c, update$mus)
    .min <- max(0, min(mm) - net); .max <- max(mm) + net
    new.mus <- generate_mus(.min, .max)
    nus <- mapply(function(G) exp(G %*% update$delta), G = G)
    nn <- do.call(c, nus); nu.vec <- sort(c(unique(nn)))
    
    # Calculate \lambda, logZ and V at possible new nu values.
    if(verbose) start.grids <- proc.time()[3]
    lambda.new.nu <- structure(gen_lambda_mat(mu_lower = .min, mu_upper = .max,
                                              nus = nu.vec, summax = summax),
                               dimnames = list(as.character(new.mus),
                                               as.character(nu.vec)))
    
    logZ.new.nu <- structure(gen_logZ_mat(mu_lower = .min, mu_upper = .max,
                                          nus = nu.vec, lambdamat = lambda.new.nu, summax = summax),
                             dimnames = list(as.character(new.mus),
                                             as.character(nu.vec)))
    
    V.new.nu <- structure(gen_V_mat(mu_lower = .min, mu_upper = .max,
                                    nus = nu.vec, lambdamat = lambda.new.nu, logZmat = logZ.new.nu, summax = summax),
                          dimnames = list(as.character(new.mus),
                                          as.character(nu.vec)))
    
    # Fit these into the lambda grid
    new.order <- order(as.numeric(dimnames(lambda.new.nu)[[2]]))
    lambda.mat <- lambda.new.nu[,new.order,drop=F] # reorder columns to be in ascending \nu values.
    # logZ
    logZ.mat <- logZ.new.nu[,new.order,drop=F]
    # V
    V.mat <- V.new.nu[,new.order,drop=F]
    # Update all.mus and all.nus
    all.mus <- new.mus; nu.vec <- nu.vec
    rm(lambda.new.nu, logZ.new.nu, V.new.nu, mm, nn, new.mus) # large data objects
    
    if(verbose){
      end.grids <- proc.time()[3]
      d <- dim(lambda.mat)
      cat(sprintf("Iteration %d: %d x %d lambda, logZ and V grids obtained in %.3f seconds.\n", iter + 1, d[1], d[2], round(end.grids - start.grids, 3)))
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
              iter = iter)
  modelInfo <- list(
    forms = formulas,
    summax = summax
  )
  if(auto.summax) modelInfo$summax.type <- 'automatic' else modelInfo$summax.type <- 'manual'
  out$modelInfo <- modelInfo
  out$hazard <- cbind(ft = sv$ft, haz = l0, nev = sv$nev)
  out$stepmat <- cbind(iter = 1:iter, time = step.times)
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time.p <- proc.time()[3]
    #' Calculating \b and \Sigma at MLEs
    #' Find b.hat and Sigma
    b.hat <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
      optim(b, joint_density, joint_density_ddb,
            X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
            S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
            gamma = gamma, zeta = zeta, 
            lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, 
            all_mus = all.mus, all_nus = round(nu.vec, 3),
            summax = summax, method = 'BFGS')$par
    }, b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
    l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
    
    
    Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
      solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
                              S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                              gamma = gamma, zeta = zeta, lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, 
                              all_mus = all.mus, all_nus = round(nu.vec, 3), summax = summax, eps = .001))
    }, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
    l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
    # The Information matrix
    I <- structure(vcov(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, num, summax,
                        all.mus, nu.vec, lambda.mat, logZ.mat, V.mat),
                   dimnames = list(names(params), names(params)))
    I.inv <- tryCatch(solve(I), error = function(e) e)
    if(inherits(I.inv, 'error')) I.inv <- structure(MASS::ginv(I),
                                                    dimnames = dimnames(I))
    out$SE <- sqrt(diag(I.inv))
    out$vcov <- I
    out$RE <- do.call(rbind, b)
    out$postprocess.time <- round(proc.time()[3]-start.time.p, 2)
    
    #' Calculate the log-likelihood.
    # Timing done separately as EMtime + postprocess.time is for EM + SEs.
    #out$logLik <- log.lik(Omega, dmats, b, surv, sv, l0u, l0i, summax)
    #out$lltime <- round(proc.time()[3] - start.time, 2) - out$postprocess.time
  }
  out$comp.time = round(proc.time()[3]-start.time, 3)
  out
}
