#' #######
#' EM.R
#' --
#' This largely based on my PBC-case-study/EM, inits, _Functions.R files
#' ######

library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
source('inits.R')
source('grid-vcov.R')
sourceCpp('grid-test.cpp')
vech <- function(x) x[lower.tri(x, T)]
# Load parameter matrices
if(!"N"%in%ls()) stop("N MUST BE DEFINED PRE-SOURCE!\n")
if(!"pete.flag"%in%ls()) stop("pete.flag MUST BE DEFINED PRE-SOURCE!\n")
if(!N%in%c(1e3,1e4)) stop("N MUST EITHER 1000 OR 10,000!\n")  
source('_Functions.R')
# Load parameter matrices
if(!'lambda.mat' %in% ls()) lambda.mat <- load.grid(N, 'lambda', pete.flag)
if(!'V.mat' %in% ls()) V.mat <- load.grid(N, 'V', pete.flag)
if(!'logZ.mat' %in% ls()) logZ.mat <- load.grid(N, 'logZ', pete.flag)
.mats <- c('lambda.mat', 'V.mat', 'logZ.mat')
.rm <- function() rm(list = setdiff(ls(), ls(pattern = '\\.mat')))
message('Parameter matrices loaded...')

EMupdate <- function(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                     sv, w, v, n, m, summax, debug, N){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; delta <- Omega$delta; gamma <- Omega$gamma; zeta <- Omega$zeta
  
  #' Find b.hat and Sigma
  b.hat <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
    optim(b, joint_density, joint_density_ddb,
          X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
          S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
          gamma = gamma, zeta = zeta, 
          lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, N = N, method = 'BFGS')$par
  }, b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
  l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)


  Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
    solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
                            S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                            gamma = gamma, zeta = zeta, lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, N = N, eps = 10/N))
  }, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
  l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  if(debug){DEBUG.b.hat <<- b.hat; DEBUG.Sigma <<- Sigma}
  
  check <- sapply(Sigma, det)
  if(any(check < 0)){
    message('Some non pos def Sigma...')
    ind <- which(check<0)
    for(i in seq_along(ind)){
      Sigma[[ind[i]]] <- as.matrix(Matrix::nearPD(Sigma[[ind[i]]], corr = T)$mat)
    }
  }
  #' #########
  #' E-step ##
  #' #########
  
  nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b.hat, SIMPLIFY = F)
  
  mus2 <- lapply(mus, mu_fix, N)
  nus2 <- lapply(nus, mu_fix, N)
  
  #' Grid lookups ----
  lambdas <- mapply(function(mu, nu){
    m <- (mu*(N/10)) - 1; n <- (nu*(N/10)) - 1
    mat_lookup(m, n, lambda.mat)
  }, mu = mus2, nu = nus2, SIMPLIFY = F)
  
  Vs <- mapply(function(mu, nu){
    m <- (mu*(N/10)) - 1; n <- (nu*(N/10)) - 1
    mat_lookup(m, n, V.mat)
  }, mu = mus2, nu = nus2, SIMPLIFY = F)
  
  logZ <- mapply(function(mu, nu){
    m <- (mu*(N/10)) - 1; n <- (nu*(N/10)) - 1
    mat_lookup(m, n, logZ.mat)
  }, mu = mus2, nu = nus2, SIMPLIFY = F)
  
  ABC <- mapply(function(b, X, Z, G, tau){
    A_B_C(b, X, Z, beta, delta, G, tau, 100, N, lambda.mat, logZ.mat)
  }, b = b.hat, X = X, Z = Z, G = G, tau = tau, SIMPLIFY = F)

  # D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  # \beta
  Sb <- mapply(Sbeta, X, Y, mus2, nus2, lambdas, Vs, SIMPLIFY = F)
  Hb <- mapply(getW1, X, mus2, nus2, lambdas, Vs, SIMPLIFY = F)
  
  # \delta
  Sd <- mapply(function(G, b, X, Z, Y, lY, tau){
    Sdelta_cdiff(delta, G, b, X, Z, Y, lY, beta, tau, w, v, N, lambda.mat, logZ.mat)
  }, G = G, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)
  Hd <- mapply(getW2, ABC, Vs, nus2, G, SIMPLIFY = F)
  
  # Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 10/N)
  }, b = b.hat, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)
  
  Hgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Hgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 10/N)
  }, b = b.hat, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  #' #########
  #' M-step ##
  #' #########
  
  # D
  D.new <- Reduce('+', D.update)/n
  
  # \beta and \delta
  beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
  delta.new <- delta - solve(Reduce('+', Hd), c(Reduce('+', Sd)))

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

EM <- function(long.formula, disp.formula, surv.formula, data, summax = 100, N, post.process = T, control = list(),
               delta.init=NULL){
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
  if(!is.null(control$summax.override)) summax.override <- control$summax.override else summax.override <- F
  
  
  if(debug){DEBUG.inits <<- params; DEBUG.dmats <<- dmats; DEBUG.sv <<- sv; DEBUG.surv <<- surv}
  #' Gauss-hermite quadrature ----
  GH <- statmod::gauss.quad.prob(.gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  #' Begin EM ----
  diff <- 100; iter <- 0
  EMstart <- proc.time()[3]
  step.times <- c()
  while(diff > tol && iter < maxit){
    update <- EMupdate(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                       sv, w, v, n, m, summax, debug, N)
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
            gamma = gamma, zeta = zeta, 
            lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, N = N, method = 'BFGS')$par
    }, b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
    l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
    
    Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
      solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
                              S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                              gamma = gamma, zeta = zeta, lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, N = N, eps = 10/N))
    }, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
    l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
    # The Information matrix
    I <- structure(vcov(Omega, dmats, surv, sv, Sigma, b, l0u, w, v, n, N, summax),
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
    #out$logLik <- log.lik(Omega, dmats, b, surv, sv, l0u, l0i, summax)
    #out$lltime <- round(proc.time()[3] - start.time, 2) - out$postprocess.time
  }
  
  out
}
