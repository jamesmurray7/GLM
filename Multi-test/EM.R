library(survival)
library(glmmTMB)
library(Rcpp)
library(RcppArmadillo)
source('_Functions.R')
source('simData.R')
source('inits.R')
source('vcov2.R')
sourceCpp('funs.cpp')
vech <- function(x) x[lower.tri(x, T)]

EMupdate <- function(Omega, family, X, Y, Z, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, sv, w, v, n, m, optimiser, hessian,
                     beta.inds, b.inds, K, q){
  
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
  beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1) # Indexed for C++ use.
  gamma.rep <- rep(gamma, sapply(b.inds, length))

  #' Find b.hat and Sigma
  if(hessian == 'auto') .hess <- T else .hess <- F
  if(optimiser == 'optim'){
    #' Posterior mode and hessian via BFGS search (optim) ----
    b.update <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
      optim(b, joint_density, joint_density_ddb,
            Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
            Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
            beta_inds = beta.inds2, b_inds = b.inds2, K = K,
            method = 'BFGS', hessian = .hess)
    }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
    Fu = Fu, l0u = l0u, SIMPLIFY = F)
    b.hat <- lapply(b.update, function(x) x$par)
    if(hessian == 'manual'){
      Sigma <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
        solve(
          joint_density_sdb(b = b, Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family,
                            Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep,
                            zeta = zeta, beta_inds = beta.inds2, b_inds = b.inds2, K = K, eps = 1e-5))
      }, b = b.hat, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
      Fu = Fu, l0u = l0u, SIMPLIFY = F)
    }else{
      Sigma <- lapply(b.update, function(x) solve(x$hessian))
    }
    }else{
    #' Posterior mode via ucminf ----
    b.hat <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
      ucminf::ucminf(b, joint_density, joint_density_ddb,
                     Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
                     Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
                     beta_inds = beta.inds2, b_inds = b.inds2, K = K)$par
    }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
    Fu = Fu, l0u = l0u, SIMPLIFY = F)
    #' And its variance ----
    Sigma <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
      solve(
        joint_density_sdb(b = b, Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family,
                          Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep,
                          zeta = zeta, beta_inds = beta.inds2, b_inds = b.inds2, K = K, eps = 1e-5))
    }, b = b.hat, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
    Fu = Fu, l0u = l0u, SIMPLIFY = F)
    }
  SigmaSplit <- lapply(Sigma, function(x) lapply(b.inds, function(y) as.matrix(x[y,y])))
  bsplit <- lapply(b.hat, function(x) lapply(b.inds, function(y) x[y])) # Needed for updates to beta.
  bmat <- lapply(bsplit, bind.bs) # Needed for E[\ell(\gamma,\zeta)|\b...|\Omega].
  
  #' #########
  #' E-step ##
  #' #########
  
  #' D ----------------------------------------
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  #' \beta ------------------------------------
  # Sb <- mapply(function(X, Y, Z, b){
  #   Sbeta(beta, X, Y, Z, b, sigma, family, beta.inds2, K)
  # }, X = X, Y = Y, Z = Z, b = bsplit, SIMPLIFY = F)
  # Hb <- mapply(function(X, Y, Z, b){
  #   Hbeta(beta, X, Y, Z, b, sigma, family, beta.inds2, K, .Machine$double.eps^(1/4))
  # }, X = X, Y = Y, Z = Z, b = bsplit, SIMPLIFY = F)
  Sb <- Sbeta2(beta, X, Y, Z, bsplit, sigma, family, beta.inds2, K)
  Hb <- Hbeta2(beta, X, Y, Z, bsplit, sigma, family, beta.inds2, K, .Machine$double.eps^(1/3))
  # Sbq <- mapply(function(X, Y, Z, b, S){
  #   Sbeta_q(beta, X, Y, Z, b, sigma, family, beta.inds2, K, w, v, S)
  # }, X = X, Y = Y, Z = Z, b = bsplit, S = SigmaSplit, SIMPLIFY = F)
  # 
  # Hbq <- mapply(function(X, Y, Z, b, S){
  #   Hbeta_q(beta, X, Y, Z, b, sigma, family, beta.inds2, K, w, v, S, .Machine$double.eps^(1/3))
  # }, X = X, Y = Y, Z = Z, b = bsplit, S = SigmaSplit, SIMPLIFY = F)
  # Hb <- Reduce('+', Hbq); Sb <- Reduce('+', Sbq)

  
  #' Dispersion ('\sigma') --------------------
  if(any(unlist(family) %in% c('gaussian', 'negative.binomial'))){
    gauss.inds <- which(unlist(family) == 'gaussian')
    nb.inds <- which(unlist(family) == 'negative.binomial')
    sigma.update <- vector('list', K)
    # Gaussian -> Residual variance update
    for(j in gauss.inds){
      tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z[[j]] %*% S[[j]], Z[[j]])))), Z = Z, S = SigmaSplit)
      sigma.update[[j]] <- sum(mapply(function(X, Y, Z, b, tau){
        vare_update(X[[j]], Y[[j]], Z[[j]], b[[j]], beta[beta.inds[[j]]], tau, w, v)
      }, X = X, Y = Y, Z = Z, b = bsplit, tau = tau))
    }
    # Negative Binomial -> Scalar dispersion update.
    for(j in nb.inds){
      sigma.upate[[j]] <- mapply(function(X, Y, Z, b, Sigma){
        out <- vector('list', 2)
        out[[1]] <- Stheta(sigma[[j]], beta[beta.inds[[j]]], X[[j]], Y[[j]], Z[[j]], b[[j]], 
                           Sigma[[j]], w, v, .Machine$double.eps^(1/3))
        out[[2]] <- Htheta(sigma[[j]], beta[beta.inds[[j]]], X[[j]], Y[[j]], Z[[j]], b[[j]], 
                           Sigma[[j]], w, v, .Machine$double.eps^(1/4))
        out
      }, X = X, Y = Y, Z = Z, b = bsplit, Sigma = SigmaSplit, SIMPLIFY = F)
    }
    for(j in setdiff(1:K, c(gauss.inds, nb.inds))) sigma.update[[j]] <- NA # Return null for all those not gaussian or negative binomial.
  }else{
    sigma.update <- as.list(rep(NA, K))
  }
  
  #' Survival parameters (\gamma, \zeta) ------
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, q, .Machine$double.eps^(1/3))
  }, b = b.hat, Sigma = SigmaSplit, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)
  
  Hgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Hgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, q, .Machine$double.eps^(1/4))
  }, b = b.hat, Sigma = SigmaSplit, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  #' #########
  #' M-step ##
  #' #########
  
  # D
  D.new <- Reduce('+', D.update)/n
  
  # beta
  # beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
  beta.new <- beta-solve(Hb,Sb)
  
  # Dispersion
  ms <- colSums(do.call(rbind, m))
  if(any(unlist(family) %in% c('gaussian', 'negative.binomial'))){
    sigma.new <- vector('list', K)
    for(j in gauss.inds){ # Gaussian residual variance.
      sigma.new[[j]] <- sigma.update[[j]]/ms[j]
    }
    for(j in nb.inds){ # Negative binomial scalar dispersion.
      sigma.new[[j]] <- sigma[[j]] - sum(unlist(lapply(sigma.update[[j]], el, 1)))/sum(unlist(sigma.update[[j]], el, 2))
    }
    for(j in setdiff(1:K, c(gauss.inds, nb.inds))) sigma.new[[j]] <- 0.0 # Return null for all those not gaussian or negative binomial.
  }else{
    sigma.new <- as.list(rep(0.0, K))
  }
  
  # Survival parameters (gamma, zeta)
  gammazeta.new <- c(gamma, zeta) - solve(Reduce('+', Hgz), rowSums(Sgz))
  
  # The baseline hazard and related objects
  lambda.update <- lambdaUpdate(sv$surv.times, do.call(cbind, sv$ft.mat), gamma, gamma.rep, zeta, S, SigmaSplit, b.hat, w, v, b.inds2, K, q)
  l0.new <- sv$nev/rowSums(lambda.update)
  l0u.new <- lapply(l0u, function(ll){
    l0.new[1:length(ll)]
  })
  l0i.new <- l0.new[match(sv$Ti, sv$ft)] 
  l0i.new[is.na(l0i.new)] <- 0
  
  # Return
  list(
    D = D.new, beta = beta.new, sigma = sigma.new,       # Yk responses
    gamma = gammazeta.new[1:K], zeta = gammazeta.new[(K+1):length(gammazeta.new)],  # Survival
    l0 = l0.new, l0u = l0u.new, l0i = as.list(l0i.new),  #   Hazard
    b = b.hat                                            #   REs.
  )
  
}


EM <- function(long.formulas, surv.formula, data, family, post.process = T, control = list()){
  
  start.time <- proc.time()[3]
  
  #' Initial parsing ----
  formulas <- lapply(long.formulas, parseFormula)
  surv <- parseCoxph(surv.formula, data)
  n <- nrow(surv$survdata); K <- length(family)
  if(K!=length(long.formulas)) stop('Mismatched lengths of "family" and "long.formulas".')
  
  #' Initial conditons ----
  if(!is.null(control$dispformula)) dispformula <- control$dispformula else dispformula <- NULL
  inits.long <- Longit.inits(long.formulas, data, family, dispformula = dispformula)
  # Suss out indices of b_k and beta_k.
  b.inds <- lapply(1:K, function(k){
    nm <- inits.long$responses[k]
    which(grepl(nm, colnames(inits.long$b)))
  })
  beta.inds <- lapply(1:K, function(k){
    nm <- inits.long$responses[k]
    which(grepl(nm, names(inits.long$beta.init)))
  })
  q <- length(do.call(c, b.inds))
  
  inits.surv <- TimeVarCox(data, inits.long$b, surv$ph, formulas, b.inds)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  sigma <- inits.long$sigma.init # dispersion / resid. variance / 0 otherwise.
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  # Survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma\\_', names(inits.surv$inits))]
  
  #' Longitudinal and survival data objects ----
  sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
  dmats <- createDataMatrices(data, formulas)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  m <- lapply(Y, function(y) sapply(y, length))
  # survival
  Fi <- sv$Fi; Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u; Delta <- surv$Delta 
  l0 <- sv$l0
  S <- sv$S; SS <- sv$SS
  
  #' Assign family to joint density and parameter updates ----
  family <- lapply(family, function(family){
    if("function"%in%class(family)) family <- family()$family
    if(family == 'nbinom2') family <- 'negative.binomial' ## Catch-all.
    family
  })
  
  # Do we want to add-in correlated random effects between responses? For large K this greatly increases
  #   computation time and model instability.
  if(!is.null(control$correlated)) correlated <- control$correlated else correlated <- T
  if(!correlated) D[inits.long$off.inds] <- 0
  
  #' Parameter vector and list ----
  Omega <- list(D=D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, gamma, zeta)
  if(any(unlist(family)=='gaussian')){
    sigma.unlist <- unlist(sigma)
    # Variance - Gaussian.
    vare.name <- paste0(inits.long$responses[which(unlist(family) == 'gaussian')], '_var.e')
    vare <- sigma.unlist[which(unlist(family) == 'gaussian')]; names(vare) <- vare.name
    params <- c(params, vare)
  }
  if(any(unlist(family) == 'negative.binomial')){
    sigma.unlist <- unlist(sigma)
    # Dispersion - Negative binomial.
    theta.name <- paste0(inits.long$responses[which(unlist(family) == 'negative.binomial')], '_theta')
    theta <- sigma.unlist[which(unlist(family) == 'negative.binomial')]; names(theta) <- theta.name
    params <- c(params, theta)
  }
  
  #' Gauss-Hermite Quadrature ----
  if(!is.null(control$gh.nodes)) gh <- control$gh.nodes else gh <- 3
  if(!is.null(control$gh.sigma)) .sigma <- control$gh.sigma else .sigma <- 1
  
  GH <- statmod::gauss.quad.prob(gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  #' Begin EM ----
  diff <- 100; iter <- 0;
  if(!is.null(control$tol)) tol <- control$tol else tol <- 1e-2
  if(!is.null(control$maxit)) maxit <- control$maxit else maxit <- 200
  if(!is.null(control$conv)) conv <- control$conv else conv <- "relative"
  if(!conv%in%c('absolute', 'relative')) stop('Only "absolute" and "relative" difference convergence criteria implemented.')
  if(!is.null(control$verbose)) verbose <- control$verbose else verbose <- F
  if(!is.null(control$optimiser)) optimiser <- control$optimiser else optimiser <- 'optim'
  if(!optimiser %in% c('ucminf', 'optim')) stop("Only optimisers 'optim' and 'ucminf' supported.")
  if(!is.null(control$hessian)) hessian <- control$hessian else hessian <- 'manual'
  if(!hessian %in% c('auto', 'manual')) stop("Argument 'hessian' needs to be either 'auto' (i.e. from optim) or 'manual' (i.e. from _sdb, the defualt).")
  if(!is.null(control$SEs)) SEs <- control$SEs else SEs <- 'appx'
  if(!SEs %in% c('appx', 'exact', 'score')) stop("Argument gamma.SE needs to be either:\n",
                                                       "'appx' (the default) --> Calculated by observed empirical information matrix calculation\n",
                                                       "'exact' --> Using direct hessian calculation for survival terms\n",
                                                       "'score' --> Calculte the score and then hessian by forward differencing.")
  if(!is.null(control$SE.D)) SE.D <- control$SE.D else SE.D <- T
  if(!is.logical(SE.D)) stop("'SE.D' must be either TRUE or FALSE.")
  
  EMstart <- proc.time()[3]
  while(diff > tol && iter < maxit){
    update <- EMupdate(Omega, family, X, Y, Z, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, sv, w, v, n, m, optimiser, hessian, beta.inds, b.inds, K, q)
    if(!correlated) update$D[inits.long$off.inds] <- 0
    params.new <- c(vech(update$D), update$beta, update$gamma, update$zeta)
    if(any(unlist(family) %in%c('gaussian', 'negative.binomial'))){
      gauss.inds <- which(unlist(family) == 'gaussian')
      nb.inds <- which(unlist(family) == 'negative.binomial')
      params.new <- c(params.new, unlist(update$sigma)[gauss.inds], unlist(update$sigma)[nb.inds])
    }
    names(params.new) <- names(params)
    # Check convergence
    diff <- difference(params, params.new, conv)
    if(verbose){
      print(sapply(params.new, round, 4))
      message("Iteration ", iter + 1, ' maximum ', conv, ' difference: ', round(diff, 4))
    }

    #' Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; sigma <- update$sigma
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    iter <- iter + 1
    Omega <- list(D=D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
    params <- params.new
    
  }
  EMend <- proc.time()[3]
  coeffs <- Omega
  coeffs$beta <- setNames(c(Omega$beta), names(inits.long$beta.init))
  out <- list(coeffs = coeffs,
              RE = do.call(rbind, b),
              iter = iter,
              EMtime = EMend-EMstart,
              long.formulas = long.formulas, 
              surv.formula = surv.formula, 
              totaltime = proc.time()[3] - start.time)
  out$hazard <- cbind(ft = sv$ft, haz = l0, nev = sv$nev)
  out$family <- family
  
  out$ResponseInfo <- sapply(1:K, function(k){
    paste0(inits.long$responses[k], ' (', family[k], ')')
  })
  
  #' Post processing ----
  if(post.process){
    message('\nCalculating SEs...')
    beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1) 
    gamma.rep <- rep(gamma, sapply(b.inds, length))
    start.time <- proc.time()[3]
    
    #' b and Sigmai at MLEs
    b.update <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
      optim(b, joint_density, joint_density_ddb,
            Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
            Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
            beta_inds = beta.inds2, b_inds = b.inds2, K = K,
            method = 'BFGS', hessian = T)
    }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
    Fu = Fu, l0u = l0u, SIMPLIFY = F)
    Sigma <- lapply(b.update, function(x) solve(x$hessian))
    b <- lapply(b.update, function(x) x$par)
    SigmaSplit <- lapply(Sigma, function(x) lapply(b.inds, function(y) as.matrix(x[y,y])))
    bsplit <- lapply(b, function(x) lapply(b.inds, function(y) x[y])) # Needed for updates to beta.
    # The Information matrix
    if(SEs != 'score'){
      I <- structure(vcov(coeffs, dmats, surv, sv, 
                          Sigma, SigmaSplit, b, bsplit, 
                          l0u, w, v, n, family, K, q, beta.inds, b.inds, SEs),
                     dimnames = list(names(params), names(params)))
    }else{
      sigmas <- lapply(1:K, function(k){
        if(!family[[k]]%in%c('gaussian', 'negative.binomial')) NULL else Omega$sigma[[k]]
      })
      Omega2 <- list(D = vech(Omega$D), beta = c(Omega$beta), gamma = Omega$gamma, zeta = Omega$zeta, sigma = sigmas)
      np <- names(params)
      if(!SE.D){
        Omega2$D <- NULL
        np <- np[!grepl('^D\\[', np)] 
      }
      
      H <- structure(.H(unlist(Omega2), dmats, surv, sv, Sigma, SigmaSplit, b, bsplit, l0u, w, v, n, family, K, q, 
                        beta.inds, b.inds, verbose = verbose, SE.D = SE.D),
                        dimnames = list(np, np))
      I <- -H
    }
    I.inv <- tryCatch(solve(I), error = function(e) e)
    if(inherits(I.inv, 'error')) I.inv <- structure(MASS::ginv(I),
                                                    dimnames = dimnames(I))
    out$SE <- sqrt(diag(I.inv))
    out$vcov <- I
    
    out$RE <- do.call(rbind, b)
    out$postprocess.time <- round(proc.time()[3]-start.time, 2)
    # Calculate log-likelihood. Done separately as EMtime + postprocess.time is for EM + SEs.
    out$logLik <- log.lik(coeffs, dmats, b, surv, sv, l0u, l0i, gamma.rep, beta.inds2, b.inds2, K, q, family)
    out$lltime <- round(proc.time()[3] - start.time, 2) - out$postprocess.time
  }
  out
}
