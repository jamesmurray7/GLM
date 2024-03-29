#' #######
#' EM.R                 ((HYBRID2))
#' ----
#' At each iteration, constructs grid, or adds to existing grid...
#' ######

# rm(list=ls())
library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
sourceCpp('hybrid.cpp')
source('_Functions.R')
source('simData.R')
source('inits.R')
source('vcov.R')
vech <- function(x) x[lower.tri(x, T)]

EMupdate <- function(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                     sv, w, v, num, m, summax, debug, all.mus, nu.vec, lambda.mat, logZ.mat, V.mat, optimiser){
  s <- proc.time()[3]
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; delta <- Omega$delta; gamma <- Omega$gamma; zeta <- Omega$zeta
  
  #' Find b.hat and Sigma
  b.update <- b.minimise(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta, Omega, 
                         lambda.mat, logZ.mat, V.mat, 
                         all.mus, nu.vec, summax, optimiser, 'joint_density')
  b.hat <- b.update$b.hat
  Sigma <- b.update$Sigma
  
  if(debug){DEBUG.b.hat <<- b.hat; DEBUG.Sigma <<- Sigma}
  
  check <- sapply(Sigma, det)
  if(any(check < 0 | is.nan(check))){
    message('Some non pos-def or NaN Sigma...')
    ind <- which(check < 0 | is.nan(check))
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
  m <- mapply(function(a) get_indices(a, all.mus, 3), a = mus, SIMPLIFY = F)
  n <- mapply(function(a) get_indices(a, round(nu.vec, 2) , 2), a = nus, SIMPLIFY = F)
  
  #' \lambdas, Vs for \beta update. (With hardcode for NA/out-of-range values).
  lambdas <- mapply(function(mu, nu, m, n){
    out <- numeric(length(mu))
    if(any(is.na(m))){
      nas <- is.na(m)
      out[nas] <- lambda_appx(mu[nas], nu[nas], summax)
      out[!nas] <- mat_lookup(m[!nas], n[!nas], lambda.mat)
    }else{
      out <- mat_lookup(m, n, lambda.mat)
    }
    out
  }, mu = mus, nu = nus, m = m, n = n, SIMPLIFY = F)
  
  Vs <- mapply(function(mu, nu, m, n){
    out <- numeric(length(mu))
    if(any(is.na(m))){
      nas <- is.na(m)
      lams <- lambda_appx(mu[nas], nu[nas], summax)
      logZs <- logZ_c(log(lams), nu[nas], summax)
      out[nas] <- calc_V_vec(mu[nas], lams, nu[nas], logZs, summax)
      out[!nas] <- mat_lookup(m[!nas], n[!nas], V.mat)
    }else{
      out <- mat_lookup(m, n, V.mat)
    }
    out
  }, mu = mus, nu = nus, m = m, n = n, SIMPLIFY = F)
  
  #' D
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  #' \beta
  Sb <- mapply(Sbeta, X, Y, mus, nus, lambdas, Vs, SIMPLIFY = F)
  Hb <- mapply(getW1, X, mus, nus, lambdas, Vs, SIMPLIFY = F)
  
  #' \delta
  Sd <- mapply(function(G, b, X, Z, Y, lY, tau){
    Sdelta_cdiff(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps = .Machine$double.eps^(1/3))
  }, G = G, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)
  
  Hd <- mapply(function(G, b, X, Z, Y, lY, tau){
    Hdelta(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps = .Machine$double.eps^(1/4))
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

EM <- function(long.formula, disp.formula, surv.formula, data, summax = 100, post.process = T, 
               control = list(), disp.control = list(),  delta.init = NULL, grid.summax = 'same', return.matrices = F,
               optimiser = 'bobyqa'){
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
  if(!is.null(control$mu.rule)) mu.rule <- control$mu.rule else mu.rule <- 'quantile'
  if(!is.null(control$mu.rule.probs)) mu.rule.probs <- control$mu.rule.probs else mu.rule.probs <- c(.25, .75)
  if(!mu.rule%in%c('range', 'quantile')) stop("'mu.rule' must be either 'quantile' or 'range'.")
  if(any(mu.rule.probs > 1) | length(mu.rule.probs) != 2) stop("Provide 'mu.rule.probs' as c(lower, upper).")
  #' Control arguments specific to dispersion estimates ----
  if(!is.null(disp.control$delta.method)) delta.method <- disp.control$delta.method else delta.method <- 'optim'
  if(!delta.method %in% c('uniroot', 'optim')) stop('delta.method must be either "optim" or "uniroot".\n')
  if(!is.null(disp.control$min.profile.length)) min.profile.length <- disp.control$min.profile.length else min.profile.length <- 1
  if(!is.null(disp.control$what)) what <- disp.control$what else what <- 'median'
  if(!what %in% c('mean', 'median')) stop("what must be either 'mean' or 'median'.")
  if(!is.null(disp.control$cut)) cut <- disp.control$cut else cut <- T
  if(!is.null(disp.control$re.maximise)) re.maximise <- disp.control$re.maximise else re.maximise <- T
  if(!is.null(disp.control$percentiles)) percentiles <- disp.control$percentiles else percentiles <- c(.4, .6)
  if(any(percentiles > 1) | length(percentiles) != 2) stop("Provide percentiles as c(lower, upper).")
  if(!is.null(disp.control$interval)) interval <- disp.control$interval else interval <- c(-2, 2)
  
  if(auto.summax){
    summax.old <- summax
    summax <- max(sapply(Y, max)) * 2
    if(verbose) cat(sprintf("Automatically setting summax to %d\n", summax))
  }
  
  # Initial conditon for delta
  if(is.null(delta.init)){
    delta.inits.raw <- get.delta.inits(dmats, beta, b, delta.method, summax, verbose, min.profile.length, percentiles, interval)  
    # Return the user-specified POINT estimate (Defaulting to cut + median)
    if(cut){
      initdelta <- if(what == 'mean') delta.inits.raw$mean.cut.estimate else delta.inits.raw$median.cut.estimate
      pc.delta <- delta.inits.raw$pc.cut
    }else{
      initdelta <- if(what == 'mean') delta.inits.raw$mean.estimate else delta.inits.raw$median.estimate
      pc.delta <- delta.inits.raw$pc.raw
    }
    delta <- setNames(initdelta, names(inits.long$delta.init))
    if(verbose) cat(sprintf('Initial condition for delta: %.3f. [%.3f, %.3f].\n', 
                            delta, pc.delta[1], pc.delta[2]))
  }else{
    delta <- setNames(delta.init, names(inits.long$delta.init))
  }
  
  #' Parameter vector and list ----
  Omega <- list(D=D, beta = beta, delta = delta, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, delta, gamma, zeta)
  
  
  if(debug){DEBUG.inits <<- params; DEBUG.dmats <<- dmats; DEBUG.sv <<- sv; DEBUG.surv <<- surv}
  #' Gauss-hermite quadrature ----
  GH <- statmod::gauss.quad.prob(.gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  #' Re-maximise Marginal wrt b for CMP rather than Poisson (obtained thus far)
  if(re.maximise){
    s <- proc.time()[3]
    if(verbose) cat('Re-maximising in light of new delta estimate')
    b <- b.minimise(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta, Omega, 
                    lambda.mat = NULL, logZ.mat = NULL, V.mat = NULL, 
                    all.mus = NULL, nu.vec = NULL,
                    summax, optimiser, 'marginal')$b.hat
    remaximisation.time <- proc.time()[3]-s
    if(verbose) cat(sprintf(', this took %.2f seconds.\n', remaximisation.time))
  }
  startup.time <- round(proc.time()[3] - start.time, 3)
    
  #' Matrices ----------------------------
  #' Define mus ----
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b)
  mm <- do.call(c, mus)
  # Generate mu vector
  if(mu.rule == 'range'){
    .min <- max(0, min(mm)); .max <- max(mm)
  }else{
    qn <- unname(quantile(mm, prob = mu.rule.probs))
    .min <- max(0, qn[1]); .max <- qn[2]
  }
  
  all.mus <- generate_mus(.min, .max)
  
  #' Define nus ----
  allG <- apply(do.call(rbind, G), 2, unique)
  if(!is.list(allG)) allG <- as.list(allG)
  allG <- do.call(cbind, allG)
  nn.lb <- min(exp(allG %*% pc.delta[1])); nn.ub <- max(exp(allG %*% pc.delta[2]))
  nu.vec <- seq(round(nn.lb, 2), round(nn.ub, 2), 1e-2)

  if(verbose)
    cat(sprintf("%d-vector of nu values generated from %.2f -- %.2f.\n", 
                length(nu.vec), min(nu.vec), max(nu.vec)))
  
  grid.times <- as.matrix(t(setNames(numeric(2), c('iteration', 'elapsed time'))))
  start.grids <- proc.time()[3]
  if(grid.summax == 'same') grid.summax <- summax else grid.summax <- max(20, floor(summax / 4))
  lambda.mat <- structure(gen_lambda_mat(.min, .max, nu.vec, grid.summax),
                          dimnames = list(as.character(all.mus),
                                          as.character(nu.vec)))
  
  logZ.mat <- structure(gen_logZ_mat(.min, .max, nu.vec, lambda.mat, grid.summax),
                         dimnames = list(as.character(all.mus),
                                         as.character(nu.vec))
  )
  
  V.mat <- structure(gen_V_mat(.min, .max, nu.vec, lambda.mat, logZ.mat, grid.summax),
                     dimnames = list(as.character(all.mus),
                                     as.character(nu.vec))
  )
  end.grids <- proc.time()[3]
  grid.times[1,] <- c(0, end.grids - start.grids)
  if(verbose){
    d <- dim(lambda.mat)
    cat(sprintf("First %d x %d lambda, logZ and V grids calculated in %.3f seconds on summax %d.\n\n", d[1], d[2], end.grids - start.grids, grid.summax))
  }
  
  #' Begin EM ----
  diff <- 100; iter <- 0
  EMstart <- proc.time()[3]
  step.times <- c()
  while(diff > tol && iter < maxit){
    update <- EMupdate(Omega, X, Y, lY, Z, G, b, S, SS, Fi, Fu, l0i, l0u, Delta, l0, 
                       sv, w, v, num, m, summax, debug, all.mus, nu.vec,
                       lambda.mat, logZ.mat, V.mat, optimiser)
    
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
    
    #' Create new matrices ------------------------------------------------
    mm.new <- do.call(c, update$mus)
    
    #' Define new mus ----
    if(mu.rule == 'range'){
      .min.new <- max(0, min(mm.new)); .max.new <- max(mm.new)
    }else{
      qn <- quantile(mm.new, prob = mu.rule.probs)
      .min.new <- max(0, qn[1]); .max.new <- qn[2]
    }
    
    #' Define new nus ----
    update.nu <- exp(allG %*% update$delta)
    
    #' Work out A/B1/B2/C matrices ----
    start.grids <- proc.time()[3]
    #' ======
    #' A ====
    #' ======
    
    #' Elements calculated at old mu values at new values of nu < nu.old.
    if(round(update.nu, 2) < min(nu.vec)){
      # If lower than smallest value, sequence from the new nu to the minimal value
      temp.nu <- unique(seq(round(update.nu, 2), min(nu.vec) - 1e-2, 1e-2)) 
      Amats <- gen_all_mats(.min, .max, temp.nu, summax)
      nu.vec <- c(temp.nu, nu.vec)
    }else{
      # Otherwise, create an empty matrix to append to left of existing matrix.
      Amats <- setNames(replicate(3, lambda.mat[,0,drop=F], simplify = F),
                    c('lambda', 'logZ', 'V'))
    }
    num.new <- unique(do.call(rbind, lapply(Amats, dim))[,2])
    if(verbose) cat(sprintf('%d new nu values < old nu.\n', num.new))
    
    #' ======
    #' C ====
    #' ======
    
    #' Elements calculated at old mu values at new values of nu > nu.old
    if(round(update.nu, 2) > max(nu.vec)){
      # If greater than largest value, sequence from maximal value to the new nu.
      temp.nu <- unique(seq(max(nu.vec) + 1e-2, round(update.nu, 2), 1e-2)) 
      Cmats <- gen_all_mats(.min, .max, temp.nu, summax)
      nu.vec <- c(nu.vec, temp.nu)
    }else{
      Cmats <- setNames(replicate(3, lambda.mat[,0,drop=F], simplify = F),
                        c('lambda', 'logZ', 'V'))
    }
    num.new <- unique(do.call(rbind, lapply(Cmats, dim))[,2])
    if(verbose) cat(sprintf('%d new nu values > old nu.\n', num.new))
    
    #' ======
    #' B ====
    #' ======

    #' Elements calculated at new mu values at old values of nu.
    if(round(.min.new, 3) < round(.min, 3)){
      # If no difference between .min.new and (.min-0.001), then only need to generate matrices for .min.new.
      if(no.diff(round(.min.new, 3), round(.min, 3) - 1e-3)){
        B1mats <- gen_all_mats(round(.min.new, 3), round(.min.new, 3), nu.vec, summax)
        all.mus <- c(round(.min.new, 3), all.mus)
      }else{
        B1mats <- gen_all_mats(round(.min.new, 3), min(all.mus) -1e-3, nu.vec, summax)
        all.mus <- c(seq(round(.min.new, 3), min(all.mus) - 1e-3, 1e-3), all.mus)
      }
      .min <- .min.new
    }else{
      B1mats <- setNames(replicate(3, matrix(0, nr = 0, nc = length(nu.vec)), simplify = F),
                         c('lambda', 'logZ', 'V'))
    }
    num.new <- unique(do.call(rbind, lapply(B1mats, dim))[,1])
    if(verbose) cat(sprintf('%d new mu values < old mu.\n', num.new))
      
    if(round(.max.new, 3) > round(.max, 3)){ 
      if(no.diff(round(.max.new, 3), round(.max, 3) + 1e-3)){ # Coding for special case where new element is one more than existing max.
        B2mats <- gen_all_mats(round(.max, 3), round(.max.new, 3), nu.vec, summax)
        all.mus <- c(all.mus, round(.max.new, 3))
      }else{
        B2mats <- gen_all_mats(max(all.mus) + 1e-3, round(.max.new, 3), nu.vec, summax)
        all.mus <- c(all.mus, seq(max(all.mus) + 1e-3, round(.max.new, 3), 1e-3))
      }
      .max <- .max.new
    }else{
      B2mats <- setNames(replicate(3, matrix(0, nr = 0, nc = length(nu.vec)), simplify = F),
                         c('lambda', 'logZ', 'V'))
    }
    num.new <- unique(do.call(rbind, lapply(B2mats, dim))[,1])
    if(verbose) cat(sprintf('%d new mu values > old mu.\n', num.new))
    
    end.grids <- proc.time()[3]
    
    #' Create new matrices ------------------------------------------------
    lambda.mat.new <- rbind(B1mats$lambda,
                            cbind(cbind(Amats$lambda, lambda.mat), Cmats$lambda),
                            B2mats$lambda)
    logZ.mat.new  <- rbind(B1mats$logZ,
                           cbind(cbind(Amats$logZ, logZ.mat), Cmats$logZ),
                           B2mats$logZ)
    V.mat.new  <- rbind(B1mats$V,
                        cbind(cbind(Amats$V, V.mat), Cmats$V),
                        B2mats$V)
    
    if(verbose){
      d <- dim(lambda.mat.new)
      if(all(dim(lambda.mat.new) == dim(lambda.mat))) 
        cat("No change from previous iteration's matrices!\n\n")
      else
        cat(sprintf("Iteration %d: %d x %d lambda, logZ and V grids obtained in %.3f seconds.\n\n", iter + 1, d[1], d[2], round(end.grids - start.grids, 3)))
    }
    grid.times <- rbind(grid.times, c(iter + 1, end.grids - start.grids))
    
    rm(Amats, B1mats, B2mats, Cmats) # (potentially) large data objects.

    #' Update parameter estimates, matrices, etc...
    # Parameters
    b <- update$b
    D <- update$D; beta <- update$beta; delta <- update$delta
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    iter <- iter + 1
    Omega <- list(D = D, beta = beta, delta = delta, gamma = gamma, zeta = zeta)
    params <- params.new
    
    # Matrix stuff
    # {all.mus, nu.vec} is done in creation of {A,B,C} matrices.
    # .min <- .min.new; .max <- .max.new
    lambda.mat <- lambda.mat.new; logZ.mat <- logZ.mat.new; V.mat <- V.mat.new
    rm(lambda.mat.new, logZ.mat.new, V.mat.new)
    
    # This iteration's time.
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
  out$modelInfo <- modelInfo
  if(is.null(delta.init)) out$modelInfo$delta.inits <- delta.inits.raw
  out$hazard <- cbind(ft = sv$ft, haz = l0, nev = sv$nev)
  out$stepmat <- cbind(iter = 1:iter, time = step.times)
  out$gridtimes <- grid.times
  
  if(post.process){
    message('\nCalculating SEs...')
    start.time.p <- proc.time()[3]
    #' Calculating \b and \Sigma at MLEs
    #' Find b.hat and Sigma
    b.update <- b.minimise(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta, Omega, 
                           lambda.mat, logZ.mat, V.mat, 
                           all.mus, nu.vec, summax, optimiser, 'joint_density')
    b.hat <- b.update$b.hat
    Sigma <- b.update$Sigma
  
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
    postprocess.time <- round(proc.time()[3]-start.time.p, 2)
    
    #' Calculate the log-likelihood.
    # Timing done separately as EMtime + postprocess.time is for EM + SEs.
    #out$logLik <- log.lik(Omega, dmats, b, surv, sv, l0u, l0i, summax)
    #out$lltime <- round(proc.time()[3] - start.time, 2) - out$postprocess.time
  }
  comp.time <-  round(proc.time()[3]-start.time, 3)
  
  out$elapsed.time <- c(`delta optimisation` = if(is.null(delta.init)) unname(delta.inits.raw$time) else NULL,
                        `re-maximisation` = if(re.maximise) unname(remaximisation.time) else NULL,
                        `startup time` = unname(startup.time),
                        `EM time` = unname(EMtime),
                        `post processing` = if(post.process) unname(postprocess.time) else NULL,
                        `Total computation time` = unname(comp.time))
  
  if(return.matrices)
    out$final.matrices <- list(`lambda` = lambda.mat, `logZ` = logZ.mat, `V` = V.mat)
  
  out
}
