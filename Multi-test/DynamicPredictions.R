#' #####
#' DynamicPredictions.R
#' #####

# Prepare longitudinal data -----------------------------------------------
prepareLongData <- function(data, formulas, responses, u = NULL){
  #' Prepares data for one subject i. That is, data argument should be data[data$id==i, ] type
  #' Formulas/responses should come from fitted aEM object
  #' If u is NULL then the whole data[,] object is partitioned and returned.
  #'     One should use the `u` argument to project from Tstart in the Tstart+delta window.
  
  formulas <- lapply(formulas, parseFormula); K <- length(formulas)
  responses <- lapply(strsplit(responses, '\\s\\('), el, 1)
  
  if(is.null(u)) data <- data[data$time <= u, ] # Truncate if necessary
  
  X <- lapply(1:K, function(k){
    .DataWorkhorse(data, data$id[1], formulas[[k]]$fixed, formulas[[k]]$random, responses[[k]], what = 'X')
  })
  Y <- lapply(1:K, function(k){
    .DataWorkhorse(data, data$id[1], formulas[[k]]$fixed, formulas[[k]]$random, responses[[k]], what = 'Y')
  })
  Z <- lapply(1:K, function(k){
    .DataWorkhorse(data, data$id[1], formulas[[k]]$fixed, formulas[[k]]$random, responses[[k]], what = 'Z')
  })
  
  list(X = X, Y = Y, Z = Z)
}


# Prepare survival data ---------------------------------------------------
prepareSurvData <- function(data, surv.formula, formulas, u = NULL, hazard){
  #' data: i.dat-type object; survival data is then derived in-formula (by taking the first row).
  #' surv.formula: formula from the fitted aEM object.
  #' formulas: list of formulae from the fitted aEM object.
  #' u: horizon time; if not supplied then this is found for the entire profile is returned.
  #' hazard: hazard from fitted aEM object.
  
  formulas <- lapply(formulas, parseFormula); K <- length(formulas)
  
  ph <- coxph(surv.formula, data, x = T) # null fit for later.
  
  #' Extract failure times and \lambda_0 from hazard
  ft <- hazard[,1]; l0 <- hazard[,2]
  l0 <- predict(lm(l0 ~ splines::ns(ft, df = 3))) # Use smoothed cubic spline in hazard for predictions
  survtime <- data$survtime[1]
  status <- as.integer(data$status[1])
  data$longtime <- data$time; max.longtime <- max(data$longtime)
  data$time <- data$survtime # potentially bad move, we'll see!
  data <- data[1, ]
  
  if(!is.null(u) && data$time > u) data$time <- u # Don't extrapolate beyond time u
  
  #' Fi
  Fi <- do.call(cbind, lapply(formulas, function(x){
    model.matrix(as.formula(paste0('~', x$random)), data)
  }))
  
  #' Fu and l0u
  # Work out survived times (and if necessary truncate at u if supplied)
  survived.times <- ft[which(ft <= data$time)]
  
  if(!is.null(u)){
    survived.times <- survived.times[which(survived.times <= u)]
  }else{
    survived.times <- survived.times[which(survived.times <= max.longtime)]
  }
  print(survived.times)
  if(length(survived.times) > 0){ # Avoid complications if they're censored before first failure time.
    Fu <- do.call(cbind, lapply(formulas, function(x){
      model.matrix(as.formula(paste0('~', x$random)), data.frame(time=survived.times))
    }))
    l0u <- l0[match(survived.times, ft)]
  }else{
    Fu <- matrix(0, nr = 1, nc = ncol(Fi))
    l0u <- 0
  }
  
  if(status == 1L){
    l0i <- l0[which(ft == survtime)]
    if(!is.null(u) && u < survtime) l0i <- l0[max(which(ft <= u))] # truncate if necessary
  }else{
    l0i <- 0
  }
  
  S <- ph$x[1,,drop=F]
  SS <- apply(S, 2, rep, nrow(Fu))
  
  list(
    S = S, SS = SS,
    l0u = l0u, l0i = l0i,
    Fi = Fi, Fu = Fu, Delta = status
  )
}


# Collate and prepare all -------------------------------------------------
prepareData <- function(data, id, fit, b.inds, beta.inds, u = NULL){ # b.inds, beta.inds?
  long <- prepareLongData(data, fit$long.formulas, fit$ResponseInfo, u = u)
  surv <- prepareSurvData(data, fit$surv.formula, fit$long.formulas, u = u, hazard = fit$hazard)
  
  b <- fit$RE[id,]; K <- length(fit$family)
  responsenames <- lapply(strsplit(fit$ResponseInfo, '\\s\\('), el , 1)
  
  gamma.rep <- rep(fit$coeffs$gamma, sapply(b.inds, length))
  
  if(is.null(u)){
    Sigma <- solve(optim(
      b, joint_density, joint_density_ddb, Y = long$Y, X = long$X, Z = long$Z,
      beta = fit$coeffs$beta, D = fit$coeffs$D, sigma = fit$coeffs$D, family = fit$family,
      Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i, SS = surv$SS, Fu = surv$Fu,
      haz = surv$l0u, gamma_rep = gamma.rep, zeta = fit$coeffs$zeta, beta_inds = beta.inds,
      b_inds = b.inds, K = K, method = 'BFGS', hessian = T)$hessian
    )
  }else{
    Sigma <- NA
  }
  
  list(
    long = long, surv = surv, b = b, Sigma = Sigma
  )
  
}

# Draws from \Omega and \b ------------------------------------------------
Omega.draw <- function(fit){
  
  # Populat mean and covariance
  sigma.unlist <- unlist(fit$coeffs$sigma)
  Omega.mean <- setNames(
    c(vech(fit$coeffs$D), 
    fit$coeffs$beta,
    fit$coeffs$gamma,
    fit$coeffs$zeta,
    sigma.unlist[which(sigma.unlist != 0.0)]),
    names(fit$SE)
  )
  if(any(names(Omega.mean)!=dimnames(fit$vcov)[[1]])) stop('Something wrong -- Check Omega mean and vcov structure.')
  Omega.vcov <- fit$vcov
  
  draw <- MASS::mvrnorm(n = 1, mu = Omega.mean, Sigma = solve(Omega.vcov))
  
  D <- vech2mat(draw[grepl('^D\\[', names(draw))], ncol(fit$RE))
  beta <- draw[match(names(fit$coeffs$beta), names(draw))]
  gamma <- draw[match(names(fit$coeffs$gamma), names(draw))]
  zeta <- draw[match(names(fit$coeffs$zeta), names(draw))]
  if(any(sigma.unlist != 0)){
    .sigma <- draw[grepl('var\\.e$|theta$', names(draw))]
  }
  
  # Get sigma into same form as fit$sigma (i.e. zero if not Gaussian/NB)
  
  if(any(unlist(fit$family) %in% c('gaussian', 'negative.binomial'))){
    gauss.inds <- which(unlist(fit$family) == 'gaussian')
    nb.inds <- which(unlist(fit$family) == 'negative.binomial')
    sigma <- vector('list', length(which(unlist(fit$family) %in% c('gaussian', 'negative.binomial'))))
    for(j in gauss.inds){ # Gaussian residual variance.
      sigma[[j]] <- .sigma[j]
    }
    for(j in nb.inds){ # Negative binomial scalar dispersion.
      sigma[[j]] <- .sigma[[j]]
    }
    for(j in setdiff(1:length(fit$family), c(gauss.inds, nb.inds))) sigma[[j]] <- 0.0 # Return null for all those not gaussian or negative binomial.
  }else{
    sigma <- as.list(rep(0.0, K))
  }
  
  list(
    D = D, beta = beta, gamma = gamma, zeta = zeta, sigma = sigma
  )
  
}

b.draw <- function(b, long, surv, O, beta.inds, b.inds, fit){
  # Unpack Omega.draw object
  beta <- O$beta; D <- O$D; gamma <- rep(O$gamma, sapply(b.inds, length)); zeta <- O$zeta; sigma <- O$sigma
  
  # Use optim to draw from f(b,Y,T,Delta;Omega^{\ell})
  b.l <- optim(
    b, joint_density, joint_density_ddb,
    Y = long$Y, X = long$X, Z = long$Z, beta = beta, D = D, sigma = sigma,
    family = fit$family, Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i,
    SS = surv$SS, Fu = surv$Fu, haz = surv$l0u, gamma_rep = gamma, zeta = zeta, beta_inds = beta.inds,
    b_inds = b.inds, K = length(b.inds), method = 'BFGS', hessian = T
  )
  
  Sigma.l <- solve(b.l$hessian)
  b.hat.l <- b.l$par
  # check positive semi-definite
  xformed <- 0
  if(any(eigen(Sigma.l)$val<0)||det(Sigma.l)<=0){
    xformed <- 1
    Sigma.l <- as.matrix(Matrix::nearPD(Sigma.l))
  }
  
  list(
    b.hat = MASS::mvrnorm(n = 1, mu = b.hat.l, Sigma = Sigma.l),
    Sigma.l = Sigma.l,
    xformed = xformed # internal use - how often is covariance matrix poorly defined?
  )
  
}

# Dynamic predictions -----------------------------------------------------
# Dynamic survial predictions for ONE subject `id` across a vector of candidate times, u.
dynSurv <- function(data, id, fit, u = NULL, nsim = 200, progress = T){
  ft <- fit$hazard[,1];tmax <- max(ft); K <- length(fit$family)
  if(!is.null(u) & any(u > tmax)) stop("Can't extrapolate beyond last failure time.")
  
  # Subset the required subject
  newdata <- data[data$id == id, ] # subset required subject
  
  # If u isn't supplied then arbitrarily find probability of surviving until next failure time.
  #              (for each longitudinal time recorded)
  if(is.null(u)){
    u <- sapply(newdata$time, function(x) ft[which(ft > x)][1])
    if(any(u > tmax)) u <- u[!which(u > tmax)] # and ensure those after T_{max} aren't included
  }
  
  # Get indices for \b and \beta
  responsenames <- lapply(strsplit(fit$ResponseInfo, '\\s\\('), el , 1)
  nm <- colnames(fit$RE)
  b.inds <- lapply(1:K, function(k){ 
    which(grepl(responsenames[[k]], nm)) - 1
  })
  nm <- names(fit$coeff$beta)
  beta.inds <- lapply(1:K, function(k){
    which(grepl(responsenames[[k]], nm)) - 1
  })                                 
  
  # Obtain 'denominator' dataset
  newdata2 <- newdata[newdata$time <= u[1], ]
  data.t <- prepareData(newdata2, id = id, fit = fit, b.inds = b.inds, beta.inds = beta.inds, u = NULL)
  
  u <- u[-1] # Don't want to find preds for first time (T_{start})...
  
  pi <- setNames(vector('list', length = length(u)), paste0('u = ', u))
  for(uu in seq_along(u)){
    data.u <- prepareData(newdata2, id = id, b.inds, beta.inds, fit = fit,  u = u[uu])
    
    pi.store <- numeric(nsim)
    if(progress) pb <- utils::txtProgressBar(max = nsim, style = 3)
    for(i in 1:nsim){
      O <- Omega.draw(fit)
      b <- b.draw(data.t$b, data.t$long, data.t$surv, O, beta.inds, b.inds, fit)$b.hat
      pi.store[i] <- Surv_ratio(data.t$surv, data.u$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b)
      if(progress) utils::setTxtProgressBar(pb, i)
    }
    pi[[uu]] <- pi.store
  }
  return(lapply(pi, quantile, probs = c(0.500, 0.025, 0.975), na.rm = T))
}
