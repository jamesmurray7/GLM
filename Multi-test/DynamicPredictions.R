#' #####
#' DynamicPredictions.R
#' ----
#' All functions which help with dynamic probabilistic survival predictions
#' - Preparing longitudinal and survival data objects from an i.dat-type data.frame
#'                                                   and response/model info from a fit.
#' - Sampling from the parameter vector \Omega and density of random effects \b|\Y, T, \Delta;\Omega.
#' - Generate many simulated survival probabilities from time vector \u.
#' #####

# Prepare longitudinal data -----------------------------------------------
prepareLongData <- function(data, formulas, responses, u = NULL){
  #' Prepares data for one subject i. That is, data argument should be data[data$id==i, ] type
  #' Formulas/responses should come from fitted aEM object
  #' If u is NULL then the whole data[,] object is partitioned and returned.
  #'     One should use the `u` argument to project from Tstart in the Tstart+delta window.
  
  formulas <- lapply(formulas, parseFormula); K <- length(formulas)
  responses <- lapply(strsplit(responses, '\\s\\('), el, 1)
  
  if(!is.null(u)) data <- data[data$time <= u, ] # Truncate if necessary
  
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
prepareSurvData <- function(data, fit, formulas, u = NULL, hazard){
  #' data: i.dat-type object; survival data is then derived in-formula (by taking the first row).
  #' fit: joint model fit
  #' formulas: list of formulae from the fitted aEM object.
  #' u: horizon time; if not supplied then this is found for the entire profile is returned.
  #' hazard: hazard from fitted aEM object.
  
  formulas <- lapply(formulas, parseFormula); K <- length(formulas)
  
  #' Extract failure times and \lambda_0 from hazard
  ft <- hazard[,1]; l0 <- hazard[,2]
  
  survtime <- data$survtime[1]
  status <- as.integer(survtime <= max(data$time)) # Do they die in the window 0 < t / 0 < u?
  data$longtime <- data$time; max.longtime <- max(data$longtime)
  if(max.longtime == 0) max.longtime <- data$survtime # prevent issues
  data$time <- data$survtime # potentially bad move, we'll see!
  data <- data[1, ]
  # l0 <- predict(lm(l0 ~ splines::ns(ft, knots = quantile(data$longtime, probs = c(0.25, 0.75))))) # Use smoothed cubic spline in hazard for predictions
  
  if(!is.null(u)) data$time <- u # Don't extrapolate beyond time u
  if(is.null(u) && data$time > max.longtime) data$time <- max.longtime
  
  #' Fi
  Fi <- do.call(cbind, lapply(formulas, function(x){
    model.matrix(as.formula(paste0('~', x$random)), data)
  }))
  
  #' Fu and l0u
  # Work out survived times (and if necessary truncate at u if supplied)
  survived.times <- ft[which(ft <= data$time)]
  
  # print(survived.times)
  if(length(survived.times) > 0){ # Avoid complications if they're censored before first failure time.
    Fu <- do.call(cbind, lapply(formulas, function(x){
      model.matrix(as.formula(paste0('~', x$random)), data.frame(time=survived.times))
    }))
    l0u <- l0[1:length(survived.times)]
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
  Sname <- gsub("^zeta\\_", '', names(fit$coeffs$zeta))
  S <- as.matrix(data[1,Sname,drop=F])
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
  surv <- prepareSurvData(data, fit, fit$long.formulas, u = u, hazard = fit$hazard)
  
  b <- fit$RE[id,]; K <- length(fit$family)
  responsenames <- lapply(strsplit(fit$ResponseInfo, '\\s\\('), el , 1)
  
  gamma.rep <- rep(fit$coeffs$gamma, sapply(b.inds, length))
  
  if(is.null(u)){
    bfit <- optim(
        rep(0, ncol(fit$RE)), joint_density, joint_density_ddb, Y = long$Y, X = long$X, Z = long$Z,
        beta = fit$coeffs$beta, D = fit$coeffs$D, sigma = fit$coeffs$sigma, family = fit$family,
        Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i, SS = surv$SS, Fu = surv$Fu,
        haz = surv$l0u, gamma_rep = gamma.rep, zeta = fit$coeffs$zeta, beta_inds = beta.inds,
        b_inds = b.inds, K = K, method = 'BFGS', hessian = T
        )
  }else{
    bfit <- NULL
  }
  
  list(
    long = long, surv = surv, b.MLE = b, b = bfit
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
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){ # xform if necessary
    D <- as.matrix(Matrix::nearPD(D)$mat)
  }
  # D <- 0.5 * (D + t(D))
  # D <- as.matrix(Matrix::nearPD(D)$mat)
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

logLik.b <- function(b, long, surv, O, beta.inds, b.inds, fit){
  beta <- O$beta; D <- O$D; gamma <- rep(O$gamma, sapply(b.inds, length)); zeta <- O$zeta; sigma <- O$sigma
  
  neg.ll <- joint_density(b = b, Y = long$Y, X = long$X, Z = long$Z, beta = beta, D = D, sigma = sigma,
                          family = fit$family, Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i,
                          SS = surv$SS, Fu = surv$Fu, haz = surv$l0u, gamma_rep = gamma, zeta = zeta, beta_inds = beta.inds,
                          b_inds = b.inds, K = length(b.inds))
  
  -neg.ll
}

b.mh <- function(b.current, shift, long, surv, O, beta.inds, b.inds, fit, Sigma, b.density, df){
  accept <- 0
  if(b.density == 'normal'){
    b.prop <- MASS::mvrnorm(n = 1, mu = shift, Sigma = Sigma) # shift = rep(0, q) in nm case?
  }else{
    b.prop <- mvtnorm::rmvt(n = 1, sigma = Sigma, df = df, delta = shift)
  }
  
  # Joint data log likelihood on current and proposed values of b
  current.ll <- logLik.b(b.current, long, surv, O, beta.inds, b.inds, fit)
  prop.ll <- logLik.b(b.prop, long, surv, O, beta.inds, b.inds, fit)
  # Difference
  jointll.diff <- current.ll - prop.ll
  # Log-likelihoods
  # Faster using BLAS-C++ implementation of dmvnorm.
  if(b.density == 'normal'){ # N(b|Sigma)
    current.ll <- dmvnrm_arma_fast(t(b.current), shift, Sigma, T)
    prop.ll <- dmvnrm_arma_fast(t(b.prop), shift, Sigma, T)
  }else{
    current.ll <- dmvt_arma_fast(x = t(b.current), mean = shift, sigma = Sigma, df = df, logd = T)
    prop.ll <- dmvt_arma_fast(x = b.prop, mean = shift, sigma = Sigma, df = df, logd = T)
  }
  
  # Difference
  dens.diff <- current.ll - prop.ll
  
  # MH accept/reject scheme
  a <- exp(jointll.diff - dens.diff)
  if(is.nan(a)){
    cat(paste0('jointll.diff:', jointll.diff, '\n',
               'dens.diff:', dens.diff,'\n'))
    cat(paste0('b.prop: ', b.prop,'\nb.current: ', b.current,'\n',
               'joint b.prop: ', logLik.b(b.prop, long, surv, O, beta.inds, b.inds, fit), '\n',
               'joint b.current: ', logLik.b(b.current, long, surv, O, beta.inds, b.inds, fit), '\n'))
  }
  a <- min(1, a)
  u <- runif(1)
  if(u <= a){
    accept <- 1
    b.current <- b.prop
  }
  
  list(b.current = c(b.current), accept = accept)
}

# Dynamic predictions -----------------------------------------------------
# Dynamic survial predictions for ONE subject `id` across a vector of candidate times, u.
dynSurv <- function(data, id, fit, u = NULL, nsim = 200, progress = T,
                    b.density = 'normal', scale = NULL, df = NULL){
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
  
  u <- u[-1] # Don't want to find preds for first time T_{start}...
  
  pi <- setNames(vector('list', length = length(u)), paste0('u = ', u))
  MH.accept <- 0
  b.current <- shift <- data.t$b$par; Sigma <- solve(data.t$b$hessian)
  if(!is.null(scale)) Sigma <- Sigma * scale
  for(uu in seq_along(u)){
    data.u <- prepareData(newdata2, id = id, b.inds, beta.inds, fit = fit,  u = u[uu])
    
    pi.store <- numeric(nsim)
    if(progress) pb <- utils::txtProgressBar(max = nsim, style = 3)
    for(i in 1:nsim){
      O <- Omega.draw(fit)
      b.sim <- b.mh(b.current, shift, data.t$long, data.t$surv, O, beta.inds, b.inds, fit, Sigma, b.density, df)
      b.current <- b.sim$b.current
      # b.sim <- b.draw(b.current, shift, Sigma, data.t$long, data.t$surv, O, beta.inds, b.inds, fit)
      # b.current <- b.sim$b.current
      # cat('\n', b.current,'\n')
      MH.accept <- MH.accept + b.sim$accept
      pi.store[i] <- Surv_ratio(data.t$surv, data.u$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)
      if(progress) utils::setTxtProgressBar(pb, i)
    }
    pi[[uu]] <- pi.store
  }
  return(
    list(pi = lapply(pi, quantile, probs = c(0.500, 0.025, 0.975), na.rm = T),
         pi.mean = lapply(pi, mean, na.rm = T),
         MH.accept = MH.accept/nsim)
  )
}

# Dynamic survial predictions for ONE subject `id` across a vector of candidate times, u.
dynSurv2 <- function(data, id, fit, u = NULL, nsim = 200, progress = T,
                    b.density = 'normal', scale = NULL, df = NULL){
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
  
  u <- u[-1] # Don't want to find preds for first time T_{start}...
  
  pi <- structure(matrix(NA, nr = nsim, nc = length(u)),
                  dimnames = list(as.character(1:nsim), paste0('u=',round(u,3))))
  MH.accept <- 0
  b.current <- shift <- data.t$b$par; Sigma <- solve(data.t$b$hessian)
  if(!is.null(scale)) Sigma <- Sigma * scale
  if(progress) pb <- utils::txtProgressBar(max = nsim, style = 3)
  for(i in 1:nsim){
    O <- Omega.draw(fit)
    # b.sim <- b.mh(b.current, shift, data.t$long, data.t$surv, O, beta.inds, b.inds, fit, Sigma, b.density, df)
    # b.current <- b.sim$b.current
    b.sim <- b.draw(b.current, shift, Sigma, data.t$long, data.t$surv, O, beta.inds, b.inds, fit, b.density, df)
    b.current <- b.sim$b.current
    MH.accept <- MH.accept + b.sim$accept
    St <- Surv_(data.t$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)
    for(uu in seq_along(u)){
      # cat('uu:', uu, '; u[uu]:', u[uu],  # uncomment for loop debugging
      #     '\nb:', b.current,'.\n')
      data.u <- prepareData(newdata2, id = id, b.inds, beta.inds, fit = fit,  u = u[uu])
      pi[i, uu] <- Surv_(data.u$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)/(St)# + 1e-6)
      # cat('pi(uu):', 
      #     Surv_(data.u$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)/(St),
      #     '\n.pi[i, uu]:', pi[i,uu], '.\n')
    }
    if(progress) utils::setTxtProgressBar(pb, i)
  }
  pi.df <- data.frame(
    u = u,
    mean = colMeans(pi),
    median = apply(pi, 2, median),
    lower  = apply(pi, 2, quantile, prob = 0.025),
    upper  = apply(pi, 2, quantile, prob = 0.975),
    id = id
  )
  row.names(pi.df) <- NULL
  return(
    list(pi = pi.df,
         pi.raw = pi,
         MH.accept = MH.accept/nsim)
  )
} 

# Reciever Operator Characteristics (ROC) ---------------------------------
ROC <- function(fit, data, Tstart, delta, control = list()){
  #' Fit: Fitted joint model
  #' data: Data to which the joint model was fit
  #' Tstart: Time from which all survival probabilities should be calculated from.
  #' delta: interval to check for failures in.
  #' control: list of control arguments to be passed to DynSurv
  #' ...: Additional arguments to dynSurv
  
  # Parse control arguments. If nothing is supplied then default to the normal case.
  if(!is.null(control$b.density)) b.density <- control$b.density else b.density <- 'normal'
  if(!b.density%in%c('normal', 't')) stop("b.density must be either 'normal' or 't'")
  if(!is.null(control$scale)) scale <- control$scale else scale <- NULL
  if(b.density == 't' & is.null(scale)){
    message('Scale not supplied for t distribution, defaulting to * 2')
    scale <- 2
  }
  if(!is.null(control$df)) df <- control$df else df <- NULL
  if(b.density == 't' & is.null(df)){
    message('df not supplied for t distribution, defaulting to df = 4')
    df <- 4
  }
  if(!is.null(control$nsim)) nsim <- control$nsim else nsim <- 25 # set quite low as we need to get through (potentially many) id's!
  
  # Set out new data and remove IDs where only one longitudinal measurement is available as this causes issues in calculation 
  newdata <- data[data$survtime > Tstart, ] # subjects who are ALIVE at Tstart.
  if('factor'%in%class(newdata$id)) newdata$id <- as.numeric(as.character(newdata$id)) # this confuses tapply later on
  bad.ids <- as.numeric(names(which(with(newdata, tapply(time, id, function(x) length(unique(x)))) == 1)))
  newdata <- newdata[!newdata$id%in%bad.ids, ]
  alive.ids <- unique(newdata$id)
  n.alive <- length(alive.ids)
  
  # Set out candidate failure times (u)
  ft <- fit$hazard[, 1]; tmax <- max(ft)
  window <- c(Tstart + 1e-6, Tstart + 1e-6 + delta)
  if(window[2] > tmax) window[2] <- tmax
  candidate.u <- c(Tstart, ft[ft > window[1] & ft <= window[2]])
  
  # Loop over ids and failure times
  probs <- acceptance <- setNames(vector('list', length = length(alive.ids)),
                         paste0('id ', alive.ids))
  pb <- utils::txtProgressBar(max = length(alive.ids), style = 3)
  for(i in seq_along(alive.ids)){
    ds <- dynSurv2(newdata, alive.ids[i], fit, u = candidate.u, progress = F, 
                   b.density = b.density, scale = scale, df = df, nsim = nsim)
    probs[[i]] <- ds$pi#do.call(rbind, ds$pi)
    acceptance[[i]] <- ds$MH.accept
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Obtaining conditional probabilities for those alove subjects at Tstart.
  infodf <- lapply(alive.ids, function(x){
    p <- as.data.frame(probs[[paste0('id ', x)]])
    p$mean <- p$mean   
    p$id <- x
    p
  })
  # pi(u|t)
  infodf <- do.call(rbind, infodf)
  pi <- 1 - with(infodf, tapply(`mean`, id, min))
  
  # Working out whether individuals failed in the window
  survtimes <- with(newdata, tapply(survtime, id, unique))
  events <- survtimes >= window[1] & survtimes <= window[2]
  status <- as.logical(with(newdata, tapply(status, id, unique))) # Check if they failed
  event <- status & events # Check if they failed in window.
  
  n.window.events <- sum(event)
  BS <- (pi-sapply(event, as.numeric))^2   # Brier score
  
  # Defining threshold and calculating performance metrics.
  t <- seq(0, 1, length = 101)
  simfail <- structure(outer(pi, t, '<'),
                       dimnames = list(names(pi) , paste0('t: ', t)))
  
  TP <- colSums(c(event) * simfail)        # True positives
  FN <- sum(event) - TP                    # False negatives
  FP <- colSums(c(!event) * simfail)        # False positives
  TN <- sum(!event) - FP                   # True negatives
  TPR <- TP/(TP + FN)                      # True positive rate (sensitivity)
  FPR <- FP/(FP + TN)                      # False positive rate (1 - specificity)
  Acc <- (TP + TN) / (TP + TN + FP + FN)   # Accuracy
  PPV <- TP/(TP + FP + 1e-6)               # Positive predictive value (precision)
  NPV <- TN/(TN + FN + 1e-6)               # Negative predictive value
  F1s <- 2*(PPV* TPR) / (PPV + TPR + 1e-6) # F1 score
  
  # Sanity checks -- mainly for internal use.
  if(!all.equal(TPR, TP/sum(event))) stop('Something wrong: TP + FN != sum(event)')
  if(!all.equal(FP / (FP + TN)  , FP/(sum(!event)))) stop('Something wrong: FP + TN != sum(!event)')
  
  # Making a nice dataframe to report
  out <- data.frame(threshold = t,
                    TP = TP, TN = TN, FP = FP, FN = FN,
                    TPR = TPR, FPR = FPR, PPV = PPV, NPV = NPV,
                    Acc = Acc, F1 = F1s)
  row.names(out) <- NULL
  
  # Flip table so if multiple thresholds have same TPR/FPR then we take the largest threshold
  out <- out[order(out$threshold, decreasing = T), ]
  # Remove duplicated TPR/FPR
  out <- out[!duplicated(out[, c('TPR', 'FPR')]),]
  
  simulation.info <- paste0('b.density: ', b.density, ' with ', nsim, ' simulations per subject-failure time.')
  if(b.density == 't') simulation.info <- paste0(simulation.info, '\nt simulated using proposal covariance scale: ', scale, ' and df: ', df, '.\n')
  
  list(
    Tstart = Tstart, delta = delta, candidate.u = candidate.u,
    window.failures = n.window.events,
    Tstart.alive = n.alive,
    metrics = out, BrierScore = mean(BS),
    MH.acceptance.bar = mean(do.call(c, acceptance)),
    simulation.info = simulation.info
  )
}

#' Functionality to plot ROC ----
plotROC <- function(ROC, legend = F){
  TPR <- ROC$metrics$TPR; FPR <- ROC$metrics$FPR;
  plot(FPR, TPR,
       xlab = '1 - Specificity', ylab = 'Sensitivity',
       main = paste0('ROC curve for interval (', ROC$Tstart, ', ', ROC$Tstart + ROC$delta, ']'),
       type = 'l')
  abline(0, 1, lty = 3)
  if(legend){
    legend('bottomright', 
           paste0(ROC$Tstart.alive, ' at risk; ', ROC$window.failures, ' failures in interval.\n',
                  'AUC: ', round(AUC(ROC), 3),', Brier score: ', round(ROC$BrierScore, 3), '.'),
           bty = 'n', cex = .75)
  }
  invisible()
}


# Area under ROC curve (AUC), from ROC ------------------------------------
AUC <- function(ROC){
  TPR <- rev(ROC$metrics$TPR); FPR <- rev(ROC$metrics$FPR);
  auc <- sum(0.5 * diff(FPR) * (TPR[-1] + TPR[-length(TPR)]), na.rm = TRUE)
  auc
}


# Brier score -------------------------------------------------------------
BrierScore <- function(fit, data, Tstart, delta, ...){
  #' Fit: Fitted joint model
  #' data: Data to which the joint model was fit
  #' Tstart: Time from which all survival probabilities should be calculated from.
  #' delta: interval to check for failures in.
  #' ...: Additional arguments to dynSurv
  #' This is currently supposing it's correct to have one per subject, instead of one contribution
  #'                                                                PER subject PER x candidate.u combo.
  
  # Set out new data and remove IDs where only one longitudinal measurement is available as this causes issues in calculation 
  newdata <- data[data$survtime > Tstart, ] # subjects who are ALIVE at Tstart.
  bad.ids <- as.numeric(names(which(with(newdata, tapply(time, id, function(x) length(unique(x)))) == 1)))
  newdata <- newdata[!newdata$id%in%bad.ids, ]
  alive.ids <- unique(newdata$id)
  n.alive <- length(alive.ids)
  
  # Set out candidate failure times (u)
  ft <- fit$hazard[, 1]; tmax <- max(ft)
  window <- c(Tstart + 1e-6, Tstart + 1e-6 + delta)
  if(window[2] > tmax) window[2] <- tmax
  candidate.u <- c(Tstart, ft[ft >= window[1] & ft <= window[2]])
  
  # Loop over ids and failure times
  probs <- setNames(vector('list', length = length(alive.ids)),
                    paste0('id ', alive.ids))
  pb <- utils::txtProgressBar(max = length(alive.ids), style = 3)
  for(i in seq_along(alive.ids)){
    ds <- dynSurv(newdata, alive.ids[i], fit, u = candidate.u, progress = F, ...)
    probs[[i]] <- do.call(rbind, ds)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Work out whether alive subjects at Tstart failed in window [Tstart, Tstart + delta]
  events <- with(newdata, tapply(survtime, id, function(x){
    x >= window[1] && x <= window[2]   # Either failure or censor in window
  }))
  status <- as.logical(with(newdata, tapply(status, id, unique))) 
  event <- status & events # Specifically experienced failure.
  n.window.event <- sum(event)
  event <- sapply(event, as.numeric)
  
  # Obtaining conditional probabilities for those alove subjects at Tstart.
  infodf <- lapply(alive.ids, function(x){
    p <- as.data.frame(probs[[paste0('id ', x)]])
    p$id <- x
    p
  })
  # pi(u|t)
  pi <- with(do.call(rbind, infodf), tapply(`50%`, id, min))
  
  if(length(pi) != length(event)) stop("Vector of probabilities pi is not the same length as observed event vector...")
 
  BS <- ((1-pi) - event)^2    # (1-pi) since this [pi] is probability of survival.
  
  mean(BS)
}


b.draw <- function(b.current, b.hat.t, Sigma.t, long, surv, O, beta.inds, b.inds, fit, 
                   b.density, df){
  #' Metropolis-Hasting scheme to draw from distribution of f(\b|T_i^*>t;\Omega)
  #' This distribution can be either normal (using approximation) or t on supplied df.

  #' Unpack \Omega
  beta <- O$beta; D <- O$D; gamma <- rep(O$gamma, sapply(b.inds, length)); zeta <- O$zeta; sigma <- O$sigma
  b.density <- match.arg(b.density, c('normal', 't'), several.ok = F)
  
  if(b.density == 'normal'){
    #' Use optim to draw from f(\b,\Y,T,\Delta;\Omega^{\ell}). Only thing changing per simulation is \Omega
    #' Posterior mode and its variance at \Omega^{\ell}
    b.hat.l <- optim(
      b.current, joint_density, joint_density_ddb,
      Y = long$Y, X = long$X, Z = long$Z, beta = beta, D = D, sigma = sigma,
      family = fit$family, Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i,
      SS = surv$SS, Fu = surv$Fu, haz = surv$l0u, gamma_rep = gamma, zeta = zeta, beta_inds = beta.inds,
      b_inds = b.inds, K = length(b.inds), method = 'BFGS'
    )$par
    Sigma.hat.l <- solve(joint_density_sdb(b = b.hat.l, Y = long$Y, X = long$X, Z = long$Z, beta = beta, D = D, sigma = sigma,
                                           family = fit$family, Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i,
                                           SS = surv$SS, Fu = surv$Fu, haz = surv$l0u, gamma_rep = gamma, zeta = zeta, beta_inds = beta.inds,
                                           b_inds = b.inds, K = length(b.inds), eps = 1e-4))
    # Draw from N(b.hat.l, Sigma.hat.l)
    b.prop.l <- MASS::mvrnorm(n = 1, mu = b.hat.l, Sigma = Sigma.hat.l)
    # DMNVORM on b draws
    current.dens <- as.double(dmvnrm_arma_fast(t(b.current), b.hat.t, Sigma.t, T))
    prop.dens <- as.double(dmvnrm_arma_fast(t(b.prop.l), b.hat.t, Sigma.t, T))
  }else{
    #' Draw from shifted t distribution at \hat{b}, \hat{\Sigma} for subject|T_i>t
    b.prop.l <- mvtnorm::rmvt(n = 1, sigma = Sigma.t, df = df, delta = b.hat.t)
    current.dens <- as.double(dmvt_arma_fast(t(b.current), b.hat.t, Sigma.t, df = df))
    prop.dens <- as.double(dmvt_arma_fast(b.prop.l, b.hat.t, Sigma.t, df = df))
  }
  diff.dens <- current.dens - prop.dens # Difference in current - proposal log-likelihood.
  
  # Joint data log likelihood on current and proposed values of b
  current.joint.ll <- logLik.b(b.current, long, surv, O, beta.inds, b.inds, fit)
  proposed.joint.ll <- logLik.b(c(b.prop.l), long, surv, O, beta.inds, b.inds, fit)
  diff.joint.ll <-  proposed.joint.ll - current.joint.ll

  # Accept/reject scheme
  accept <- 0
  a <- exp(diff.joint.ll - diff.dens)
  if(is.nan(a)){
    cat(paste0('jointll.diff:', diff.joint.ll, '\n',
               'dens.diff:', diff.dens,'\n'))
    cat(paste0('b.prop: ', b.prop.l,'\nb.current: ', b.current,'\n',
               'joint b.prop: ', logLik.b(b.prop.l, long, surv, O, beta.inds, b.inds, fit), '\n',
               'joint b.current: ', logLik.b(b.current, long, surv, O, beta.inds, b.inds, fit), '\n'))
  }
  a <- min(1, a)
  u <- runif(1)
  if(u <= a){
    accept <- 1
    b.current <- b.prop.l
  }
  
  list(b.current = c(b.current), accept = accept)
}


# Predictions + plot for one individual -----------------------------------

dynPlot <- function(data, id, fit, Tstart, delta, nsim = 200, progress = T,
                    b.density = 't', scale = 2, df = 4){
  
  ft <- fit$hazard[,1]; tmax <- max(ft)
  candidate.u <- c(Tstart, ft[ft>Tstart & ft <= (Tstart + delta)])
  # Obtain pi matrix for subject
  ds <- dynSurv2(data = data, id = id, fit = fit, u = candidate.u, nsim = nsim,
                 progress = progress,
                 b.density = b.density, scale = scale, df = df)
  
  # Available longitudinal profile.
  lhs.data <- data[data$id == id & data$time <= Tstart, ]
  
  
  
  
}

