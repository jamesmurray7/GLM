#' ######
#' _Functions.R
#' ######

#' ########################################################################
# Data functions ----------------------------------------------------------
#' ########################################################################
 
#' Parsing the longitudinal formula
parseFormula <- function(formula){ 
  split <- glmmTMB:::splitForm(formula, allowFixedOnly = F)
  fixed <- split$fixedFormula
  random <- el(split$reTrmFormulas)
  
  #' Parse fixed effects
  response <- as.character(fixed)[2]
  fixed <- as.character(fixed)[3]
  if(grepl('splines\\:\\:|ns\\(|bs\\(', fixed)){
    attr(fixed, 'special') <- 'spline'
    if(grepl('ns\\(', fixed)) attr(fixed, 'spline type') <- 'natural'
    else if(grepl('bs\\(', fixed)) attr(fixed, 'spline type') <- 'basis'
    else stop('Unknown spline type')
  }else{
    attr(fixed, 'special') <- 'none'
  }
  # Flag splines as will likely be useful late wrhen setting up F_u, F_i, ftmat etc.^^^
  
  #' Parse random effects
  random.by <- as.character(random)[3]
  random <- as.character(random)[2]
  if(grepl('splines\\:\\:|ns\\(|bs\\(', random)){
    attr(random, 'special') <- 'spline'
    if(grepl('ns\\(', random)) attr(random, 'spline type') <- 'natural'
    else if(grepl('bs\\(', random)) attr(random, 'spline type') <- 'basis'
    else stop('Unknown spline type')
  }else{
    attr(random, 'special') <- 'none'
  }
  # Flag splines as will likely be useful later when setting up F_u, F_i, ftmat etc. ^^^
  
  return(list(
    response = response,
    fixed = fixed,
    random = random,
    random.by = random.by
  ))
}

#' Obtaining the data matrices X, Y, Z ------------------------------------
createDataMatrices <- function(data, formulas, disp.formula){
  fixed <- as.formula(paste0('~', formulas$fixed))
  response <- formulas$response
  random <- as.formula(paste0('~', formulas$random))
  
  #' Create dataset update if splines exit in model formula
  if(!is.null(attr(formulas$fixed, 'special')) & (attr(formulas$fixed, 'special') == 'spline' | attr(formulas$random, 'special') == 'spline')){
    newData <- as.data.frame(cbind(id = data$id, model.matrix(fixed, data)))
    newDataRandom <- as.data.frame(cbind(id = data$id, model.matrix(random, data)))
    n <- length(unique(data$id))
    X <- Z <- Y <- W <- vector('list', n)
    for(i in 1:n){
      X[[i]] <- as.matrix(newData[newData$id == i, -1, drop = F])
      Y[[i]] <- matrix(data[data$id == i, response], nc = 1)
      Z[[i]] <- as.matrix(newDataRandom[newDataRandom$id == i, -1, drop = F])
      W[[i]] <- as.matrix(model.matrix(disp.formula, newData))
    }
  }else{ #' Otherwise proceed as usual.
    n <- length(unique(data$id))
    X <- Z <- Y <- W <- vector('list', n)
    for(i in 1:n){
      i.dat <- subset(data, id == i)
      X[[i]] <- model.matrix(fixed, i.dat)
      Y[[i]] <- matrix(i.dat[, response], nc = 1)
      Z[[i]] <- model.matrix(random, i.dat)
      W[[i]] <- as.matrix(model.matrix(disp.formula, i.dat))
    }
    newData <- newDataRandom <- NULL
  }
  
  list(
    X = X, Y = Y, Z = Z, W = W,
    w = ncol(W[[1]]), p = ncol(X[[1]]), q = ncol(Z[[1]]),     # Return dimensions + the matrices.
    nw = colnames(W[[1]]), np = colnames(X[[1]]), nq = colnames(Z[[1]])
  )
}

#' Parsing the survival formula -------------------------------------------
parseCoxph <- function(surv.formula, data){
  survdata <- data[!duplicated(data[, 'id']), ]; n <- nrow(survdata)
  ph <- coxph(surv.formula, survdata, x = T)
  
  survdata <- data.frame(id = survdata$id, ph$x, survtime = ph$y[, 1], status = ph$y[, 2])
  pS <- ncol(ph$x)
  sf <- survfit(ph)
  ft <- sf$time[sf$n.event >= 1] # failure times
  nev <- sf$n.event[sf$n.event >= 1] # Number of failures per failure time

  Delta <- as.list(survdata$status)
  
  #' Return ----
  list(
    survdata = survdata, ph = ph, n = n, Delta = Delta
  )
}

#' Creation of survival objects -------------------------------------------
surv.mod <- function(ph, survdata, formulas, l0.init){
  uids <- unique(survdata$id); n <- length(uids)
  l0 <- l0.init; splines <- F
  
  # Failure times
  ft <- coxph.detail(ph)$time
  
  #' Unpack ph, data, formulas and construct data-related objects. ----
  #' Case 1:: No spline terms
  if((!is.null(attr(formulas$random, 'special'))) & attr(formulas$random, 'special') == 'spline'){
    splines <- T
    survdata$time <- unname(ph$y[,1])
    newSurvdata <- as.data.frame(cbind(id = survdata$id, model.matrix(as.formula(paste0('~', formulas$random)), survdata)))
    q <- ncol(newSurvdata) - 1 # (Intercept) INCLUDED as of now (22)
    .getFi <- function(time, q = q){
      id <- which(survdata$time == time)[1] # Take the first one in case of ties -- same design matrix regardless
      as.matrix(newSurvdata[id, -1, drop = F])
    } 
    .getFu <- function(times, q = q){
      as.matrix(newSurvdata[match(times, survdata$time), -1, drop = F])
    }
    Fi <- lapply(survdata$survtime, .getFi)
    ftmat <- .getFu(ft)
  }else{ #' Case 2:: Anything else // Maybe more specials in future(?)
    random <- formulas$random
    q <- length(el(strsplit(random, '\\+|\\*|\\:')))
    #' Define two functions to construct F_i and F_u ----
    .getFi <- function(time, q = q){
      Fi <- matrix(NA, nr = 1, nc = q)
      for(i in 1:q) Fi[, i] <- time^(i - 1)
      Fi
    }
    .getFu <- function(times, q = q){
      out <- sapply(1:q, function(i) times ^ (i - 1))
      if(!"matrix"%in%class(out)) out <- t(out)
      out
    }
    # Generate F_i and ftmat
    Fi <- lapply(survdata$survtime, .getFi, q)
    ftmat <- .getFu(ft, q)
  }
  
  # initialise empty stores
  l0i <- vector('numeric', n)
  Fu <- l0u <- surv.times <- vector('list', n)
  
  #' loop over subjects -----
  for(i in as.numeric(uids)){
    # slice this id
    survtime <- survdata[survdata$id == i, 'survtime']
    status <- survdata[survdata$id == i, 'status']
    # INDICES of survived times
    surv.times[[i]] <- which(ft <= survtime)
    # Fu, and l0u
    st <- ft[which(ft <= survtime)] # survived times
    if(length(st) > 0){
      Fu[[i]] <- .getFu(st, q)
      l0u[[i]] <- l0[which(ft <= survtime)]
    }else{ # Those who are censored before first failure time
      l0u[[i]] <- 0
      Fu[[i]] <- do.call(cbind, replicate(q, 0, simplify = F))
    }
    if(status == 1) l0i[i] <- l0[which(ft == survtime)] else l0i[i] <- 0
  }
  
  # Design matrix SS and rowvec S
  S <- lapply(1:n, function(i) ph$x[i, , drop = F])
  SS <- lapply(1:n, function(i){
    out <- apply(S[[i]], 2, rep, nrow(Fu[[i]]))
    if("numeric"%in%class(out)) out <- t(out)
    out
  })
  
  #' Return list ----
  return(list(
    ft = ft, ft.mat = ftmat, nev = coxph.detail(ph)$nevent, surv.times = surv.times,
    l0 = l0, l0i = as.list(l0i), l0u = l0u, Fi = Fi, Fu = Fu, Tis = survdata$survtime,
    S = S, SS = SS, q = q
  ))
}


# Taking difference -------------------------------------------------------
difference <- function(params.old, params.new, type){
  if(type == 'absolute'){
    rtn <- abs(params.new - params.old)
  }else if(type == 'relative'){
    rtn <- abs(params.new - params.old)/(abs(params.old) + 1e-3)
  }else{
    rtn <- NA
  }
  rtn
}

# Calculating the log-likelihood ------------------------------------------
log.lik <- function(Omega, dmats, b, surv, sv, l0u, l0i){
  # joint density - REs; Marginal ll of Y(s) added to f(T_i,\Delta_i|...).
  Y <- dmats$Y; X <- dmats$X; Z <- dmats$Z; lY <- lapply(Y, lfactorial); 
  W <- dmats$W
  
  beta <- Omega$beta; D <- Omega$D; zeta <- Omega$zeta; gamma <- Omega$gamma
  phi <- Omega$phi
  S <- sv$S; SS <- sv$SS; Delta <- surv$Delta
  Fu <- sv$Fu; Fi <- sv$Fi; 
  q <- dmats$q; w <- dmats$w; P <- dmats$p; n <- surv$n
  
  #' f(Y|...; Omega)
  ll.GP1 <- mapply(function(Y, X, Z, W, b){
    mu <- exp(X %*% beta + Z %*% b)
    phiv <- W %*% phi
    ll_genpois(mu, phiv, Y)
  }, Y = Y, X = X, Z = Z, W = W, b = b)
  
  #' f(T, Delta|...; Omega)
  ll.ft <- mapply(function(S, SS, Fi, Fu, b, Delta, l0i, l0u){
    temp <- if(Delta == 1) log(l0i) else 0.0
    # temp + Delta * (S %*% zeta + Fi %*% (gamma * b)) - 
    # crossprod(l0u, exp(SS %*% zeta + Fu %*% (gamma * b)))
    temp + Delta * S %*% zeta - crossprod(l0u, exp(SS %*% zeta))
  }, S = S, SS = SS, Fi = Fi, Fu = Fu, b = b, Delta = Delta, 
  l0i = l0i, l0u = l0u)
  
  #' log likelihood.
  ll.b <- mapply(function(b) mvtnorm::dmvnorm(b, sigma = D, log = T), b = b)
  ll <- sum(ll.GP1 + ll.ft + ll.b)
  
  # Number of observations
  N <- sum(sapply(Y, length))
  # Number of estimated parameters (n = # dispersion params)
  Ps <- length(gamma) + length(surv$ph$coefficients) + length(vech(D))
  df <- P + Ps + w                             #' All __estimated__ parameters
  df.residual <- N - (df + n * q)              #'     + all random effects
  
  # AIC and BIC
  aic <- -2 * ll + 2 * df
  bic <- -2 * ll + log(N) * df
  
  # Output
  structure(ll,
            nobs = N, n = n, 
            AIC = aic, BIC = bic, 
            df = df, df.residual = df.residual)
}

#' ########################################################################
# Misc functions ----------------------------------------------------------
#' ########################################################################

.any <- function(x, f) any(f(x))

plot.stepmat <- function(fit){
  plot(fit$stepmat, type = 'o', col = 'blue', ylab = 'Time (s)', xlab = 'Iteration #', pch = 20,
       main = paste0("EM took ", round(fit$EMtime + fit$postprocess.time, 2), 's'))
}

.print.loglik <- function(fit){
  ll <- fit$logLik
  df <- attr(ll, 'df')
  aic <- attr(ll, 'AIC')
  cat(sprintf("logLik: %.3f (df: %d), AIC: %.3f\n", ll, df, aic))
}

# Tabulate summary --------------------------------------------------------

my.summary <- function(myfit, printD = F){
  if(is.null(myfit$SE)) stop('Need to run EM with post.process = T')
  qz <- qnorm(.975)
  .to3dp <- function(x) round(x, 3)
  # Model fit info
  K <- 1
  responses <- myfit$modelInfo$forms$response
  families <- "Generalised Poisson"
  cat(paste0(responses, ' (', families, ')\n'))
  
  # Log-likelihood
  #' Print loglikelihood
  cat(.print.loglik(myfit))
  cat('\n')
  
  # Standard errors and parameter estimates.
  SE <- myfit$SE
  D <- myfit$co$D
  betas <- myfit$co$beta
  sigmas <- unlist(myfit$co$sigma)
  gammas <- myfit$co$gamma
  phi <- myfit$co$phi
  zetas <- myfit$co$zeta
  
  #' Random effects matrix
  if(printD){
    cat(paste0('Random effects variance-covariance matrix: \n'))
    print(.to3dp(D))
    cat('\n')
  }
  
  beta <- setNames(c(betas), row.names(betas))
  my <- c(beta, phi)
  
  rSE <- SE[grepl('^beta|^phi', names(SE))]
  lb <- my - qz * rSE; ub <- my + qz * rSE
  z <- my/rSE
  p <- 2 * (pnorm(abs(z), lower.tail = F))
  
  
  this.out <- setNames(data.frame(.to3dp(my), .to3dp(rSE), .to3dp(lb), .to3dp(ub), round(p, 3)),
                       c('Estimate', 'SE', '2.5%', '97.5%', 'p-value'))
  cat("Fixed effects:\n")
  print(this.out)
  cat('\n')
  
  # Longitudinal parts
  
  #' Survival
  cat('Event-time sub-model: \n')
  # Rename gammas?
  survs <- c(zetas, gammas)
  surv.SE <- SE[match(names(survs), names(SE))]
  
  lb <- survs - qz * surv.SE; ub <- survs + qz * surv.SE
  z <- survs/surv.SE
  p <- 2 * (pnorm(abs(z), lower.tail = F))
  
  surv.out <- setNames(data.frame(.to3dp(survs), .to3dp(surv.SE), .to3dp(lb), .to3dp(ub), round(p, 3)),
                       c('Estimate', 'SE', '2.5%', '97.5%', 'p-value'))
  print(surv.out)
  
  #' Elapsed times
  cat(sprintf('\nElapsed times as follows (%d iterations):\n', myfit$iter))
  print(.to3dp(myfit$elapsed.time))
  
  
  
  invisible(1+1)
}


