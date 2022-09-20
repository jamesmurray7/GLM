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
createDataMatrices <- function(data, formulas){
  fixed <- as.formula(paste0('~', formulas$fixed))
  response <- formulas$response
  random <- as.formula(paste0('~', formulas$random))
  
  #' Create dataset update if splines exit in model formula
  if(!is.null(attr(formulas$fixed, 'special')) & (attr(formulas$fixed, 'special') == 'spline' | attr(formulas$random, 'special') == 'spline')){
    newData <- as.data.frame(cbind(id = data$id, model.matrix(fixed, data)))
    newDataRandom <- as.data.frame(cbind(id = data$id, model.matrix(random, data)))
    n <- length(unique(data$id))
    X <- Z <- Y <- G <- vector('list', n)
    for(i in 1:n){
      X[[i]] <- as.matrix(newData[newData$id == i, -1, drop = F])
      Y[[i]] <- matrix(data[data$id == i, response], nc = 1)
      Z[[i]] <- as.matrix(newDataRandom[newDataRandom$id == i, -1, drop = F])
    }
  }else{ #' Otherwise proceed as usual.
    n <- length(unique(data$id))
    X <- Z <- Y <- G <- vector('list', n)
    for(i in 1:n){
      i.dat <- subset(data, id == i)
      X[[i]] <- model.matrix(fixed, i.dat)
      Y[[i]] <- matrix(i.dat[, response], nc = 1)
      Z[[i]] <- model.matrix(random, i.dat)
    }
    newData <- newDataRandom <- NULL
  }
  
  list(
    X = X, Y = Y, Z = Z
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


# Functions for finding \lambda from \mu and \nu --------------------------
# DEFUNCT AS OF JUNE 2022, APART FROM USE IN SIMDATA.
.getlambda <- function(mui, nui, summax){
  tryCatch(uniroot(mu_lambdaZ_eq, interval = c(1e-6, 1e3), mu = mui, nu = nui, summax = summax)$root,
           error = function(e) NA)
}

getlambda <- function(mu, nu, summax){
  # i. An approximation 
  loglambdas.appx <- suppressWarnings(
    nu * log(mu + (nu - 1) / (2 * nu))
  )
  
  lambdas.appx <- exp(loglambdas.appx)
  
  # ii. Find solutions to mean constraint (Huang (2017)) and clean so NaN/Infs/NAs not in output.
  lambdas <- mapply(function(mu, nu, lambdas.appx){
    out <- .getlambda(mu, nu, summax)
    # If uniroot fails to find a root, set it as the approximation above
    if((is.na(out) | is.nan(out)) & (!is.nan(lambdas.appx) & !is.na(lambdas.appx))) out <- lambdas.appx
    # And if this still NA/NaN/Inf, simply set as mean
    if(is.na(out)) out <- mu
    out
  }, mu = mu, nu = nu, lambdas.appx = lambdas.appx, SIMPLIFY = T)
  
  # Print how many rate parameters simply used the mean.
  sprintf('%.2f%% values used mean', length(which(lambdas == mu))/length(mu) * 100)
  lambdas
}

Sbeta <- function(X, Y, mu, nu, lambda, V){
  if(length(mu) > 1) lhs <- diag(c(mu)) else lhs <- diag(mu)
  crossprod(lhs %*% X, ((Y-mu) / V))
}

getW1 <- function(X, mu, nu, lambda, V){
  if(length(mu) > 1) lhs <- diag(c(mu^2)/c(V)) else lhs <- diag(mu^2/V)
  -crossprod(X, lhs) %*% X
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
# NYI!!!!!!!!!!!!
log.lik <- function(coeffs, dmats, b, surv, sv, l0u, l0i, summax){
  # joint density - REs; Marginal ll of Y(s) added to f(T_i,\Delta_i|...).
  Y <- dmats$Y; X <- dmats$X; Z <- dmats$Z; G <- dmats$G; lY <- lapply(Y, lfactorial)
  beta <- coeffs$beta; D <- coeffs$D; delta <- coeffs$delta; zeta <- coeffs$zeta; gamma <- coeffs$gamma
  S <- sv$S; SS <- sv$SS; Delta <- surv$Delta
  Fu <- sv$Fu; Fi <- sv$Fi; 
  q <- ncol(Z[[1]]); K <- 1; n <- length(Z)
  
  #' "Full" joint density ----
  ll1 <- numeric(n)
  for(i in 1:n){ # For some reason mapply() not liking this as it was in Multi-test(?)
    ll1[i] <- -joint_density(b[[i]],X[[i]],Y[[i]],lY[[i]],Z[[i]],G[[i]],beta, delta,D,S[[i]],SS[[i]],Fi[[i]], Fu[[i]],l0i[[i]],l0u[[i]],Delta[[i]],gamma,zeta,100)
  }
  
  #' Only the joint density associated with REs ----
  ll2 <- mapply(function(b) mvtnorm::dmvnorm(b, rep(0, q), D, log = T), b = b)
  #' Calculate the marginal ll of the Yks and the survival density.
  out <- sum(ll1 - ll2)
  
  #' Calculate AIC and BIC
  N <- sum(sapply(Y, length))
  P <- ncol(X[[1]]-1)
  Ps <- ncol(S[[1]])
  df <- P + Ps + 2 * K + (q * (q + 1)/2)
  
  attr(out, 'df') <- df; attr(out, 'N') <- N
  attr(out, 'AIC') <- -2 * c(out) + 2 * df
  attr(out, 'BIC') <- -2 * c(out) + log(N) * df
  out
}


# Wrapper for minimisation of f(Y|.) wrt b --------------------------------
b.minimise <- function(b, X, Y, lY, Z, delta, S, SS, Fi, Fu, l0i, l0u, Delta, 
                       Omega, summax, method = 'optim', obj = 'joint_density', Hessian = 'grad', Hessian.epsilon){
  if(!obj%in%c('joint_density', 'marginal')) 
    stop("'obj' must be either the 'joint_density' or the 'marginal' Y|~CMP(.).\n")
  if(!method%in%c('optim', 'bobyqa'))
    stop("'method' must be either 'optim' or 'bobyqa'.\n")
  if(!Hessian%in%c('grad', 'obj'))
    stop("'Hessian' must be either 'grad' (forward differencing on gradient function) or 'obj' (forward differencing on objective function).")
  
  beta <- Omega$beta; D <- Omega$D; gamma <- Omega$gamma; zeta <- Omega$zeta; 
  
  if(obj == 'joint_density'){
    if(method == 'optim'){
      b.hat <- mapply(function(b, X, Y, lY, Z, delta, S, SS, Fi, Fu, l0i, l0u, Delta, summax){
        optim(b, joint_density, joint_density_ddb,
              X = X, Y = Y, lY = lY, Z = Z, beta = beta, delta = delta, D = D,
              S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
              gamma = gamma, zeta = zeta, summax = summax, method = 'BFGS', hessian = F)$par
      }, b = b, X = X, Y = Y, lY = lY, Z = Z, delta = delta, S = S, SS = SS, Fi = Fi, Fu = Fu,
      l0i = l0i, l0u = l0u, Delta = Delta, summax = summax, SIMPLIFY = F)
    }else{
      b.hat <- mapply(function(b, X, Y, lY, Z, delta, S, SS, Fi, Fu, l0i, l0u, Delta, summax){
        nloptr::bobyqa(b, joint_density,
                       X = X, Y = Y, lY = lY, Z = Z, beta = beta, delta = delta, D = D,
                       S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                       gamma = gamma, zeta = zeta, summax = summax)$par
      }, b = b, X = X, Y = Y, lY = lY, Z = Z, delta = delta, S = S, SS = SS, Fi = Fi, Fu = Fu,
      l0i = l0i, l0u = l0u, Delta = Delta, summax = summax, SIMPLIFY = F)
    }
    
    #' Calculate the inverse of the second derivative of the negative log-likelihood ----
    if(Hessian == 'obj'){
      Sigma <- mapply(function(b, X, Y, lY, Z, delta, S, SS, Fi, Fu, l0i, l0u, Delta, summax){
        out <- solve(H_joint_density(b = b, X = X, Y = Y, lY = lY, Z = Z, beta = beta, delta = delta, D = D,
                              S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                              gamma = gamma, zeta = zeta,
                              summax = summax, eps = Hessian.epsilon))
        if(det(out) <= 0 && Hessian.epsilon != .Machine$double.eps^(1/3)){ 
          # 1/4 __SHOULD__ work majoirty of time; fail-safe coded in case it doesn't.
          out <- solve(H_joint_density(b = b, X = X, Y = Y, lY = lY, Z = Z, beta = beta, delta = delta, D = D,
                                       S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                                       gamma = gamma, zeta = zeta,
                                       summax = summax, eps = .Machine$double.eps^(1/3)))
        }
        out
      }, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, delta = delta, S = S, SS = SS, Fi = Fi, Fu = Fu,
      l0i = l0i, l0u = l0u, Delta = Delta, summax = summax, SIMPLIFY = F)
    }else{
      Sigma <- mapply(function(b, X, Y, lY, Z, delta, S, SS, Fi, Fu, l0i, l0u, Delta, summax){
        solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, beta = beta, delta = delta, D = D,
                                S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                                gamma = gamma, zeta = zeta, summax = summax, eps = Hessian.epsilon))
      }, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, delta = delta, S = S, SS = SS, Fi = Fi, Fu = Fu,
      l0i = l0i, l0u = l0u, Delta = Delta, summax = summax, SIMPLIFY = F)
    }
    
    
  }else{
    if(method == 'optim'){
      b.hat <- mapply(function(b, X, Y, lY, Z, delta, summax){
        optim(b, marginal_Y, marginal_Y_db,
              X = X, Y = Y, lY = lY, Z = Z, beta = beta, delta = delta, D = D,
              summax = summax, method = 'BFGS')$par
      }, b = b, X = X, Y = Y, lY = lY, Z = Z, delta = delta, summax = summax, SIMPLIFY = F)
    }else{
      b.hat <- mapply(function(b, X, Y, lY, Z, delta, summax){
        nloptr::bobyqa(b, marginal_Y,
                       X = X, Y = Y, lY = lY, Z = Z, beta = beta, delta = delta, 
                       D = D,summax = summax)$par
      }, b = b, X = X, Y = Y, lY = lY, Z = Z, delta = delta, summax = summax, SIMPLIFY = F)
    }
  }
  
  list(
    b.hat = b.hat,
    Sigma = if(obj == 'joint_density') Sigma else NULL
  )
  
}


#' ########################################################################
# Misc functions ----------------------------------------------------------
#' ########################################################################

.safevar <- function(x) ifelse(length(x)>1, var(x, na.rm =T), 1)
.summax <- function(x) ceiling(max(x) + 20 * sqrt(.safevar(x)))
.any <- function(x, f) any(f(x))

plot.stepmat <- function(fit){
  plot(fit$stepmat, type = 'o', col = 'blue', ylab = 'Time (s)', xlab = 'Iteration #', pch = 20,
       main = paste0("EM took ", round(fit$EMtime + fit$postprocess.time, 2), 's'))
}

plot.delta.inits <- function(fit, show.medians = F, show.means = F){
  
  # Extract delta info
  x <- fit$modelInfo$delta.init
  s <- x$subject.estimates
  # Make 'cut' object
  cuts <- s[round(abs(s), 3) < 2 & !is.na(s)]
  
  # Densities and plot limits
  d1 <- density(s); d2 <- density(cuts)
  xlims <- c(min(d1$x, d2$x), max(d1$x, d2$x))
  ylims <- c(min(d1$y, d2$y), max(d1$y, d2$y))
  
  # The plot
  plot(d1, main = '', xlab = expression(delta), xlim = xlims, ylim = ylims)
  lines(d2, col = 'red')
  
  # (optional vertical lines)
  if(show.means)
    abline(v = c(x$mean.estimate, x$mean.cut.estimate), col = c('black', 'red'))
  if(show.medians)
    abline(v = c(x$median.estimate, x$median.cut.estimate), col = c('black', 'red'))
  
  # The legend
  legend('topleft', col = c('black', 'red'), lty = c(1, 1),
         legend = c('No cuts', 'Cuts'), bty = 'n')
}


# Tabulate summary --------------------------------------------------------

my.summary <- function(myfit, printD = F){
  if(is.null(myfit$SE)) stop('Need to run EM with post.process = T')
  qz <- qnorm(.975)
  .to3dp <- function(x) round(x, 3)
  # Model fit info
  K <- 1
  responses <- myfit$modelInfo$forms$response
  families <- "mean parameterised Conway-Maxwell Poisson"
  # Standard errors and parameter estimates.
  SE <- myfit$SE
  D <- myfit$co$D
  betas <- myfit$co$beta
  sigmas <- unlist(myfit$co$sigma)
  gammas <- myfit$co$gamma
  delta <- myfit$co$delta
  zetas <- myfit$co$zeta
  
  #' Random effects matrix
  if(printD){
    cat(paste0('Random effects variance-covariance matrix: \n'))
    print(.to3dp(D))
    cat('\n')
  }
  
  beta <- setNames(c(betas), row.names(betas))
  my <- c(beta, delta)
  
  rSE <- SE[grepl('^beta|^delta', names(SE))]
  lb <- my - qz * rSE; ub <- my + qz * rSE
  z <- my/rSE
  p <- 2 * (pnorm(abs(z), lower.tail = F))
  
  cat(paste0(responses, ' (', families, ')\n'))
  this.out <- setNames(data.frame(.to3dp(my), .to3dp(rSE), .to3dp(lb), .to3dp(ub), round(p, 3)),
                       c('Estimate', 'SE', '2.5%', '97.5%', 'p-value'))
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
  
  #' Optimiser info
  cat('\n')
  a <- myfit$modelInfo$optim
  cat(sprintf('Optimiser used: %s, Hessian appraised: %s with epsilon %s.\n', a[1], a[2], a[3]))
  invisible(1+1)
}


# Compare fitted delta with true ------------------------------------------
# (obviously only used in sim. studies!)
fitted.disps <- function(myfit, S, return.df = F, plot = T){
  # S: simdata_joint(2) object list
  Sdf <- data.frame(id = as.numeric(names(S$true.deltas)), true = S$true.deltas)
  # myfit: EM object
  df <- myfit$delta.df
  # merge
  df2 <- merge(df, Sdf, 'id')
  
  # 1) Check that the SIGN matches.
  df2$same.sign <- with(df2, sign(delta) == sign(true))
  # 2) Check coverage of true value by lower/upper bounds.
  df2$coverage <- with(df2, lb <= true & ub >= true)
  
  cat(sprintf("Sign matches %.2f%%, 95%% coverage %.2f.\n",
              sum(df2$same.sign)/nrow(df2) * 100,
              sum(df2$coverage)/nrow(df2)))
  
  if(plot){
    par(mfrow = c(1,2))
    uq.deltas <- unique(Sdf$true)
    for(i in seq_along(uq.deltas)){
      dfi <- df2[df2$true == uq.deltas[i],]
      plot(dfi$delta, pch = 20,
           xaxt = 'n', xlab = '', ylab = 'Estimate', 
           main = bquote('True delta:' ~ .(uq.deltas[i])))
      abline(h = 0, col = 'lightgrey', lty = 5)
      abline(h = uq.deltas[i], col = 'red', lty = 5)
    }
    par(mfrow = c(1,1))
  }
 
  if(return.df) return(df2) else return(invisible(1+1))
}




