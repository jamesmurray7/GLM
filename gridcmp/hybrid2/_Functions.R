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
    X <- Z <- Y <- G <- vector('list', n)
    for(i in 1:n){
      X[[i]] <- as.matrix(newData[newData$id == i, -1, drop = F])
      Y[[i]] <- matrix(data[data$id == i, response], nc = 1)
      Z[[i]] <- as.matrix(newDataRandom[newDataRandom$id == i, -1, drop = F])
      G[[i]] <- as.matrix(model.matrix(disp.formula, data[data$id==i,]))
    }
  }else{ #' Otherwise proceed as usual.
    n <- length(unique(data$id))
    X <- Z <- Y <- G <- vector('list', n)
    for(i in 1:n){
      i.dat <- subset(data, id == i)
      X[[i]] <- model.matrix(fixed, i.dat)
      Y[[i]] <- matrix(i.dat[, response], nc = 1)
      Z[[i]] <- model.matrix(random, i.dat)
      G[[i]] <- as.matrix(model.matrix(disp.formula, i.dat))
    }
    newData <- newDataRandom <- NULL
  }
  
  list(
    X = X, Y = Y, Z = Z, G = G# , 
    #newData = newData, newDataRandom = newDataRandom
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

# Functions for expectations ----------------------------------------------
E.lfactorialY <- function(lambda, nu, Z, summax){ # mu, nu, vectors
  out <- matrix(0, nr = length(lambda), nc = summax)
  for(j in 1:summax){
    out[, j] <- lgamma(j) * exp((j-1) * log(lambda) - nu * lgamma(j) - Z)
  }
  apply(out, 1, sum)
}
#  # very diff. results for t->oo as summax is increased...
E.YlfactorialY <- function(lambda, nu, Z, summax){
  out <- matrix(0, nr = length(lambda), nc = summax)
  for(j in 1:summax){
    out[, j] <- exp(
      log(j - 1) + log(lgamma(j)) + (j - 1) * log(lambda) - nu * lgamma(j) - Z
    )
  }
  apply(out, 1, sum)
}

V.lfactorialY <- function(lambda, nu, Z, summax, B){
  out <- matrix(0, nr = length(lambda), nc = summax)
  for(j in 1:summax){
    out[, j] <- lgamma(j)^2 * exp((j-1)*log(lambda) - nu * lgamma(j) - Z)
  }
  apply(out, 1, sum) - B^2
}

calc.ABC <- function(mu, nu, lambda, summax){
  # lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
  Z <- calcZ(log(lambda), nu, summax)
  B <- E.lfactorialY(lambda, nu, Z, summax)
  A <- E.YlfactorialY(lambda, nu, Z, summax) - mu * B
  C <- V.lfactorialY(lambda, nu, Z, summax, B) # c is potentially needed in W2 matrix creation, remove if not!
  list(A = A, B = B, C = C)
}

calc2.ABC <- function(mu, nu, lambda, summax, tau, w, v){
  # lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
  Z <- calcZ(log(lambda), nu, summax)
  B <- E.lfactorialY(lambda, nu, Z, summax)
  rhs <- numeric(length(mu))
  for(l in 1:length(w)) rhs <- rhs + w[l] * mu * exp(tau * v[l]) * B
  A <- E.YlfactorialY(lambda, nu, Z, summax) - rhs # mu * B
  C <- V.lfactorialY(lambda, nu, Z, summax, B) # c is potentially needed in W2 matrix creation, remove if not!
  list(A = A, B = B, C = C)
}

# Score for \beta and \delta ----------------------------------------------
# Sbeta <- function(beta, X, Y, Z, G, b, delta, summax){
#   mu <- exp(X %*% beta + Z %*% b)
#   nu <- exp(G %*% delta)
#   lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
#   V <- V_mu_lambda(mu, lambda, nu, summax)
#   if(length(mu) > 1) lhs <- diag(c(mu)) else lhs <- diag(mu)
#   crossprod(lhs %*% X, ((Y-mu) / V))
# }

Sbeta <- function(X, Y, mu, nu, lambda, V){
  if(length(mu) > 1) lhs <- diag(c(mu)) else lhs <- diag(mu)
  crossprod(lhs %*% X, ((Y-mu) / V))
}

# This the same as forward differencing, and probably slightly more preferrable given it's analytic!
# getW1 <- function(X, Z, G, b, beta, delta, summax){
#   mu <- exp(X %*% beta + Z %*% b)
#   nu <- exp(G %*% delta)
#   lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
#   V <- V_mu_lambda(mu, lambda, nu, summax)
#   if(length(mu) > 1) lhs <- diag(c(mu^2)/c(V)) else lhs <- diag(mu^2/V)
#   -crossprod(X, lhs) %*% X
# }

getW1 <- function(X, mu, nu, lambda, V){
  if(length(mu) > 1) lhs <- diag(c(mu^2)/c(V)) else lhs <- diag(mu^2/V)
  -crossprod(X, lhs) %*% X
}


# Correct, but no longer in use...
Sdelta <- function(delta, X, Y, lY, Z, b, G, beta, summax){
  mu <- exp(X %*% beta + Z %*% b)
  nu <- exp(G %*% delta)
  lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
  V <- V_mu_lambda(mu, lambda, nu, summax)
  AB <- calc.ABC(mu, nu, lambda, summax)
  Snu <- AB$A * (Y-mu) / V - (lY-AB$B)
  if(length(nu) > 1) lhs <- diag(c(nu)) else lhs <- diag(nu)
  crossprod(lhs %*% G, Snu)
}

# Taking difference -------------------------------------------------------
difference <- function(params.old, params.new, type){
  if(type == 'absolute'){
    rtn <- max(abs(params.new - params.old))
  }else if(type == 'relative'){
    rtn <- max(
      abs(params.new - params.old)/(abs(params.old) + 1e-3)
    )
  }else{
    rtn <- NA
  }
  rtn
}


# Calculating the log-likelihood ------------------------------------------
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

# delta initial conditions ------------------------------------------------


#' ########################################################################
# Misc functions ----------------------------------------------------------
#' ########################################################################

.safevar <- function(x) ifelse(length(x)>1, var(x, na.rm =T), 1)
.summax <- function(x) ceiling(max(x) + 20 * sqrt(.safevar(x)))
.any <- function(x, f) any(f(x))

plot.stepmat <- function(fit){
  plot(fit$stepmat, type = 'o', col = 'blue', ylab = 'Time (s)', xlab = 'Iteration #', pch = 20,
       main = paste0("EM took ", round(fit$EMtime + fit$postprocess.time, 2), 's'))
  abline(h = median(fit$stepmat), lty = 5, col = 'red')
}



# Expectation and Variance of Y -------------------------------------------
  # (https://documentation.sas.com/doc/en/pgmsascdc/v_017/ds2ref/n1bwtppmvqmdlgn1rx99e8qchvof.htm)

Ey <- function(mu, nu, summax){
  lam <- lambda_appx(mu, nu, summax)
  logZ <- logZ_c(log(lam), nu, summax)
  rhs <- matrix(0, nr = length(mu), nc = summax)
  for(j in 1:summax){
    rhs[,j] <- exp(log(j-1) + (j-1) * log(lam) - nu * lgamma(j) - logZ)
  }
  apply(rhs, 1, sum)
}

VarY <- function(mu, nu, summax, Ey){
  lam <- lambda_appx(mu, nu, summax)
  logZ <- logZ_c(log(lam), nu, summax)
  rhs <- matrix(0, nr = length(mu), nc = summax)
  for(j in 1:summax){
    rhs[,j] <- exp(
      2 * log(j - 1) + (j - 1) * log(lam) - nu * lgamma(j) - logZ
    )
  }
  apply(rhs, 1, sum) - (Ey^2) # in cases where delta = 0, this will just equal E[Y].
}

#  Seek to minimise difference between real and expected variance.
ff <- function(delta, Y, G, mu, summax){
  nu <- exp(G %*% delta)
  EY <- Ey(mu, nu, summax)
  VY <- VarY(mu, nu, summax, EY)
  var(Y) - mean(VY)
}

ff2 <- function(delta, Y, G, mu, summax){
  nu <- exp(G %*% delta)
  # EY <- Ey(mu, nu, summax)
  # VY <- VarY(mu, nu, summax, EY)
  VY <- expected_variance(mu, nu, summax)
  sum(c(var(Y)) - VY)
}

# Determine interval in the 1-d case.
determine.interval <- function(mu, summax){
  m <- length(mu)
  lower <- -3; upper <- 3 # Nominally choose lower and upper values of interval to check over.
  test.lower <- lambda_appx(mu, rep(exp(lower), m), summax)
  test.upper <- lambda_appx(mu, rep(exp(upper), m), summax)
  
  flag.lower <- any(is.nan(test.lower)); flag.upper <- any(is.nan(test.upper))
  while(flag.lower){
    lower <- lower + 0.2
    flag.lower <- any(is.nan(lambda_appx(mu, rep(exp(lower), m), summax)))
  }
  while(flag.upper){
    upper <- upper - 0.2
    flag.upper <- any(is.nan(lambda_appx(mu, rep(exp(upper), m), summax)))
  }
  c(lower = lower, upper = upper)
}


delta.uniroot <- function(Y, G, mu, summax){
  .interval <- determine.interval(mu, summax)
  tryCatch(uniroot(ff2, interval = c(.interval[1], .interval[2]), Y = Y, G = G, mu = mu, summax = summax, extendInt = 'yes')$root,
           error = function(e) NA)
}

find.deltas <- function(Ylist, Glist, mulist, summax, verbose = F){
  candidateY <- sapply(Ylist, length) > 1
  numY <- sum(candidateY); inds <- which(candidateY)
  out <- numeric(numY); p <- 1
  for(i in inds){
    out[p] <- delta.uniroot(Ylist[[i]], Glist[[i]], mulist[[i]], summax)
    if(verbose) cat(sprintf('%d/%d\r', p, numY))
    p <- p + 1
  }
  out
}

# Optim version ----
ff3 <- function(delta, Y, G, mu, summax){
  nu <- exp(G %*% delta)
  VY <- expected_variance(mu, nu, summax)
  abs(sum(c(var(Y)) - VY))
}

# Determine interval for the n-dimensional case...
# ... (wip)
determine.interval.nd <- function(Y, G, mu, summax){
  g <- ncol(G)
  
  # we want g permutations of seq(-3,3,.2)
  candidates <- seq(-3, 3, .2)
  candidates.mat <- structure(t(combn(candidates, g)),
                              dimnames = list(as.character(1:choose(length(candidates),g)),
                                              paste0('delta', 1:g))
                              )

  vars <- sort(apply(candidates.mat, 1, function(i){
    ff2(c(i), Y[[1]], G[[1]], mus[[1]], summax)
  }))
}

delta.optim <- function(Y, G, mu, summax){
  g <- ncol(G)
  if(g == 1){
    .interval <- determine.interval(mu, summax)
    out <- tryCatch(optim(c(0), ff3, NULL,
                   Y = Y, G = G, mu = mu, summax = summax,
                   method = 'Brent', lower = .interval[1], upper = .interval[2])$par,
             error = function(e) NA)
  }else{
    out <- tryCatch(optim(rep(0, g), ff3, NULL,
                   Y = Y, G = G, mu = mu, summax = summax,
                   method = 'L-BFGS-B', lower = rep(-3, g), upper = rep(3, g))$par,
             error = function(e) NA)
  }
  out
}

find.deltas.optim <- function(Ylist, Glist, mulist, summax, verbose = F){
  candidateY <- sapply(Ylist, length) > 1
  numY <- sum(candidateY); inds <- which(candidateY)
  out <- numeric(numY); p <- 1
  for(i in inds){
    out[p] <- delta.optim(Ylist[[i]], Glist[[i]], mulist[[i]], summax)
    if(verbose) cat(sprintf('%d/%d\r', p, numY))
    p <- p + 1
  }
  out
}
