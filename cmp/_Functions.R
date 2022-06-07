#' ######
#' _Functions.R
#' ---
#' Important functions not pertaining to another specific 
#' file
#' ######


# Functions for finding \lambda from \mu and \nu --------------------------

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


# Functions for taking score via quadrature -------------------------------

Score_delta <- function(delta, Y, X, Z, b, G, tau){
  nu <- exp(G %*% delta)
  lhs <- numeric(length(Y))
  for(l in 1:.ngh){
    mu <- exp(X %*% beta + Z %*% b + v[l] * tau)
    lambda <- getlambda(mu, nu, 10)
    log.lambda <- log(lambda)
    logZ_ <- logZ_c(log.lambda, nu, 10)
    lfY <- E_logfactY(lambda, nu, logZ_, 10)
    A <- E_YlogfactY(lambda, nu, logZ_, 10) - lfY * mu
    B <- lfY
    V <- V_mu_lambda(mu, lambda, nu, 10)
    lhs <- lhs + w[l] * Score_nu(A, Y, mu, V, B)
  }
  crossprod(lhs * nu, G)
}

Score_beta_quad <- function(beta, Y, X, Z, b, nu, tau){
  lhs <- numeric(length(Y))
  for(l in 1:.ngh){
    mu <- exp(X %*% beta + Z %*% b + tau * v[l])
    lam <- getlambda(mu, nu, 10)
    V <- V_mu_lambda(mu, lam, nu, 10)
    lhs <- lhs + w[l] * mu * (Y - mu) / V
  }
  crossprod(lhs, X) 
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

calc.AB <- function(mu, nu, lambda, summax){
  # lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
  Z <- calcZ(log(lambda), nu, summax)
  B <- E.lfactorialY(lambda, nu, Z, summax)
  A <- E.YlfactorialY(lambda, nu, Z, summax) - mu * B
  list(A = A, B = B)
}


