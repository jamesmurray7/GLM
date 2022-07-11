#' Copies from Thomas Fung's mpcmp package.
#' These are copied as tha package is not available for my version of R
library(Rcpp)
library(RcppArmadillo)
cppFunction('NumericVector logZ_c(NumericVector log_lambda, NumericVector nu, int summax) {
  // Control loop
  // int maxiter = 1e4;
  double log_epsilon = std::log(1e-10);
  // Output vector
  int n = log_lambda.size();
  NumericVector out(n);
  // Compute logz
  for (int i = 0; i < n; ++i) {
    double logz  = 0;
    double logz_ = 0;
    for (int j = 1; j < summax; ++j) {
      logz_ += log_lambda[i] - nu[i] * log((double)j);
      logz = R::logspace_add(logz, logz_);
      if (logz_ - logz < log_epsilon) break;
    }
    out[i] = logz;
  }
  return out;
}')

logZ <- function(log_lambda, nu, summax = 100) {
  # approximates normalizing constant for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(log_lambda = log_lambda, nu = nu)
  log_lambda <- df[, 1]
  nu <- df[, 2]
  return(logZ_c(log_lambda, nu, summax = summax))
}

CBIND <- function(..., deparse.level = 1) {
  dots <- list(...)
  len <- sapply(dots, length)
  dots <- lapply(seq_along(dots),
                 function(i, x, len) rep(x[[i]], length.out = len),
                 x = dots, len = max(len)
  )
  do.call(cbind, c(dots, deparse.level = deparse.level))
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

dcomp <- function(x, mu, nu = 1, lambda, log.p = FALSE, lambdalb = 1e-10,
                  lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6, summax) {
  # compute the pmf/density for COMP distirbution with mean mu and dispersion nu
  # x, mu, nu are recycled to match the length of each other.
  # lambdaub will be scaled down/up  if there is
  # over-/under-dispersion so that the correct lambda can be found
  if (missing(mu) && missing(lambda)) {
    stop('argument "mu" is missing, with no default')
  }
  if (!missing(mu) && !missing(lambda)) {
    stop("specify 'mu' or 'lambda' but not both")
  }
  not.miss.mu <- !missing(mu)
  if (missing(mu)) {
    mu <- Inf
  }
  if (missing(lambda)) {
    lambda <- Inf
  }
  df <- CBIND(x = x, mu = mu, nu = nu, lambda = lambda)
  x <- df[, 1]
  mu <- df[, 2]
  nu <- df[, 3]
  lambda <- df[, 4]
  warn <- FALSE
  if (not.miss.mu) {
    if (missing(summax)) {
      summax <- ceiling(max(c(mu + 20 * sqrt(mu / nu), 100)))
    }
    mu.ok.ind <- which(mu > 0)
    mu.err.ind <- which(mu <= 0)
    if (length(mu.err.ind) > 0) {
      lambda[mu.err.ind] <- mu[mu.err.ind]
    }
    if (length(mu.ok.ind) > 0) {
      lambda.ok <- comp_lambdas(mu[mu.ok.ind], nu[mu.ok.ind],
                                lambdalb = lambdalb,
                                lambdaub = lambdaub,
                                # lambdaub = min(lambdaub,2*max(lambdaold))
                                maxlambdaiter = maxlambdaiter, tol = tol,
                                summax = summax
      )
      lambda[mu.ok.ind] <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
    }
  } else {
    # A <- (8*nu^2+12*nu+3)/(96*nu^2*lambda^(1/nu))
    # B <- (1+6*nu)/(144*nu^3*lambda^(2/nu))
    # D <- 1+(nu-1)*(A+B)
    # mu <- rep(max(lambda^(1/nu)-(nu-1)/(2*nu)+1/D*((nu-1)*(-A/nu+2*B/nu))),
    # length(lambda))
    # mu_error <- which(is.nan(mu)>0 | mu< 0)
    mu <- comp_means(lambda, nu, summax = 500)
    if (missing(summax)) {
      summax <- ceiling(max(c(mu + 20 * sqrt(mu / nu), 100)))
      cat("As you do not specify mu nor summax, summax will be calculated based on\n")
      cat("mu which is calcualted by truncated sum at 500.\n")
      cat("If you believe the mean of the distribution is somewhat close to 500,\n")
      cat("you may want to do some experiment with the comp_means() and\n")
      cat("specify summax instead to improve the accuracy.\n")
    }
  }
  # at a vector of yvalues
  pmf <- rep(0, length(x))
  for (i in 1:length(x)) {
    if ((mu[i] == 0 || lambda[i] == 0) && x[i] == 0) {
      pmf[i] <- 0 # log(1), 1 as the distribution is degenerated at 0
    } else if (mu[i] < 0 | lambda[i] < 0 | nu[i] <= 0) {
      pmf[i] <- NaN
      warn <- TRUE
    } else {
      if (!is.wholenumber(x[i])) {
        warning(paste("non-integer x =", x[i]))
        pmf[i] <- -Inf # log(0)
      } else {
        if (x[i] < 0) {
          pmf[i] <- -Inf
        } else { # log(0)
          # pmf <- log(density)
          pmf[i] <- x[i] * log(lambda[i]) - (nu[i] * lfactorial(x[i])) -
            logZ(log(lambda[i]), nu[i], summax)
        }
      }
    }
  }
  if (!log.p) {
    pmf <- exp(pmf)
  }
  if (warn) {
    warning("NaN(s) produced")
  }
  return(pmf)
}

rcomp <- function(n, mu, nu = 1, lambda, lambdalb = 1e-10,
                  lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6,
                  summax) {
  # generates random deviates of CMP variables with mean mu and dispersion nu
  # test to see at least one of mu and lambda is missing
  # mu, nu, lambda are recycled to give vectors length n
  # lambdaub will be scaled down/up if there is
  # over-/under-dispersion so that the correct lambda can be found
  if (length(n) > 1) {
    n <- length(n)
  }
  if (missing(mu) && missing(lambda)) {
    stop('argument "mu" is missing, with no default')
  }
  if (!missing(mu) && !missing(lambda)) {
    stop("specify 'mu' or 'lambda' but not both")
  }
  not.miss.mu <- !missing(mu)
  if (missing(mu)) {
    mu <- Inf
  }
  if (missing(lambda)) {
    lambda <- Inf
  }
  if (n < max(length(mu), length(nu), length(lambda))) {
    stop("unused argument in mu or nu or lambda")
  }
  df <- CBIND(x = rep(0, n), mu = mu, nu = nu, lambda = lambda)
  x <- df[, 1]
  mu <- df[, 2]
  nu <- df[, 3]
  lambda <- df[, 4]
  unif <- runif(n)
  warn <- FALSE
  if (not.miss.mu) {
    if (missing(summax)) {
      summax <- ceiling(max(c(mu + 20 * sqrt(mu / nu), 100)))
    }
    mu.ok.ind <- which(mu > 0)
    mu.err.ind <- which(mu <= 0)
    if (length(mu.err.ind) > 0) {
      lambda[mu.err.ind] <- mu[mu.err.ind]
    }
    if (length(mu.ok.ind) > 0) {
      lambda.ok <- comp_lambdas(mu[mu.ok.ind], nu[mu.ok.ind],
                                lambdalb = lambdalb, lambdaub = lambdaub,
                                # lambdaub = min(lambdaub,2*max(lambdaold))
                                maxlambdaiter = maxlambdaiter, tol = tol,
                                summax = summax
      )
      lambda[mu.ok.ind] <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
    }
  } else {
    # A <- (8*nu^2+12*nu+3)/(96*nu^2*lambda^(1/nu))
    # B <- (1+6*nu)/(144*nu^3*lambda^(2/nu))
    # D <- 1+(nu-1)*(A+B)
    # mu <- rep(max(lambda^(1/nu)-(nu-1)/(2*nu)+1/D*((nu-1)*(-A/nu+2*B/nu))),
    # length(lambda))
    # mu_error <- which(is.nan(mu)>0 | mu< 0)
    if (missing(summax)) {
      mu <- comp_means(lambda, nu, summax = 500)
      summax <- ceiling(max(c(mu + 20 * sqrt(mu / nu), 100)))
      cat("As you do not specify mu nor summax, summax will be calculated based on\n")
      cat("mu which is calcualted by truncated sum at 500.\n")
      cat("If you believe the mean of the distribution is somewhat close to 500,\n")
      cat("you may want to do some experiment with the comp_means() and\n")
      cat("specify summax instead to improve the accuracy.\n")
    }
  }
  for (i in 1:n) {
    if (mu[i] == 0 | lambda[i] == 0) {
      x[i] <- 0
    } else if (mu[i] < 0 | lambda[i] < 0 | nu[i] <= 0) {
      x[i] <- NA
      warn <- TRUE
    } else {
      y <- 0
      dc <- dcomp(0:summax, nu = nu[i], lambda = lambda[i], summax = summax)
      py <- dc[y + 1]
      while (py <= unif[i]) {
        y <- y + 1
        py <- py + dc[y + 1]
      }
      x[i] <- y
    }
  }
  if (warn) {
    warning("NAs produced")
  }
  return(x)
}

comp_lambdas <- function(mu, nu, lambdalb = 1e-10, lambdaub = 1000,
                         maxlambdaiter = 1e3, tol = 1e-6, lambdaint = 1, summax = 100) {
  df <- CBIND(mu = mu, nu = nu, lambda = lambdaint, lb = lambdalb, ub = lambdaub)
  mu <- df[, 1]
  nu <- df[, 2]
  lambda <- df[, 3]
  lambdalb <- df[, 4]
  lambdaub <- df[, 5]
  lambda.ok <-
    comp_lambdas_fixed_ub(mu, nu,
                          lambdalb = lambdalb, lambdaub = lambdaub,
                          maxlambdaiter = maxlambdaiter, tol = tol,
                          lambdaint = lambda, summax = summax
    )
  lambda <- lambda.ok$lambda
  lambdaub <- lambda.ok$lambdaub
  lambdaub.err.ind <- which(is.nan(lambda))
  sub_iter1 <- 1
  while (length(lambdaub.err.ind) > 0 && sub_iter1 <= 100) {
    lambdaub[lambdaub.err.ind] <- 0.5 * lambdaub[lambdaub.err.ind]
    lambda.ok <-
      comp_lambdas_fixed_ub(mu[lambdaub.err.ind], nu[lambdaub.err.ind],
                            lambdalb = lambdalb[lambdaub.err.ind],
                            lambdaub = lambdaub[lambdaub.err.ind],
                            maxlambdaiter = maxlambdaiter, tol = tol,
                            lambdaint = lambda[lambdaub.err.ind], summax = summax
      )
    lambda[lambdaub.err.ind] <- lambda.ok$lambda
    lambdaub[lambdaub.err.ind] <- lambda.ok$lambdaub
    sub_iter1 <- sub_iter1 + 1
    lambdaub.err.ind <- which(is.nan(lambda))
  }
  lambdaub.err.ind <- which(lambda / lambdaub >= 1 - tol)
  sub_iter1 <- 1
  while (length(lambdaub.err.ind) > 0 && sub_iter1 <= 100) {
    lambdaub[lambdaub.err.ind] <- 2^(sub_iter1) * lambdaub[lambdaub.err.ind]
    lambda.ok <-
      comp_lambdas_fixed_ub(mu[lambdaub.err.ind], nu[lambdaub.err.ind],
                            lambdalb = lambdalb[lambdaub.err.ind],
                            lambdaub = lambdaub[lambdaub.err.ind],
                            maxlambdaiter = maxlambdaiter, tol = tol,
                            lambdaint = lambda[lambdaub.err.ind],
                            summax = summax
      )
    lambda[lambdaub.err.ind] <- lambda.ok$lambda
    lambdaub[lambdaub.err.ind] <- lambda.ok$lambdaub
    sub_iter1 <- sub_iter1 + 1
    lambdaub.err.ind <- which(lambda / lambdaub >= 1 - tol)
  }
  out <- list()
  out$lambda <- lambda
  out$lambdaub <- max(lambdaub)
  return(out)
}

#' @rdname comp_lambdas
#' @export
comp_lambdas_fixed_ub <- function(mu, nu, lambdalb = 1e-10, lambdaub = 1000,
                                  maxlambdaiter = 1e3, tol = 1e-6, lambdaint = 1,
                                  summax = 100) {
  # df <- CBIND(mu=mu, nu=nu, lambda = lambdaint, lb = lambdalb, ub = lambdaub)
  # mu <- df[,1]
  # nu <- df[,2]
  # lambda <- df[,3]
  # lb <- df[,4]
  # ub <- df[,5]
  lambda <- lambdaint
  lb <- lambdalb
  ub <- lambdaub
  iter <- 1
  log.Z <- logZ(log(lambda), nu, summax = summax)
  mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
  not.converge.ind <- which(abs(mean1 - mu) > tol)
  while (length(not.converge.ind) > 0 && iter < 200) {
    still.above.target.ind <- which((mean1[not.converge.ind]
                                     > mu[not.converge.ind]))
    still.below.target.ind <- which(mean1[not.converge.ind] < mu[not.converge.ind])
    lb[not.converge.ind[still.below.target.ind]] <-
      lambda[not.converge.ind[still.below.target.ind]]
    ub[not.converge.ind[still.above.target.ind]] <-
      lambda[not.converge.ind[still.above.target.ind]]
    lambda <- (lb + ub) / 2
    log.Z <- logZ(log(lambda), nu, summax = summax)
    mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
    while (sum(mean1 == 0) > 0) {
      ub[not.converge.ind[mean1 == 0]] <- ub[not.converge.ind[mean1 == 0]] / 2
      lambdaub <- lambdaub / 2
      lambda <- (lb + ub) / 2
      log.Z <- logZ(log(lambda), nu, summax = summax)
      mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
    }
    not.converge.ind <- which((1 - (((abs(mean1 - mu) <= tol) + (lambda == lb) + (lambda == ub)
                                     + (ub == lb)) >= 1)) == 1)
    iter <- iter + 1
  }
  while (length(not.converge.ind) > 0 && iter < maxlambdaiter) {
    # basically the content of comp_variances without recalculating Z and mean1
    term <- matrix(0, nrow = length(mu), ncol = summax)
    for (y in 1:summax) {
      term[, y] <- exp(2 * log(y - 1) + (y - 1) * log(lambda) - nu * lgamma(y) - log.Z)
    }
    var1 <- apply(term, 1, sum) - mean1^2
    ## newton raphson update
    newtonsteps <- -lambda[not.converge.ind] * mean1[not.converge.ind] /
      (var1[not.converge.ind])^2 * (log(mean1[not.converge.ind]) - log(mu[not.converge.ind]))
    
    
    lambda.new <- lambda[not.converge.ind] + newtonsteps
    ## if newton raphson steps out of bound, use bisection method
    out.of.bound.ind <- which((lambda.new < lb[not.converge.ind])
                              + (lambda.new > ub[not.converge.ind]) == 1)
    if (length(out.of.bound.ind > 0)) {
      lambda.new[out.of.bound.ind] <-
        (lb[not.converge.ind[out.of.bound.ind]] + ub[not.converge.ind[out.of.bound.ind]]) / 2
      # any out of bound updates are replaced with mid-point of ub and lb
    }
    lambda[not.converge.ind] <- lambda.new
    log.Z <- logZ(log(lambda), nu, summax)
    term <- matrix(0, nrow = length(mu), ncol = summax)
    for (y in 1:summax) {
      term[, y] <- exp(log(y - 1) + (y - 1) * log(lambda) - nu * lgamma(y) - log.Z)
    }
    mean1 <- apply(term, 1, sum)
    if (length(out.of.bound.ind) > 0) {
      still.above.target.ind <- which(mean1[not.converge.ind[out.of.bound.ind]]
                                      > mu[not.converge.ind[out.of.bound.ind]])
      still.below.target.ind <- which(mean1[out.of.bound.ind] < mu[out.of.bound.ind])
      if (length(still.below.target.ind) > 0) {
        lb[not.converge.ind[out.of.bound.ind[still.below.target.ind]]] <-
          lambda[not.converge.ind[out.of.bound.ind[still.below.target.ind]]]
      }
      if (length(still.above.target.ind) > 0) {
        ub[not.converge.ind[out.of.bound.ind[still.above.target.ind]]] <-
          lambda[not.converge.ind[out.of.bound.ind[still.above.target.ind]]]
        # any out of bound updates are replaced with mid-point of ub and lb
      }
    }
    not.converge.ind <- which((1 - (((abs(mean1 - mu) <= tol) + (lambda == lb) + (lambda == ub)
                                     + (ub == lb)) >= 1)) == 1)
    iter <- iter + 1
  }
  out <- list()
  out$lambda <- lambda
  out$lambdaub <- lambdaub
  return(out)
}

comp_means <- function(lambda, nu, log.Z, summax = 100) {
  # approximates mean by truncation of COMP distributions
  # lambda, nu, mu.bd are recycled to match the length of each other.
  if (missing(log.Z)) {
    df <- CBIND(lambda = lambda, nu = nu)
    log.Z <- logZ(log(df[, 1]), df[, 2], summax)
  }
  df <- CBIND(lambda = lambda, nu = nu, log.Z = log.Z)
  lambda <- df[, 1]
  nu <- df[, 2]
  log.Z <- df[, 3]
  term <- matrix(0, nrow = length(lambda), ncol = summax)
  for (y in 1:summax) {
    term[, y] <- exp(log(y - 1) + (y - 1) * log(lambda) - nu * lgamma(y) - log.Z)
  }
  mean1 <- apply(term, 1, sum)
  return(mean1)
}