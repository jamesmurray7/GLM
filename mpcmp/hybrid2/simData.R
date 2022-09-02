#' #######
#' Simulates data under a multivariate negative binomial model which is then linked with a 
#' survival sub-model by its random effects, which are assumed to be Gaussian.
#' ----
#' Bespoke functions to find pmf of, and then simulate random deviates from a
#' mean-parameterised specification of the Conway-Maxwell Poisson Distn.
#' #######
# rm(list=ls())
# old.ls <- ls()
source('../mpcmp-copy/mpcmp-copy.R')
# to.remove <- setdiff(ls(), old.ls)

simData_joint <- function(n = 250, ntms = 10, summax = 100,  fup = 5,
                          beta = c(0.0, -0.1, 0.05, -0.1),           
                          delta = c(0.05, -0.1),     
                          D = matrix(c(0.16, 0, 0, 0.04), 2, 2), 
                          gamma = 0.3, zeta = c(0.05, -0.20), theta = c(-3, 0.25),
                          cens.rate = exp(-3.5)){
  #' Necessary parameters & data generation ----
  time <- seq(0, fup, length.out = ntms); tau <- fup + 0.1
  # time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 # time variable and truncation time tau.
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  b <- MASS::mvrnorm(n, mu=c(0, 0), Sigma = D)
  #' Survival ----
  theta0 <- theta[1]; theta1 <- theta[2]
  Keta <- cbind(cont, bin) %*% zeta
  U <- runif(n); 
  b0 <- b[, 1, drop = F]; b1 <- b[, 2, drop = F]
  # Generate survival times (Austin et al 2012)
  denom <- theta1 + b1 %*% gamma  
  rhs <- (theta1 + b1 %*% gamma) * log(U)/(exp(theta0 + Keta + b0 %*% gamma))
  t <- suppressWarnings(log(1-rhs)/denom)
  t[is.nan(t)] <- tau
  
  # Collect survival times, and generate censor times.
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  #' Join onto Longtudinal part and truncate ----
  surv.data <- data.frame(id = 1:n, survtime, status)
  
  df <- dplyr::left_join(df, surv.data, by = 'id')
  df <- df[df$time < df$survtime, ]
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  # Design matrices 
  X <- model.matrix(~ time + cont + bin, data = df)
  G <- Z <- model.matrix(~ time, data = df)
  
  #' Simulate CMP response //
  eta <- X %*% beta + rowSums(Z * b[df$id, ]) # Linear predictor
  nu <- exp(G %*% delta)                      # Dispersion (time-varying, set delta[2]=0 for intercept only.)
  
  mu <- exp(eta)                              # Mean parameterisation
  
  Y <- vector('list', n)
  pb <- utils::txtProgressBar(max = n, style = 3)
  for(i in 1:n){
    ind <- which(df$id==i)
    Y[[i]] <- rcomp(n = length(ind), mu = mu[ind], nu = nu[ind], summax = 100)
    utils::setTxtProgressBar(pb, i)
  }
  df$Y <- do.call(c, Y)
  cat('\n')
  list(data = df, 
       surv.data =  dplyr::distinct(df, id, survtime, status, cont, bin))
}


# Quicker version ---------------------------------------------------------
.rcomp <- function(lambda, nu, summax){
  x <- vector('numeric', length(lambda))
  U <- runif(length(lambda))
  warn <- F
  for(i in 1:length(lambda)){
    if (lambda[i] == 0) {
      x[i] <- 0
    } else if (lambda[i] < 0 | nu[i] <= 0) {
      x[i] <- NA
      warn <- TRUE
    } else {
      y <- 0
      dc <- cmp_pmf_scalar(0:(summax+1), lambda[i], nu[i], summax)
      py <- dc[y + 1]
      while (py <= U[i]) {
        y <- y + 1
        py <- py + dc[y + 1]; 
      }
      x[i] <- y
    }
  }
  if (warn) {
    warning("NAs produced")
  }
  return(x)
}


simData_joint2 <- function(n = 250, ntms = 10, summax = 100,  fup = 5,
                           beta = c(0.0, -0.1, 0.05, -0.1),           
                           delta = c(0.05, -0.1),     
                           D = matrix(c(0.16, 0, 0, 0.04), 2, 2), 
                           gamma = 0.3, zeta = c(0.05, -0.20), theta = c(-3, 0.25),
                           cens.rate = exp(-3.5)){
  #' Necessary parameters & data generation ----
  time <- seq(0, fup, length.out = ntms); tau <- fup + 0.1
  # time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 # time variable and truncation time tau.
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  b <- MASS::mvrnorm(n, mu=c(0, 0), Sigma = D)
  #' Survival ----
  theta0 <- theta[1]; theta1 <- theta[2]
  Keta <- cbind(cont, bin) %*% zeta
  U <- runif(n); 
  b0 <- b[, 1, drop = F]; b1 <- b[, 2, drop = F]
  # Generate survival times (Austin et al 2012)
  denom <- theta1 + b1 %*% gamma  
  rhs <- (theta1 + b1 %*% gamma) * log(U)/(exp(theta0 + Keta + b0 %*% gamma))
  t <- suppressWarnings(log(1-rhs)/denom)
  t[is.nan(t)] <- tau
  
  # Collect survival times, and generate censor times.
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  #' Join onto Longtudinal part and truncate ----
  surv.data <- data.frame(id = 1:n, survtime, status)
  
  df <- dplyr::left_join(df, surv.data, by = 'id')
  df <- df[df$time < df$survtime, ]
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  # Design matrices 
  X <- model.matrix(~ time + cont + bin, data = df)
  G <- Z <- model.matrix(~ time, data = df)
  
  #' Simulate CMP response //
  eta <- X %*% beta + rowSums(Z * b[df$id, ]) # Linear predictor
  nu <- exp(G %*% delta)                      # Dispersion (time-varying, set delta[2]=0 for intercept only.)
  
  mu <- exp(eta)                              # Mean parameterisation
  
  Y <- vector('list', n)
  pb <- utils::txtProgressBar(max = n, style = 3)
  for(i in 1:n){
    ind <- which(df$id==i)
    lam <- lambda_appx(mu[ind], nu[ind], summax)
    Y[[i]] <- .rcomp(lam, nu[ind], summax)
    utils::setTxtProgressBar(pb, i)
  }
  df$Y <- do.call(c, Y)
  cat('\n')
  list(data = df, 
       surv.data =  dplyr::distinct(df, id, survtime, status, cont, bin))
}