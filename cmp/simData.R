#' #######
#' Simulates data under a multivariate negative binomial model which is then linked with a 
#' survival sub-model by its random effects, which are assumed to be Gaussian.
#' ----
#' Bespoke functions to find pmf of and then simulate random deviates from a
#' mean-parameterised specification of the Conway-Maxwell Poisson Distn.
#' #######
sourceCpp('test.cpp')

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

simData <- function(n = 250, ntms = 10, summax = 100,
                    beta = c(1, 0.10, 0.33, -0.50),           # Coeffs FE lin. pred
                    delta = c(-0.6, -0.1),                    # Coeffs FE Dispersion
                    D = matrix(c(0.25, 0, 0, 0.04), 2, 2)){   # RE covariance matrix.
  
  N <- n * ntms
  b <- MASS::mvrnorm(n, mu=c(0, 0), Sigma = D)
  
  df <- data.frame(
    id = rep(1:n, each = ntms),
    time = rep(0:(ntms-1), n),
    cont = rep(rnorm(n), each = ntms),
    bin = rep(rbinom(n, 1, 0.5), each = ntms)
  ) 
  
  X <- model.matrix(~time+cont+bin, df)
  Z <- G <- model.matrix(~time, df)
  eta <- X %*% beta + rowSums(Z * b[df$id, ]) # Linear predictor
  nu <- exp(G %*% delta)
  
  mu <- exp(eta) # Define the mean parameterisation
  
  # Working out lambda, the rate parameter from mu
  # i. An approximation 
  loglambdas.appx <- suppressWarnings(
    nu * log(mu + (nu - 1) / (2 * nu))
  )
  lambdas.appx <- exp(loglambdas.appx)
  
  # ii. Find solutions to mean constraint (Huang (2017)) and clean so NaN/Infs/NAs not in output.
  lambdas <- sapply(1:nrow(eta), function(i){
    out <- tryCatch(uniroot(mu_lambdaZ_eq, interval = c(1e-6, 1e3), mu = exp(eta[i]), nu = nu[i],
                            summax = summax)$root,
                    error = function(e) NA)
    # If uniroot fails to find a root, set it as the approximation above
    if((is.na(out) | is.nan(out)) & (!is.nan(lambdas.appx[i]) & !is.na(lambdas.appx[i]))) out <- lambdas.appx[i]
    # And if this still NA/NaN/Inf, simply set as mean
    if(is.na(out)) out <- mu[i]
    out
  })
  
  # Print how many rate parameters simply used the mean.
  sprintf('%.2f%% values used mean', length(which(lambdas == mu))/N * 100)
  
  sf <- 1
  Y <- try(.rcomp(lambdas, nu, summax), silent = T)
  while('try-error' %in% class(Y)){
    sf <- sf * .90
    cat(sprintf("Attempting with new scale factor: %.1f\n", sf * summax))
    Y <- try(.rcomp(lambdas, nu, summax * sf), silent = T)
  }
  df$Y <- Y
  df
}

# fit <- glmmTMB::glmmTMB(Y~time+cont+bin+(1+time|id), data=test, family = glmmTMB::nbinom2, dispformula = ~1)

simData_joint <- function(n = 250, ntms = 10, summax = 100,
                          beta = c(1, 0.10, 0.33, -0.50),           # Coeffs FE lin. pred
                          delta = c(-0.6, -0.1),     
                          D = matrix(c(0.5^2, 0, 0, 0.2^2), 2, 2), 
                          gamma = 0.5, surv.eta = c(0.05, -0.30), theta = c(-4, 0.2),
                          cens.rate = exp(-3.5)){
  #' Necessary parameters & data generation ----
  time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 # time variable and truncation time tau.
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  b <- MASS::mvrnorm(n, mu=c(0, 0), Sigma = D)
  #' Survival ----
  theta0 <- theta[1]; theta1 <- theta[2]
  Keta <- cbind(cont, bin) %*% surv.eta
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
  
  # Design matrices 
  X <- model.matrix(~ time + cont + bin, data = df)
  G <- Z <- model.matrix(~ time, data = df)
  
  #' Simulate CMP response //
  eta <- X %*% beta + rowSums(Z * b[df$id, ]) # Linear predictor
  nu <- exp(G %*% delta)
  
  mu <- exp(eta) # Define the mean parameterisation
  
  # Working out lambda, the rate parameter from mu
  # i. An approximation 
  loglambdas.appx <- suppressWarnings(
    nu * log(mu + (nu - 1) / (2 * nu))
  )
  lambdas.appx <- exp(loglambdas.appx)
  
  # ii. Find solutions to mean constraint (Huang (2017)) and clean so NaN/Infs/NAs not in output.
  lambdas <- sapply(1:nrow(eta), function(i){
    out <- tryCatch(uniroot(mu_lambdaZ_eq, interval = c(1e-6, 1e3), mu = mu[i], nu = nu[i],
                            summax = summax)$root,
                    error = function(e) NA)
    # If uniroot fails to find a root, set it as the approximation above
    if((is.na(out) | is.nan(out)) & (!is.nan(lambdas.appx[i]) & !is.na(lambdas.appx[i]))) out <- lambdas.appx[i]
    # And if this still NA/NaN/Inf, simply set as mean
    if(is.na(out)) out <- mu[i]
    out
  })
  
  # Print how many rate parameters simply used the mean.
  sprintf('%.2f%% values used mean', length(which(lambdas == mu))/length(lambdas) * 100)
  
  sf <- 1
  Y <- try(.rcomp(lambdas, nu, summax), silent = T)
  while('try-error' %in% class(Y)){
    sf <- sf * .90
    cat(sprintf("Attempting with new scale factor: %.1f\n", sf * summax))
    Y <- try(.rcomp(lambdas, nu, summax * sf), silent = T)
  }
  df$Y <- Y
  df
  
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  list(data = df, 
       surv.data =  dplyr::distinct(df, id, survtime, status, cont, bin))
}

