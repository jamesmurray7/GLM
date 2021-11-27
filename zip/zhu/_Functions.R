#' #####
#' Helper functions for a (+ zero-inflated) Poisson 
#' for mimicking of of Zhu et al., (2018) work.
#' ----
#' For time being these just focussed on ZIP part, survival added in due course!
#' #####

# Data simulation ---------------------------------------------------------
# 1. ZIP. Assumes RE are random intercept only.
simData_zip <- function(n, ntms, beta, alpha, D = diag(2)){
  time <- 0:(ntms - 1)
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  # Design matrices
  # Design matrices for the fixed and random effects in the non-zero inflated part.
  X <- model.matrix(~ time + cont + bin, data = df)
  Z <- model.matrix(~ 1, data = df)
  # Design matrices for the fixed and random effects in the zero inflated (ZI) part.
  Xzi <- model.matrix(~ time, data = df)
  Zzi <- model.matrix(~ 1, data = df)
  
  # Random Effects,
  b <- MASS::mvrnorm(n, rep(0, 2), Sigma = D);
  # bzi <- MASS::mvrnorm(n, rep(0, 2), Sigma = D[[2]])
  
  # Linear predictors
  eta <- X %*% beta + rowSums(Z * b[df$id, 1, drop=F])
  etazi <- Xzi %*% alpha + rowSums(Zzi * b[df$id, 2, drop=F])
  
  # Outcome y
  Y <- rpois(n * ntms, lambda = exp(eta))
  Y[as.logical(rbinom(n * ntms, 1, prob = plogis(etazi)))] <- 0
  df$Y <- Y
  
  message(sum(df$Y == 0)/length(df$Y)*100,'% zeroes')
  
  df
}

# Simulate under a joint modelling framework (specifically, the second simulation in Zhu et al. (2018)).
zhuD <- function() return(matrix(c(1.2, -0.5 * sqrt(1.2) * sqrt(0.6), -0.5 * sqrt(1.2) * sqrt(0.6), 0.6), 2, 2))
zhubeta <- function() return(c(0.5, -0.2, 0.4))
zhualpha <- function() return(c(0.8, 0.2, -0.3))
simData_zip_joint <- function(n, ntms = 6, beta = zhubeta(), alpha =zhualpha(), D = zhuD(),
                              theta = c(-4, 0.2), surv.eta = c(-1.5), 
                              gamma = c(-0.6, 0.4), cens.rate = exp(-3.5)){
  
  # Necessary parameters
  id <- 1:n
  time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 # time variable and truncation time tau.
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
               #   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  # Design matrices
  # Design matrices for the fixed and random effects in the non-zero inflated part.
  X <- model.matrix(~ time + bin, data = df)
  Z <- model.matrix(~ 1, data = df)
  # Design matrices for the fixed and random effects in the zero inflated (ZI) part.
  Xzi <- model.matrix(~ time + bin, data = df)
  Zzi <- model.matrix(~ 1, data = df)
  
  # Random Effects,
  b <- MASS::mvrnorm(n, rep(0, 2), Sigma = D);
  # bzi <- MASS::mvrnorm(n, rep(0, 2), Sigma = D[[2]])
  
  # Linear predictors
  eta <- X %*% beta + rowSums(Z * b[df$id, 1, drop=F])
  etazi <- Xzi %*% alpha + rowSums(Zzi * b[df$id, 2, drop=F])
  
  # Outcome y
  Y <- rpois(n * ntms, lambda = exp(eta))
  Y[as.logical(rbinom(n * ntms, 1, prob = plogis(etazi)))] <- 0
  df$Y <- Y; message(sum(df$Y == 0)/length(df$Y)*100,'% zeroes')
  
  # Simulating survival times
  theta0 <- theta[1]; theta1 <- theta[2]
  K <- as.matrix(bin)
  Keta <- K %*% surv.eta
  U <- runif(n)
  
  # denom <- theta1 + b1 %*% gamma                     # these for if I put intercept + slope back in!
  # rhs <- (theta1 + b1 %*% gamma) * log(U)/(exp(theta0 + Keta + b0 %*% gamma))   
  denom <- theta1
  rhs <- theta1 * log(U)/(exp(theta0 + Keta + b %*% gamma))
  
  t <- suppressWarnings(log((1 - rhs))/denom)
  t[is.nan(t)] <- tau
  
  # Collect survival times, and generate cenors times.
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  # Output Dataframes
  surv.data <- data.frame(id, survtime, status)
  long.data <- df
  
  out.data <- dplyr::left_join(df, surv.data, by = 'id')
  out.data <- out.data[out.data$time < out.data$survtime, ]
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  list(data =  out.data, 
       surv.data =  dplyr::distinct(out.data, id, survtime, status, bin))
  
}

