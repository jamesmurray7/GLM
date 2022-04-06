#' #######   mix/_03/
#' Simulates data under a Mixture of continuous, binary and count, which are linked to a 
#' survival sub-model by the vector of random effects, which are assumed to be Multivariate normal.
#' ---
#' This directory an implementation of Rustand et al. Simulation (ArXiv 220306256).
#' #######

# Establish 'true' parameters
# beta
true.beta <- rbind(c(0.20, -0.10, 0.10, -0.20),    # Gaussian
                   c(1.00, -1.00, 1.00, -1.00),    # Binomial
                   c(3.00, -0.10, 0.10, -0.20))    # Count (Poisson only (for now))
# Covariance matrix D
D <- diag(x = c(0.16, 0.09, 0.25, 0.25, 0.16))
D[lower.tri(D)] <- c(0.03, 0.00, 0.02, 0.04, -0.06, 0.03, 0.00, 0.05, 0.04, 0.08)
D[upper.tri(D)] <- t(D)[upper.tri(D)]
true.D <- D;rm(D)
# Survival
true.gamma <- c(0.50, 0.30, -0.20)

simData <- function(n = 250, ntms = 10, beta = true.beta, D = true.D, var.e = 0.4^2){
  nK <- nrow(beta)
  q <- ncol(D)
  if(any(eigen(D)$val < 0) | det(D) <= 0) stop('D not positive semi-definite :(')
  b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
  
  df <- data.frame(
    id = rep(1:n, each = ntms),
    time = rep(0:(ntms-1), n),
    cont = rep(rnorm(n, mean = 1, sd = 0.5), each = ntms),
    bin = rep(rbinom(n, 1, 0.5), each = ntms)
  ) 
  
  X <- model.matrix(~time+cont+bin, df)
  Z <- model.matrix(~time, df)
  eta1 <- X %*% beta[1,] + rowSums(Z * b[df$id, 1:2])
  eta2 <- X %*% beta[2,] + rowSums(Z[,1,drop=F] * b[df$id, 3, drop=F])
  eta3 <- X %*% beta[3,] + rowSums(Z * b[df$id, 4:5])
  
  df$Y.1 <- rnorm(n * ntms, eta1, sd = sqrt(var.e))
  df$Y.2 <- rbinom(n * ntms, 1, plogis(eta2))
  df$Y.3 <- rpois(n * ntms, lambda = exp(eta3))
  df
}

simData_joint <- function(n = 250, ntms = 10, beta = true.beta, var.e = 0.4^2,
                          D = true.D, gamma = true.gamma, theta = c(-4, 0.2),
                          cens.rate = exp(-3.5), fup = 5){
  nK <- nrow(beta)
  q <- ncol(D)
  if(any(eigen(D)$val < 0) | det(D) <= 0) stop('D not positive semi-definite :(')
  
  #' Necessary parameters & data generation ----
  # time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 
  # time variable and truncation time tau.
  tau <- fup + 0.1
  time <- seq(0, fup, length.out = ntms)
  cont <- rnorm(n, mean = 1, sd = 0.5); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  #' Design matrices ----
  X <- model.matrix(~ time + cont + bin, data = df)
  Z <- model.matrix(~ 1 + time, data = df)
  
  #' Linear predictors & response generation ----
  b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
  eta1 <- X %*% beta[1, ] + rowSums(Z * b[df$id, 1:2])
  eta2 <- X %*% beta[2, ] + rowSums(Z[, 1, drop = F] * b[df$id, 3, drop = F])
  eta3 <- X %*% beta[3, ] + rowSums(Z * b[df$id, 4:5])
  
  df$Y.1 <- c(eta1 + rnorm(n * ntms, sd = sqrt(var.e)))
  df$Y.2 <- rbinom(n * ntms, 1, plogis(eta2))
  df$Y.3 <- rpois(n * ntms, lambda = exp(eta3))
  
  if(any(is.na(df$Y.3))){
    message('NANs for\n')
    print(eta3[is.na(df$Y.3)])
    message('Which gave\n')
    print(exp(eta3)[is.na(df$Y.3)])
    message('end\n')
  }

  #' Survival ----
  theta0 <- theta[1]; theta1 <- theta[2]
  # Keta <- cbind(cont, bin) %*% surv.eta
  U <- runif(n); 
  b0 <- b[, c(1, 3, 4), drop = F]; b1 <- b[, c(2, 5), drop = F]
  # Generate survival times (Austin et al 2012)
  denom <- theta1 + b1 %*% gamma[c(1,3)]  # only slopes of 1 (gaussian) and 3 (count) affect the hazard.
  rhs <- (theta1 + b1 %*% gamma[c(1,3)]) * log(U)/(exp(theta0 + b0 %*% gamma)) # intercepts of all 3 affect hazard.
  t <- suppressWarnings(log(1-rhs)/denom)
  t[is.nan(t)] <- tau
  
  # Collect survival times, and generate cenors times.
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  #' Output Dataframes ----
  surv.data <- data.frame(id = 1:n, survtime, status)
  long.data <- df
  
  out.data <- dplyr::left_join(df, surv.data, by = 'id')
  out.data <- out.data[out.data$time < out.data$survtime, ]
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  list(data =  out.data, 
       surv.data =  dplyr::distinct(out.data, id, survtime, status))
}
