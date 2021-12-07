#' #######
#' Simulates data under a Poisson sub-model which is then linked with a 
#' survival sub-model by its random effects, which are assumed to be Gaussian.
#' #######

simData <- function(n = 250, ntms = 10, beta = c(1, 0.10, 0.33, -0.50),
                    D = matrix(c(0.5, 0, 0, 0.1), 2, 2), theta = 1.5){
  b <- MASS::mvrnorm(n, mu=c(0, 0), Sigma = D)
  
  df <- data.frame(
    id = rep(1:n, each = ntms),
    time = rep(0:(ntms-1), n),
    cont = rep(rnorm(n), each = ntms),
    bin = rep(rbinom(n, 1, 0.5), each = ntms)
  ) 
  
  X <- model.matrix(~time+cont+bin, df)
  Z <- model.matrix(~time, df)
  eta <- X %*% beta + rowSums(Z * b[df$id, ])
  
  df$Y <- rpois(n * ntms, lambda = exp(eta))
  df
}

simData_joint <- function(n = 250, ntms = 10, beta = c(1, 0.10, 0.33, -0.50),
                          D = matrix(c(0.5, 0, 0, 0.1), 2, 2), thetaDisp = 1.5,
                          gamma = 0.5, surv.eta = c(0.05, -0.30), theta = c(-4, 0.2),
                          cens.rate = exp(-3.5)){
  #' Necessary parameters & data generation ----
  time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 # time variable and truncation time tau.
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  #' Design matrices ----
  X <- model.matrix(~ time + cont + bin, data = df)
  Z <- model.matrix(~ 1 + time, data = df)
  
  #' Linear predictor & response generation ----
  b <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = D)
  eta <- X %*% beta + rowSums(Z * b[df$id, ])
  df$Y <- rpois(n * ntms, lambda = exp(eta))
  
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
       surv.data =  dplyr::distinct(out.data, id, survtime, status, cont, bin))
}
