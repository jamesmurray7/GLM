#' #######
#' Simulates data under a Mixture of continuous, binary and count, which are linked to a 
#' survival sub-model by the vector of random effects, which are assumed to be Multivariate normal.
#' #######

# Establish 'true' parameters
# beta
true.beta <- rbind(c(0.00, 1.00, 0.50, -0.25),    # Gaussian
                   c(1.00, 0.10, 0.33, -0.50),    # Binomial
                   c(0.00, 0.05, 0.01, 0.30))     # Count (poiss/negbinom)
# Covariance matrix D
true.D <- as.matrix(Matrix::bdiag(
  matrix(c(0.5^2, 0.1, 0.1, 0.25^2), 2, 2),
  matrix(c(0.6, -0.05, -0.05, 0.1), 2, 2),
  matrix(c(0.25, 0.01, 0.01, 0.05), 2, 2)
))
true.D[5:6, 1:2] <- true.D[1:2, 5:6] <- matrix(c(0.125, -0.050, -0.050, 0.200))
true.D[3:4, 1:2] <- true.D[1:2, 3:4] <- matrix(c(0.005, 0.000, 0.000, -0.005))
true.D <- as.matrix(Matrix::nearPD(true.D, keepDiag = T)$mat)
# Survival
true.gamma <- c(0.50, -0.25, 0.40)
true.eta <- c(0.05, -0.30)

simData <- function(n = 250, ntms = 10, beta = true.beta,
                    D = true.D, Disptheta = NULL){
  
  nK <- nrow(beta)
  q <- ncol(D)
  if(any(eigen(D)$val < 0) | det(D) <= 0) stop('D not positive semi-definite :(')
  if(nK != q/2) stop('nrow(beta) must = nrow(D)/2; assuming int-slope on each K response')
  b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
  
  df <- data.frame(
    id = rep(1:n, each = ntms),
    time = rep(0:(ntms-1), n),
    cont = rep(rnorm(n), each = ntms),
    bin = rep(rbinom(n, 1, 0.5), each = ntms)
  ) 
  
  X <- model.matrix(~time+cont+bin, df)
  Z <- model.matrix(~time, df)
  eta1 <- X %*% beta[1,] + rowSums(Z * b[df$id, 1:2])
  eta2 <- X %*% beta[2,] + rowSums(Z * b[df$id, 3:4])
  eta3 <- X %*% beta[3,] + rowSums(Z * b[df$id, 5:6])
  
  df$Y.1 <- eta1
  df$Y.2 <- rbinom(n * ntms, 1, plogis(eta2))
  if(is.null(Disptheta)){
    df$Y.3 <- rpois(n * ntms, lambda = exp(eta3))
  }else{
    df$Y.3 <- MASS::rnegbin(n * ntms, mu = exp(eta3), theta = rep(Disptheta, n * ntms))
  }
  df
}

simData_joint <- function(n = 250, ntms = 10, beta = true.beta, var.e = 0.25,
                          D = true.D, Disptheta = NULL,
                          gamma = true.gamma, surv.eta = true.eta, theta = c(-4, 0.2),
                          cens.rate = exp(-3.5)){
  nK <- nrow(beta)
  q <- ncol(D)
  if(any(eigen(D)$val < 0) | det(D) <= 0) stop('D not positive semi-definite :(')
  if(nK != q/2) stop('nrow(beta) must = nrow(D)/2; assuming int-slope on each K response')
  
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
  
  #' Linear predictors & response generation ----
  b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
  eta1 <- X %*% beta[1, ] + rowSums(Z * b[df$id, 1:2])
  eta2 <- X %*% beta[2, ] + rowSums(Z * b[df$id, 3:4])
  eta3 <- X %*% beta[3, ] + rowSums(Z * b[df$id, 5:6])
  
  df$Y.1 <- eta1 + rnorm(n * ntms, sd = sqrt(var.e))
  df$Y.2 <- rbinom(n * ntms, 1, plogis(eta2))
  if(is.null(Disptheta)){
    df$Y.3 <- rpois(n * ntms, lambda = exp(eta3))
  }else{
    df$Y.3 <- MASS::rnegbin(n * ntms, mu = exp(eta3), theta = rep(Disptheta, n * ntms))
  }
  
  #' Survival ----
  theta0 <- theta[1]; theta1 <- theta[2]
  Keta <- cbind(cont, bin) %*% surv.eta
  U <- runif(n); 
  b0 <- b[, c(1, 3, 5), drop = F]; b1 <- b[, c(2, 4, 6), drop = F]
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
