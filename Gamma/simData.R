#' #######
#' Simulates data under a Gamma regression model which is then linked with a 
#' survival sub-model by its random effects, which are assumed to be Gaussian.
#' ----
#' `Shape` is provided as a double, scale is then h^-1(\eta)/`shape`.
#' A log-link is used on the linear predictor \eta to simulate data.
#' #######

simData_joint <- function(n = 250, ntms = 10, fup = 3, 
                          beta = c(2, -0.1, 0.1, -0.2),               # (int, time, cont, bin)
                          shape = c(2),                               # (double)
                          D = matrix(c(0.25, 0, 0, 0.00), 2, 2), 
                          gamma = 0.6, zeta = c(0.00, -0.20), theta = c(-2, 0.10),
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
  Z <- model.matrix(~ time, data = df)
  
  #' Simulate CMP response //
  eta <- X %*% beta + rowSums(Z * b[df$id, ]) # Linear predictor
  mu <- exp(eta)                              # Mean(log link!)
  
  df$Y <- rgamma(length(mu), shape = shape, scale = mu/shape)
  
  cat('\n')
  list(data = df, 
       surv.data =  dplyr::distinct(df, id, survtime, status, cont, bin),
       true.RE = b)
}

# glmmTMB(Y ~ time + cont + bin + (1|id), data = data, family = Gamma(link = 'log'))