#' #######
#' Simulates data under a the following joint modelling specification:
#'   Longitudinal process: where marginal distribution of Y given random effects is
#'                             assumed to belong to a member of the exponential family
#'   Survival process: Survival sub-model which is then linked with the above longitudinal process by some
#'                        shared random effects, which are assumed to be Gaussian.
#'                        
#'  Family: Can be "gaussian"; "binomial"; "poisson"; "negative.binomial"
#' #######

simData <- function(n = 250, ntms = 10, fup = 5, family = 'gaussian', disp = NULL, 
                    beta = c(1, 0.10, 0.33, -0.50), D = matrix(c(0.5, 0, 0, 0.1), 2, 2), var.e = 0.16,
                    gamma = 0.5, zeta = c(0.05, -0.30), theta = c(-4, 0.2), cens.rate = exp(-3.5)){
  
  #' Checks --------------
  # Check family is valid option
  if("function"%in%class(family)) family <- family()$family
  if(!family%in%c("gaussian", "binomial", "poisson", "negative.binomial")) stop('Family must be one of "gaussian", "binomial", "poisson" or "negative.binomial"')
  # Check dispersion supplied if neg. binom.
  if(family == 'negative.binomial' & is.null(disp)) stop('Must define disp (scalar) if negative.binomial model chosen')
  family <- match.arg(family, c('gaussian', 'binomial', 'poisson', 'negative.binomial'), several.ok = F)
  # Check covariance matrix D is positive semi-definite.
  if(any(eigen(D)$value <= 0) || det(D) <= 0) stop('D must be positive semi-definite')
  
  #' Necessary parameters & data generation ----
  time <- seq(0, fup, length.out = ntms); tau <- fup + 0.1
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
  
  switch(family, 
         gaussian = Y <- rnorm(n * ntms, mean = eta, sd = sqrt(var.e)),
         binomial = Y <- rbinom(n * ntms, 1, plogis(eta)),
         poisson = Y <- rpois(n * ntms, exp(eta)),
         negative.binomial = Y <- MASS::rnegbin(n * ntms, mu = exp(eta), theta = disp)
         )
  
  df$Y <- Y
  
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


