#' #######
#' Simulates data under a multivariate Poisson model which is then linked with a 
#' survival sub-model by its random effects, which are assumed to be Gaussian.
#' Should simulate data under K = 1+ responses, which will ALL be poisson.
#' #######

simData <- function(n, ntms, beta, D, gamma, eta, 
                    theta = c(-6, 0.15), cens.rate = exp(-3.5)){
  # RE & general stuff & checks pre-simulation ----
  nK <- nrow(beta)
  q <- length(diag(D))
  if(nK^2 != q) stop('Only intercept and slope models fit, nK != q')
  if(!isSymmetric(D)) stop('D must be symmetric')
  if(any(eigen(D)$values < 0) || (det(D) <= 0)) stop("Covariance matrix must be positive semi-definite")
  
  b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
  # Split-out some helpful indices
  qk <- split(seq(q), cut(seq_along(seq(q)), nK, labels = F)) 
  if(nK > 1){
    int.slope.inds <- split(seq(nK * 2), rep(1:nK, 2))
  }else{
    int.slope.inds <- list(1, 2)
  }
  b0 <- b[, int.slope.inds[[1]]]; b1 <- b[, int.slope.inds[[2]]] # Pull-out intercept(s) and slope(s) for later use...
  
  # Necessary parameters
  id <- 1:n
  time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 # time variable and truncation time tau.
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
    
  # Simulate longitudinal outcome
  dfs <- list(); 
  for(i in 1:n){
    Yijk <- lambda_ijks <- matrix(NA, nr = length(time), nc = nK) # i sub; j time; k longit response
    for(k in 1:nK){
      lambda_ijks[, k] <- exp(cbind(1, time, cont[i], bin[i]) %*% beta[k, ] + cbind(1, time) %*% b[i, qk[[k]]])
      for(j in seq_along(lambda_ijks[, k])){
        Yijk[j, k] <- rpois(1, lambda_ijks[j, k])
      }
    }
    colnames(Yijk) <- paste0('Y.', 1:nK)
    dfs[[i]] <- data.frame(id = id[i], time = time, cont = cont[i], bin = bin[i], Yijk)
  }
  
  # Simulating survival times
  theta0 <- theta[1]; theta1 <- theta[2]
  K <- cbind(cont, bin)
  Keta <- K %*% eta
  U <- runif(n)
  
  denom <- theta1 + b1 %*% gamma
  rhs <- (theta1 + b1 %*% gamma) * log(U)/(exp(theta0 + Keta + b0 %*% gamma))
  t <- suppressWarnings(log((1 - rhs))/denom)

  t[is.nan(t)] <- tau
  
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  # Output Dataframes
  surv.data <- data.frame(id, survtime, status)
  long.data <- do.call(rbind, dfs)
  
  dat <- dplyr::left_join(long.data, surv.data, 'id')
  dat <- dat[dat$time <= dat$survtime, ]
  
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  dat
}