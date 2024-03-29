#' #######
#' Simulates data under a Poisson model which is then linked with a 
#' survival sub-model by its random effects, which are assumed to be Gaussian.
#' #######

n <- 100; id <- 1:n
D <- matrix(c(0.5, 0, 0, 0.1), 2, 2)
b <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = D)
time <- 0:5
cont <- rnorm(n)
bin <- rbinom(n, 1, .5)
beta <- c(3.25, 0.15, -0.15, 0.6)
dfs <- list(); lambda_is <- numeric(n)
for(i in 1:n){
  Yi <- numeric(length(time))
  lambda_i <- cbind(1, time, cont[i], bin[i]) %*% beta + cbind(1, time) %*% b[i, ]
  for(j in seq_along(lambda_i)){       # This probably very inefficient, but unsure of behaviour of rpois with a vector of lambdas...
    Yi[j] <- rpois(1, lambda_i[j])
  }
  dfs[[i]] <- data.frame(id = id[i], time = time, cont = cont[i], bin = bin[i], Y = Yi)
}

K <- cbind(cont, bin)
eta <- c(-0.2, .5)
# Simulating surviva; times
Keta <- K %*% eta
U <- runif(n)

theta0 <- -6; theta1 <- .15

denom <- theta1 * b[, 2, drop = F] %*% gamma
rhs <- ((theta1 + b[, 1, drop = F] %*% gamma) * log(U))/(exp(theta0 + Keta + b[, 1, drop = F] %*% gamma))
t <- suppressWarnings(log((1 - rhs)/denom))

t[is.nan(t)] <- 5.1

cens.time <- rexp(n, exp(-3))
surv.time <- pmin(t, cens.time)
surv.time[surv.time >= 5.1] <- 5.1

# Status flag
status <- rep(1, n)
is.censored <- cens.time < surv.time
status[which(surv.time == 5.1 | is.censored | surv.time == cens.time)] <- 0 # Failure flag

# Output Dataframes
surv.data <- data.frame(id, surv.time, status)
long.data <- do.call(rbind, dfs)

dat <- dplyr::left_join(long.data, surv.data, 'id')
dat <- dat[dat$time <= dat$surv.time, ]

# Make a function as we go ------------------------------------------------

# UNIVARIATE
simData <- function(n, ntms, beta, D, gamma, eta, 
                    theta = c(-6, 0.15), cens.rate = exp(-3.5)){
  # RE stuff
  q <- length(diag(D))
  if(!isSymmetric(D)) stop('D must be symmetric')
  if(any(eigen(D)$values < 0) || (det(D) <= 0)) stop("Covariance matrix must be positive semi-definite")
  b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
  b0 <- b[, 1, drop = F]; b1 <- b[, 2, drop = F] # Pull-out intercept and slope for later use...
  
  # Necessary parameters
  id <- 1:n
  time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 # time variable and truncation time tau.
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)
    
  # Simulate longitudinal outcome
  dfs <- list(); lambda_is <- numeric(n)
  for(i in 1:n){
    Yi <- numeric(length(time))
    lambda_i <- cbind(1, time, cont[i], bin[i]) %*% beta + cbind(1, time) %*% b[i, ]
    for(j in seq_along(lambda_i)){       # This probably very inefficient, but unsure of behaviour of rpois with a vector of lambdas...
      Yi[j] <- rpois(1, lambda_i[j])
    }
    dfs[[i]] <- data.frame(id = id[i], time = time, cont = cont[i], bin = bin[i], Y = Yi)
  }
  
  # Simulating survival times
  theta0 <- theta[1]; theta1 <- theta[2]
  K <- cbind(cont, bin)
  Keta <- K %*% eta
  U <- runif(n)
  
  denom <- theta1 * b1 %*% gamma
  rhs <- ((theta1 + b0 %*% gamma) * log(U))/(exp(theta0 + Keta + b0 %*% gamma))
  t <- suppressWarnings(log((1 - rhs)/denom))
  
  t[is.nan(t)] <- tau
  
  cens.time <- rexp(n, cens.rate)
  surv.time <- pmin(t, cens.time)
  surv.time[surv.time >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < surv.time
  status[which(surv.time == tau | is.censored | surv.time == cens.time)] <- 0 # Failure flag
  
  # Output Dataframes
  surv.data <- data.frame(id, surv.time, status)
  long.data <- do.call(rbind, dfs)
  
  dat <- dplyr::left_join(long.data, surv.data, 'id')
  dat <- dat[dat$time <= dat$surv.time, ]
  
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  dat
}

rm(list=ls())

simData(150, 10, c(15, -1, 0.5, 3), matrix(c(.5, 0, 0, .05), 2, 2), 0.5, c(0.05, 0.33))

