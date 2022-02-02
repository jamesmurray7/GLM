#' ####
#' Simulating data under a quadratic(+) random effects structure (to be fitted using splines)
#' ----
#' Arguments
#' n: Number of subjects
#' ntms: maximum profile length
#' order: Order of quadratic time structure
#' beta: Longitudinal coefficients (nK x 5); {intercept, time, time^2, cont, bin}
#' eta: Survival parameters
#' gamma: Association parameters (proportional assoc. only)
#' var.e: nK-length vector of variances
#' D: Covariance matrix for RE, of dim (nKx3) * (nKx3)
#' theta: Baseline hazard control: theta0 scale; theta1 shape.
#' cens.rate: Underlying hazard for censoring to occur.
#' dt: Number of grid steps to use.
#' ####
library(dplyr)

# Making time structure via separate function.
make.design.matrix <- function(timevec, order){
  times <- matrix(NA, nr = length(timevec), nc = order)
  for(i in seq(order)){
    times[,i] <- timevec^i
  }
  structure(cbind(1, times),
            dimnames = list(as.character(1:nrow(times)),
                            c('(Intercept)', paste0('time', 1:order)))
  )
}

poly.simData <- function(n = 150, ntms = 5, order = 2, beta = matrix(1, nr = 3, nc = 5),
                         eta = c(0.05, -0.3), gamma = c(0.50, -0.25, 0.40), var.e = rep(0.25, 3),
                         D = NULL, theta0 = -6, theta1 = 0.25, cens.rate = exp(-3), dt = 0.01){
  
  nK <- nrow(beta) # K
  q <- nK * (order + 1)   # total number of random effects (Intercept) + order.
  tau <- (ntms - 1) + 0.1 # Truncation time
  
  # Checks ----
  if(ncol(beta) != (1 + order + 2)){
    stop('Incorrect number of columns for beta: (Intercept), time1, ..., time^order, cont, bin needed')
  }
  if(!is.null(D)){
    if(ncol(D) != q){
      stop('D of insufficient dimension for the number of random effects: nrow(D) = ', nrow(D), ', q = ', q)
    }
    if(any(eigen(D)$values < 0) || (det(D) <= 0)){
      stop("Covariance matrix must be positive semi-definite")
    }
  }
  if(is.null(D)) D <- diag(q)
  
  if(nK != length(var.e)){
    stop("Dimension mismatch between beta, dimension K = ", nK, " and variance terms, length = ", length(var.e))
  }
  if(nK != length(gamma)){
    stop("Dimension mismatch between beta, dimension K = ", nK, " and gamma terms, length = ", length(gamma))
  }
  

  # Baseline covariates -----------------------------------------------------
  id <- 1:n
  time <- rep(0:(ntms-1), n)
  cont <- rnorm(n, 0, 1)
  bin <- rbinom(n, 1, 0.5)
  # Data matrices
  Z <- make.design.matrix(time, order)
  X <- cbind(Z, rep(cont, each = ntms), rep(bin, each = ntms))
  K <- cbind(cont, bin)

  # Random effects
  b <- MASS::mvrnorm(n, rep(0, q), Sigma = D)
  bl <- b[rep(1:n, each = ntms), ]
  
  # Longitudinal ------------------------------------------------------------
  Zik <- make.design.matrix(0:(ntms-1), order) # always the same
  qk <- split(seq(q), cut(seq_along(seq(q)), nK, labels = F))
  Ytest <- list()
  for(i in 1:n){
    Ytest[[i]] <- matrix(NA, nc = nK, nr = ntms)
    Xi <- cbind(Zik, rep(cont[i], each = ntms), rep(bin[i], each = ntms))
    for(kk in 1:nK){
      Ytest[[i]][, kk] <- Xi %*% beta[kk,] + Zik %*% b[i, qk[[kk]]] + rnorm(ntms, 0, sqrt(var.e[kk]))
    }
  }
  Yk <- do.call(rbind, Ytest)
  colnames(Yk) <- paste0("Y.", 1:nK)
  
  # Longitudinal ------------------------------------------------------------
  # qk <- split(seq(q), cut(seq_along(seq(q)), nK, labels = F))
  # Yk <- matrix(NA, nrow = n*ntms, ncol = nK)
  # for(kk in 1:nK){
  #   Xb <- X %*% beta[kk, ]
  #   Zb <- tcrossprod(Z, bl[, qk[[kk]]])
  #   Yk[, kk] <- Xb + colSums(Zb) + rnorm(n*ntms, 0, sqrt(var.e[kk]))
  # }
  # colnames(Yk) <- paste0("Y.", 1:nK)
  
  # Survival ----------------------------------------------------------------
  Keta <- K %*% eta
  # Define gridsteps
  grid.steps <- seq(0, tau, dt)[-1]
  Zdt <- make.design.matrix(grid.steps, order)
  # Hazard
  bl.haz <- exp(Keta) %*% exp(theta0 + theta1 * grid.steps) # subj x grid
  # gamma terms
  gamma.b <- list()
  for(i in 1:nK){
    Zdtb <- tcrossprod(b[, qk[[i]]], Zdt) # subj x grid
    gamma.b[[i]] <- gamma[i] * Zdtb
  }
  gamma.b <- exp(Reduce('+', gamma.b))
  # lambda
  l0 <- (gamma.b * bl.haz) * dt
  
  # Matrix of candidate times
  dims <- dim(l0)
  candidate.times <- matrix(grid.steps, nr = dims[1], nc = dims[2], byrow = T)
  U <- matrix(runif(prod(dims)), nr = dims[1], nc = dims[2])
  
  # Generate survival times
  candidate.times[l0 < U] <- tau
  surv.time <- apply(candidate.times, 1, min)

  # Censoring
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(surv.time, cens.time) # Define 'final' survtime
  
  # Status flag
  status <- rep(0, n)
  is.censored <- cens.time < surv.time
  status[which(surv.time != tau & !is.censored)] <- 1 # Failure flag
  
  # Output Dataframes
  surv.data <- data.frame(id, cont, bin, survtime, status)
  long.data <- data.frame(id = rep(id, each = ntms), time,
                          cont = rep(cont, each = ntms),
                          bin = rep(bin, each = ntms),
                          Yk)
  
  # Check for repeats
  rep.counter <- 0
  repeated <- any(count(surv.data %>% filter(status == 1) %>% distinct(id, survtime), survtime)$n > 1)
  while(repeated){
    surv.data <- noreps(surv.data, dt)
    repeated <- any(count(surv.data %>% filter(status == 1) %>% distinct(id, survtime), survtime)$n > 1)
    rep.counter <- rep.counter + 1
  }
  message("\nSurvival times staggered ", rep.counter, " time(s)...")
  
  dat <- left_join(long.data, surv.data %>% select(id, survtime, status), "id") %>% 
    filter(time < survtime)
  
  message('Time order ', order, ', total number of random effects: ', q)
  message(round(sum(surv.data$status)/n*100, 2), " % failure rate")
  
  return(list(
    dat = dat, survdat = surv.data
  ))
}

noreps <- function(survdata, dt){
  survdata.fail <- survdata %>% 
    filter(status == 1) %>% 
    distinct(id, survtime) %>% 
    arrange(survtime) %>% 
    mutate(diffs = diff(c(0, survtime))) %>% 
    mutate(survtime.new = ifelse(diffs == 0, survtime + dt, survtime)) 
  
  survdata.out <- left_join(survdata, survdata.fail %>% select(id, survtime.new), "id") %>% 
    mutate(survtime = ifelse(!is.na(survtime.new), survtime.new, survtime)) %>%
    arrange(id) %>% 
    select(-survtime.new)
  survdata.out 
}

