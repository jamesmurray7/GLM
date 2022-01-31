#' ####
#' Simulating data under a quadratic(+) random effects structure
#' ----
#' Arguments
#' n: Number of subjects
#' ntms: maximum profile length
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

quad.simData <- function(n = 150, ntms = 5, beta = matrix(1, nr = 3, nc = 5),
                         eta = c(0, 1), gamma = c(0.5, 1, -0.5), var.e = rep(0.25, 3),
                         D = NULL, theta0 = -6, theta1 = 0.25, cens.rate = exp(-3), dt = 0.01){
  
  nK <- nrow(beta) # K
  q <- nK * 3      # total number of random effects.
  tau <- (ntms - 1) + 0.1 # Truncation time
  
  # Checks ----
  if(is.null(D)) D <- diag(q)
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){
    stop("Covariance matrix must be positive semi-definite")
  }
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
  X <- cbind(1, time, time^2, rep(cont, each = ntms), rep(bin, each = ntms))
  Z <- cbind(1, time, time ^ 2)
  K <- cbind(cont, bin)
  
  # Get block diagonal matrices - not used in rest of program -- couldnt work it out!
  # Zk <- as.matrix(Matrix::bdiag(replicate(nK, Z, simplify = F)))
  # Xk <- as.matrix(Matrix::bdiag(replicate(nK, X, simplify = F)))
  
  # Random effects
  b <- MASS::mvrnorm(n, rep(0, q), Sigma = D)
  bl <- b[rep(1:n, each = ntms), ]
  
  # Longitudinal ------------------------------------------------------------
  Zik <- cbind(1, 0:(ntms-1), (0:(ntms-1))^2)
  qk <- split(seq(q), cut(seq_along(seq(q)), nK, labels = F))
  Ytest <- list()
  for(i in 1:n){
    Ytest[[i]] <- matrix(NA, nc = nK, nr = ntms)
    Xi <- cbind(1, 0:(ntms-1), (0:(ntms-1))^2, rep(cont[i], each = ntms), rep(bin[i], each = ntms))
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
  Zdt <- cbind(1, grid.steps, grid.steps^2)
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
  message("\nSurvival times staggered ", rep.counter, " time(s)...\n")
  
  dat <- left_join(long.data, surv.data %>% select(id, survtime, status), "id") %>% 
    filter(time < survtime)
  
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

# beta <- rbind(c(0, 1, 1, 1),
#               c(0, -1, 0, 0.5),
#               c(0, 0, 0.5, -0.5))
# D <- diag(9)
# D[1, 1] <- D[4, 4] <- D[7, 7] <- 0.5^2
# D[2, 2] <- D[5, 5] <- D[8, 8] <- 0.2^2
# D[3, 3] <- D[6, 6] <- D[9, 9] <- 0.1^2
# 
# testdata <- quad.simData(D = D, beta = beta, theta0 = -3)
# testfit <- lmer(Y.1~time+cont+bin+(1+time+I(time^2)|id), testdata$dat)          # These two are the same
# testpoly <- lmer(Y.1~time+cont+bin+(1+poly(time, 2, raw = T)|id), testdata$dat) # ----------------------
# 
# source("../joineRML-compare/funcs.R")
# 
# mjoint(
#   formLongFixed = getFormLongFixed(3),
#   formLongRandom = list(
#     "Y.1" = ~ 1 + time + I(time^2)|id,
#     "Y.2" = ~ 1 + time + I(time^2)|id,
#     "Y.3" = ~ 1 + time + I(time^2)|id
#   ),
#   formSurv = Surv(survtime, status) ~ cont + bin,
#   timeVar = "time", data = testdata$dat,
#   survData = testdata$survdat, pfs = F, verbose = T,
#   control = list(type = "sobol", convCrit = "abs", tol0 = 5e-3, tol.em = 5e-3)
# ) -> jmtest
# 
# summary(jmtest)

# test$dat %>% 
#   select(id, time, Y.1:Y.3) %>% 
#   pivot_longer(Y.1:Y.3) %>% 
#   ggplot(aes(x = time, y = value, group = id)) + 
#   geom_line(alpha = .33) + 
#   facet_wrap(~name, scales = 'free')
