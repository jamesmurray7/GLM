#' #######
#' Bespoke functions to find pmf of, and then simulate random deviates from a
#' mean-parameterised specification of the Conway-Maxwell Poisson Distn.
#' ---
#' Approximately half of n subjects are simulated with \delta[1] and the other with \delta[2]
#' #######
# rm(list=ls())
# old.ls <- ls()
source('../mpcmp-copy/mpcmp-copy.R')
# to.remove <- setdiff(ls(), old.ls)

simData_joint <- function(n = 250, ntms = 10, summax = 100,  fup = 3, 
                          beta = c(2, -0.1, 0.1, -0.2),               # (int, time, cont, bin)
                          delta = c(0.8, -0.3),                       # (x, y), roughly half get dispersion exp(x) and the other exp(y).
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
  
  b <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = D)
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
  
  df <- merge(df, surv.data, by = 'id')
  df <- df[df$time < df$survtime, ]
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  # Design matrices 
  X <- model.matrix(~ time + cont + bin, data = df)
  Z <- model.matrix(~ time, data = df)
  
  #' Simulate CMP response //
  eta <- X %*% beta + rowSums(Z * b[df$id, ]) # Linear predictor
  deltas <- sample(delta, n, T)               # Sample dispersion from c(x, y))
  
  mu <- exp(eta)                              # Mean parameterisation
  
  Y <- vector('list', n)
  pb <- utils::txtProgressBar(max = n, style = 3)
  for(i in 1:n){
    ind <- which(df$id==i)
    nu <- rep(exp(deltas[i]), length(ind))
    Y[[i]] <- rcomp(n = length(ind), mu = mu[ind], nu = nu, summax = summax)
    utils::setTxtProgressBar(pb, i)
  }
  df$Y <- do.call(c, Y)
  cat('\n')
  list(data = df, 
       surv.data = df[!duplicated(df[,'id']), c('id', 'cont', 'bin', 'survtime', 'status')],
       true.deltas = setNames(deltas, 1:n))
}


# Quicker version ---------------------------------------------------------
.rcomp <- function(lambda, nu, summax){
  x <- vector('numeric', length(lambda))
  U <- runif(length(lambda))
  warn <- F
  for(i in 1:length(lambda)){
    if (lambda[i] == 0) {
      x[i] <- 0
    } else if (lambda[i] < 0 | nu[i] <= 0) {
      x[i] <- NA
      warn <- TRUE
    } else {
      y <- 0
      dc <- cmp_pmf_scalar(0:(summax), lambda[i], nu[i], summax)
      py <- dc[y + 1]
      while (py <= U[i]) {
        y <- y + 1
        py <- py + dc[y + 1]; 
      }
      x[i] <- y
    }
  }
  if (warn) {
    warning("NAs produced")
  }
  return(x)
}


simData_joint2 <- function(n = 250, ntms = 10, summax = 100,  fup = 3, 
                           beta = c(2, -0.1, 0.1, -0.2),               # (int, time, cont, bin)
                           delta = c(0.8, -0.3),                       # (x, y), roughly half get dispersion exp(x) and the other exp(y).
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
  
  b <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = D)
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
  
  df <- merge(df, surv.data, by = 'id')
  df <- df[df$time < df$survtime, ]
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  # Design matrices 
  X <- model.matrix(~ time + cont + bin, data = df)
  Z <- model.matrix(~ time, data = df)
  
  #' Simulate CMP response //
  eta <- X %*% beta + rowSums(Z * b[df$id, ]) # Linear predictor
  deltas <- sample(delta, n, T)               # Sample dispersion from c(x, y))
  
  mu <- exp(eta)                              # Mean parameterisation
  
  Y <- vector('list', n)
  pb <- utils::txtProgressBar(max = n, style = 3)
  for(i in 1:n){
    ind <- which(df$id==i)
    nu <- rep(exp(deltas[i]), length(ind))
    lam <- lambda_appx(mu[ind], nu, summax)
    Y[[i]] <- .rcomp(lam, nu, summax)
    # Y[[1]] <- rcomp(n = length(ind), nu = nu, lambda = lam, summax = summax)
    utils::setTxtProgressBar(pb, i)
  }
  df$Y <- do.call(c, Y)
  cat('\n')
  list(data = df, 
       surv.data = df[!duplicated(df[,'id']), c('id', 'cont', 'bin', 'survtime', 'status')],
       true.deltas = setNames(deltas, 1:n))
}


# Check raw subject-specific dispersion -----------------------------------
check.disps <- function(S){ # S a half.and.half (S)imulated data list.
  true.disps <- S$true.deltas
  # Tabulate and print.
  tt <- table(true.disps)
  cat(sprintf("%d with true dispersion %s and %d with %s.\n", 
              tt[1], names(tt[1]), tt[2], names(tt[2])))
  
  # Join true dispersion on  
  d <- S$data
  true.disps2 <- data.frame(id = 1:length(true.disps), true = true.disps)
  d <- merge(d, true.disps2, 'id')
  
  ids <- with(d, tapply(id, true, unique))
  
  # Work out variance/mean
  calculated.disps <- lapply(ids, function(iii){
    out <- numeric(length(iii))
    for(i in seq_along(iii)){
      d_i <- d[d$id == iii[i], ]
      out[i] <- var(d_i$Y)/mean(d_i$Y)
    }
    out[!(is.na(out) | is.infinite(out))]
  })
  
  cat('Summaries for calculated var(Y)/mean(Y) on true delta values:\n')
  
  print(lapply(calculated.disps, summary))
}

check.deltas.against.disps <- function(S, d){ # S a half.and.half (S)imulated data set
  true.disps <- S$true.deltas                 # d subject estimates from get.delta.inits
  ests <- d$subject.estimates
  
  inds <- which(ests != 0)
  true.inds <- true.disps[inds]
  ests.inds <- ests[inds]
  
  # Partition into two true deltas and plot visual
  aaa <- unique(true.inds)
  par(mfrow=c(1,2))
  plot(ests.inds[true.inds == aaa[1]], main = bquote('True delta: ' ~ .(aaa[1])),
       ylab = 'estimate', xaxt = 'n', xlab = '')
  abline(h = aaa[1], col = 'red', lty = 5)
  plot(ests.inds[true.inds = aaa[2]], main = bquote('True delta: ' ~ .(aaa[2])),
       ylab = 'estimate', xaxt = 'n', xlab = '')
  abline(h = aaa[2], col = 'red', lty = 5)
  par(mfrow=c(1,1))
  
  return(true.inds - ests.inds)
}
