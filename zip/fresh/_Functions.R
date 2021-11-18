#' #####
#' Helper functions for a (+ zero-inflated) Poisson 
#' ----
#' For time being these just focussed on ZIP part, survival added in due course!
#' #####

# Data simulation ---------------------------------------------------------
# 1. ZIP. Assumes RE are random intercept only.
simData_zip <- function(n, ntms, beta, alpha, D = diag(2)){
  time <- 0:(ntms - 1)
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  # Design matrices
  # Design matrices for the fixed and random effects in the non-zero inflated part.
  X <- model.matrix(~ time + cont + bin, data = df)
  Z <- model.matrix(~ 1, data = df)
  # Design matrices for the fixed and random effects in the zero inflated (ZI) part.
  Xzi <- model.matrix(~ time, data = df)
  Zzi <- model.matrix(~ 1, data = df)
  
  # Random Effects,
  b <- MASS::mvrnorm(n, rep(0, 2), Sigma = D);
  # bzi <- MASS::mvrnorm(n, rep(0, 2), Sigma = D[[2]])
  
  # Linear predictors
  eta <- X %*% beta + rowSums(Z * b[df$id, 1, drop=F])
  etazi <- Xzi %*% alpha + rowSums(Zzi * b[df$id, 2, drop=F])
  
  # Outcome y
  y <- rpois(n * ntms, lambda = exp(eta))
  y[as.logical(rbinom(n * ntms, 1, prob = plogis(etazi)))] <- 0
  df$y <- y
  
  message(sum(df$y == 0)/length(df$y)*100,'% zeroes')
  
  df
}

# Simulate under a (very simple) joint model
simData_zip_joint <- function(n, ntms, beta, alpha, D = diag(2),
                              theta=c(-6, 0.2), surv.eta = c(0, 1), 
                              gamma = 0.5, cens.rate = exp(-3.5)){
  
  # Necessary parameters
  id <- 1:n
  time <- 0:(ntms-1); tau <- (ntms - 1) + 0.1 # time variable and truncation time tau.
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  # Design matrices
  # Design matrices for the fixed and random effects in the non-zero inflated part.
  X <- model.matrix(~ time + cont + bin, data = df)
  Z <- model.matrix(~ 1, data = df)
  # Design matrices for the fixed and random effects in the zero inflated (ZI) part.
  Xzi <- model.matrix(~ time, data = df)
  Zzi <- model.matrix(~ 1, data = df)
  
  # Random Effects,
  b <- MASS::mvrnorm(n, rep(0, 2), Sigma = D);
  # bzi <- MASS::mvrnorm(n, rep(0, 2), Sigma = D[[2]])
  
  # Linear predictors
  eta <- X %*% beta + rowSums(Z * b[df$id, 1, drop=F])
  etazi <- Xzi %*% alpha + rowSums(Zzi * b[df$id, 2, drop=F])
  
  # Outcome y
  y <- rpois(n * ntms, lambda = exp(eta))
  y[as.logical(rbinom(n * ntms, 1, prob = plogis(etazi)))] <- 0
  df$y <- y; message(sum(df$y == 0)/length(df$y)*100,'% zeroes')
  
  # Simulating survival times
  theta0 <- theta[1]; theta1 <- theta[2]
  K <- cbind(cont, bin)
  Keta <- K %*% surv.eta
  U <- runif(n)
  
  # denom <- theta1 + b1 %*% gamma                     # these for if I put intercept + slope back in!
  # rhs <- (theta1 + b1 %*% gamma) * log(U)/(exp(theta0 + Keta + b0 %*% gamma))   
  denom <- theta1
  rhs <- theta1 * log(U)/(exp(theta0 + Keta + rowSums(b * gamma)))
  
  t <- suppressWarnings(log((1 - rhs))/denom)
  t[is.nan(t)] <- tau
  
  # Collect survival times, and generate cenors times.
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  # Output Dataframes
  surv.data <- data.frame(id, survtime, status)
  long.data <- df
  
  out.data <- dplyr::left_join(df, surv.data, by = 'id')
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  return(list(data =  out.data, 
         surv.data = surv.data))
  
}

# Families ----------------------------------------------------------------

# 1. Poisson
po <- function(){
  link <- stats::make.link('log')
  log.dens <- function(y, eta){
    mu <- link$linkinv(eta)
    y * eta - mu - lfactorial(y)
  }
  Seta <- function(y, eta){
    mu <- link$linkinv(eta)
    y - mu
  }
  structure(list(family = 'Poisson', 
                 link = link$name,
                 linkfun = link$linkfun,
                 linkinv = link$linkinv,
                 ldens = log.dens,
                 Seta = Seta), class = 'family')
}

# 2. Zero-inflated Poisson ('ZIP')

zip <- function(){
  link <- stats::make.link('log')
  log.dens <- function(y, eta, etazi){ # log density ZIP
    y0 <- y == 0
    y1 <- y > 0
    mu <- matrix(link$linkinv(eta), nc = 1)
    muzi <- matrix(link$linkinv(etazi), nc = 1)
    mu0 <- mu[y0, ]; muzi0 <- muzi[y0, ]
    mu1 <- mu[y1, ]
    out <- matrix(eta, nc = 1)
    out[y0, ] <- log(muzi0 + exp(-mu0))
    out[y1, ] <- y[y1] * log(mu1) - mu1 - lfactorial(y[y1])
    out <- out - log(1 + muzi)
    out
  }
  Seta <- function(y, eta, etazi){     # Score for \eta (count part)
    y0 <- y == 0; y1 <- y > 0
    out <- matrix(eta, nc = 1)
    eta0 <- eta[y0,]; eta1 <- eta[y1,]
    etazi0 <- etazi[y0,]
    out[y0, ] <- -exp(-exp(eta0)) * exp(eta0) / (exp(-exp(eta0)) + exp(etazi0))
    out[y1, ] <- y[y1] - exp(eta1)
    out
  }
  Setazi <- function(y, eta, etazi){   # Score for \eta_{zi} (ZI part)
    y0 <- y == 0; y1 <- y > 0
    out <- matrix(eta, nc = 1)
    eta0 <- eta[y0,];
    etazi0 <- etazi[y0,]; etazi1 <- etazi[y1,]
    out[y0, ] <- exp(etazi0) / (exp(etazi0)+exp(-exp(eta0))) - exp(etazi0) / (1 + exp(etazi0)) 
    out[y1, ] <- -exp(etazi1) / (1 + exp(etazi1))
    out
  }
  structure(list(family = 'Zero-inflated Poisson', 
                 link = link$name,
                 linkfun = link$linkfun,
                 linkinv = link$linkinv,
                 ldens = log.dens,
                 Seta = Seta, 
                 Setazi = Setazi), class = 'family')
}

# Functions for parameter updates -----------------------------------------
beta_alpha <- function(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, ldens, zi.inds, S, Sz, gh, Sigmai){
  gh <- statmod::gauss.quad.prob(gh, 'normal')
  v <- gh$n; w <- gh$w
  # Split Sigmai into two matrices
  SS <- lapply(split(seq(2), c(1,2)), function(x) Sigmai[x,x])
  # get tau
  tau <- sqrt(diag(tcrossprod(Z %*% SS[[1]], Z)))
  tauzi <- sqrt(diag(tcrossprod(Zzi %*% SS[[2]], Zzi)))
  Sbeta_l <- function(beta, b, Y, X, Z, Xzi, Zzi, alpha, zi.inds, tau, tauzi, w, v){  # Define a function for the beta-score (for central-differencing)
    y0 <- Y == 0; y1 <- Y > 0
    eta <- X %*% beta + Z %*% b[-zi.inds] + tau * v
    etazi <- Xzi %*% alpha + Zzi %*% b[zi.inds] + tauzi * v
    eta0 <- eta[y0,]; eta1 <- eta[y1,]
    etazi0 <- etazi[y0,]
    out <- matrix(eta, nc = 1)
    out[y0, ] <- -exp(-exp(eta0)) * exp(eta0) / (exp(-exp(eta0)) + exp(etazi0))
    out[y1, ] <- Y[y1] - exp(eta1)
    w * crossprod(X, out)
  }
  Salpha_l <- function(alpha, b, Y, X, Z, Xzi, Zzi, beta, zi.inds, tau, tauzi, w, v){  # Define a function for the alpha-score (for central-differencing)
    y0 <- Y == 0; y1 <- Y > 0
    eta <- X %*% beta + Z %*% b[-zi.inds] + tau * v
    etazi <- Xzi %*% alpha + Zzi %*% b[zi.inds] + tauzi * v
    eta0 <- eta[y0,];
    etazi0 <- etazi[y0,]; etazi1 <- etazi[y1,]
    out <- matrix(eta, nc = 1)
    out[y0, ] <- exp(etazi0) / (exp(etazi0)+exp(-exp(eta0))) - exp(etazi0) / (1 + exp(etazi0)) 
    out[y1, ] <- -exp(etazi1) / (1 + exp(etazi1))
    w * crossprod(Xzi, out)
  }
  # linear predictors and score, hessians
  etazi <- eta <- Seta <- Setazi <- Hbeta <- Halpha <- list()
  for(l in 1:length(w)){
    eta[[l]] <- X %*% beta + Z %*% b[-zi.inds] + tau * v[l]
    etazi[[l]] <- Xzi %*% alpha + Zzi %*% b[zi.inds] + tauzi * v[l]
    Seta[[l]] <- w[l] * S(Y, eta[[l]], etazi[[l]])
    Setazi[[l]] <- w[l] * Sz(Y, eta[[l]], etazi[[l]])
    Hbeta[[l]] <- GLMMadaptive:::cd_vec(beta, Sbeta_l, b, Y, X, Z, Xzi, Zzi, alpha, zi.inds, tau, tauzi, w[l], v[l])
    Halpha[[l]] <- GLMMadaptive:::cd_vec(alpha, Salpha_l, b, Y, X, Z, Xzi, Zzi, beta, zi.inds, tau, tauzi, w[l], v[l])
  }
  return(list(Sbeta = crossprod(X, Reduce('+', Seta)),
              Salpha = crossprod(Xzi, Reduce('+', Setazi)),
              Hbeta = Reduce('+', Hbeta),
              Halpha = Reduce('+', Halpha)))
}


# Random effects, b -------------------------------------------------------
.b <- function(){
  logfb <- function(b, Y, X, Z, Xzi, Zzi, 
                    beta, D, alpha, ldens, indzi, S, Sz){
    eta <- X %*% beta + Z %*% b[-indzi]
    etazi <- Xzi %*% alpha + Zzi %*% b[indzi]
    -sum(ldens(Y, eta, etazi)) + 0.5 * tcrossprod(b %*% solve(D), b)
  }
  score.logfb <- function(b, Y, X, Z, Xzi, Zzi, 
                          beta, D, alpha, ldens, indzi, S, Sz){
    eta <- X %*% beta + Z %*% b[-indzi]
    etazi <- Xzi %*% alpha + Zzi %*% b[indzi]
    out <- -crossprod(Z, S(Y, eta, etazi))              # Poisson process
    out <- c(out, -crossprod(Zzi, Sz(Y, eta, etazi)))   # ZIP process
    out + b %*% solve(D) 
  }
  hess <- function(b, Y, X, Z, Xzi, Zzi, 
                   beta, D, alpha, ldens, indzi, S, Sz){
    optim(b, logfb, score.logfb, Y = Y, X = X, Z = Z, Xzi = Xzi,
          Zzi = Zzi, beta = beta, D = D, alpha = alpha, ldens = ldens, indzi = indzi, S = S, Sz = Sz, method = 'BFGS', hessian = T,
          control = list(reltol = 1e-3))$hessian
  }
  list(logfb = logfb,
       Sb = score.logfb,
       H = hess)
}

beta_alpha_new <- function(b, Y, X, Z, Xz, Zz, beta, D, alpha, Sigmai, gh){
  gh <- statmod::gauss.quad.prob(gh, 'normal')
  w <- gh$w; v <- gh$n
  SS <- lapply(split(seq(2), c(1,2)), function(x) as.matrix(Sigmai[x,x]))
  # Scores
  out <- Sbeta_alpha(b, Y, X, Z, Xz, Zz, beta, D, alpha, 2, SS, w, v)
  # Hessians
  Ha <- Hb <- list()
  for(l in 1:length(w)){
    Ha[[l]] <- cd(alpha, Salphal, Y, X, Z, Xz, Zz, beta, D, b, 2, SS, w[l], v[l])
    Hb[[l]] <- cd(beta , Sbetal,  Y, X, Z, Xz, Zz, b, D, alpha, 2, SS, w[l], v[l])
  }
  out[["Hbeta"]] <- Reduce('+', Hb)
  out[["Halpha"]] <- Reduce('+', Ha)
  out
}

