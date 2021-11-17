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


