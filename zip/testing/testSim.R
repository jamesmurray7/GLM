#' ####
#' Simulating simple Zero-Inflated Poisson (ZIP) with random effects and fitting using glmmTMB
#' ###

rm(list=ls())
n <- 250 # no. subjects
ntms <- 5 # Max. follow-up time.

# Baseline covariates
time <- 0:(ntms - 1)
cont <- rnorm(n)
bin <- rbinom(n, 1, 0.5)

df <- data.frame(id = rep(1:n, each = ntms),
                 time = rep(time, n),
                 cont = rep(cont, each = ntms),
                 bin = rep(bin, each = ntms))

# Design matrices for the fixed and random effects in the non-zero inflated part.
X <- model.matrix(~ time + cont + bin, data = df)
Z <- model.matrix(~ time, data = df)
# Design matrices for the fixed and random effects in the zero inflated (ZI) part.
Xzi <- model.matrix(~ time, data = df)
Zzi <- model.matrix(~ time, data = df)

# Coefficients
beta <- c(1, 0.05, 0.3, -0.2)
alpha <- c(-.5, .2)
# Covariance matrix D - 2x2, intercept ONLY e.g. Zhu et al (2018)
D11 <- 0.5 # variance of random intercepts non-zero part
D22 <- 0.4 # variance of random intercepts zero part
# Random effects (longit and ZI part)
Dl <- matrix(c(.5^2, 0.05, 0.05, 0.2^2), 2, 2)
Dzi <- matrix(c(.3^2, 0, 0, .1^2), 2, 2)
b <- MASS::mvrnorm(n, rep(0, 2), Sigma = Dl)
bstar <- MASS::mvrnorm(n, rep(0, 2), Sigma = Dzi)

# Linear predictor for the Poisson process
eta.y <- X %*% beta + rowSums(Z * b[df$id,,drop=F])
# Linear predictor for ZI process
eta.zi <- Xzi %*% alpha + rowSums(Zzi * bstar[df$id, , drop = F])
# Simulate Poisson process

df$y <- rpois(n * ntms, lambda = exp(eta.y))
df$y[as.logical(rbinom(n * ntms, 1, prob = plogis(eta.zi)))] <- 0   # Set y=0 with probablility logit(eta.zi)

as.integer(table(df$y)[1])

glmmTMB(y ~ time + cont + bin + (1 + time|id), 
        data = df, 
        family = poisson, 
        ziformula = ~ time + (1 + time|id))

# Quick simulation --------------------------------------------------------
simData <- function(){
  n <- 250 # no. subjects
  ntms <- 5 # Max. follow-up time.
  
  # Baseline covariates
  time <- 0:(ntms - 1)
  cont <- rnorm(n)
  bin <- rbinom(n, 1, 0.5)
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  # Design matrices for the fixed and random effects in the non-zero inflated part.
  X <- model.matrix(~ time + cont + bin, data = df)
  Z <- model.matrix(~ time, data = df)
  # Design matrices for the fixed and random effects in the zero inflated (ZI) part.
  Xzi <- model.matrix(~ time, data = df)
  Zzi <- model.matrix(~ time, data = df)
  
  # Coefficients
  beta <- c(1, 0.05, 0.3, -0.2)
  alpha <- c(-.5, .2)
  # Covariance matrix D - 2x2, intercept ONLY e.g. Zhu et al (2018)
  D11 <- 0.5 # variance of random intercepts non-zero part
  D22 <- 0.4 # variance of random intercepts zero part
  # Random effects (longit and ZI part)
  Dl <- matrix(c(.5^2, 0.05, 0.05, 0.2^2), 2, 2)
  Dzi <- matrix(c(.3^2, 0, 0, .1^2), 2, 2)
  b <- MASS::mvrnorm(n, rep(0, 2), Sigma = Dl)
  bstar <- MASS::mvrnorm(n, rep(0, 2), Sigma = Dzi)
  
  # Linear predictor for the Poisson process
  eta.y <- X %*% beta + rowSums(Z * b[df$id,,drop=F])
  # Linear predictor for ZI process
  eta.zi <- Xzi %*% alpha + rowSums(Zzi * bstar[df$id, , drop = F])
  # Simulate Poisson process
  
  df$y <- rpois(n * ntms, lambda = exp(eta.y))
  df$y[as.logical(rbinom(n * ntms, 1, prob = plogis(eta.zi)))] <- 0   # Set y=0 with probablility logit(eta.zi)
  df
}

fitData <- function(){
  df <- simData()
  fit <- glmmTMB(y ~ time + cont + bin + (1 + time|id), 
          data = df, 
          family = poisson, 
          ziformula = ~ time + (1 + time|id))
  return(fit)
}

fits <- replicate(100, fitData(), simplify = F)

ee <- function(fit){
  # Fixed effects
  beta <- fixef(fit)$cond; names(beta) <- paste0('beta_', names(beta))
  alpha <- fixef(fit)$zi; names(alpha) <- paste0('alpha_', names(alpha))
  # VarCorr
  Dl <- diag(matrix(VarCorr(fit)$cond$id,2,2)); names(Dl) <- paste0('Dl_', c(11, 22))
  Dz <- diag(matrix(VarCorr(fit)$zi$id,2,2)); names(Dz) <- paste0('Dz_', c(11, 22))
  c(beta, alpha, Dl, Dz)
}


library(tidyverse)
coeffs <- do.call(rbind, lapply(fits, ee))

targets <- data.frame(
  parameter = colnames(coeffs),
  target = c(1, 0.05, 0.3, -0.2,
             -.5, .2, .5^2, .2^2, .3^2, .1^2)
)

coeffs %>% 
  as_tibble %>% 
  pivot_longer(cols = `beta_(Intercept)`:`Dz_22`, names_to = 'parameter', values_to = 'estimate') %>% 
  left_join(., targets, 'parameter') %>% 
  ggplot(aes(x = estimate)) + 
  geom_vline(aes(xintercept = target), lty = 3, colour = 'grey50') + 
  geom_density() + 
  facet_wrap(~parameter, scales = 'free') +
  theme_bw()
# Looks okay
ggsave('~/Downloads/ZIP-simplots.png')
