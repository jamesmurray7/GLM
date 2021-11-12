#' ####
#' Simulating simple Zero-Inflated Poisson (ZIP) with random effects and fitting using glmmTMB
#' ###

rm(list=ls())
# data characteristics
n <- 200
mi <- 5
N <- n * mi
# x1, x2, x3 (time cont bin)
x1 <- 0:(mi - 1)
x2 <- rnorm(n)
x3 <- rbinom(n, 1, .5)  
# Parameters
# linear predictor coefficients
beta <- c(.1, .05, -.1, .1)
# coefficient for probability of zero occurring
alpha <- c(-5, .5)
# Covariance matrix D
D <- diag(4)
diag(D) <- c(.5^2, .2^2, .15^2, .1^2)
eigen(D)
# Simulate REs
REs <- MASS::mvrnorm(n, c(0, 0, 0, 0), D)
b <- REs[,1:2]; bstar <- REs[,3:4]

Yi <- ww <- list()
for(i in 1:n){
  # data matrices
  # bstar[i,] <- 0
  W <- Z <- cbind(1, x1)
  X <- cbind(1, x1, x2[i], x3[i])
  wi <- exp(W %*% alpha + Z %*% bstar[i,])/(exp(W %*% alpha + Z %*% bstar[i,]) + 1)
  lambdai <- exp(X %*% beta + Z %*% b[i,])
  Yij <- numeric(length(x1))
  for(j in seq_along(wi)){
    Yij[j] <- ifelse(rbinom(1, 1, wi[j]) == 1, 0, rpois(1, lambdai[j]))
  }
  Yi[[i]] <- Yij
  ww[[i]] <- wi
}
Y <- do.call(c, Yi)

df <- data.frame(
  id = rep(1:n, each = mi),
  x1 = rep(x1, n),
  x2 = rep(x2, each = mi),
  x3 = rep(x3, each = mi),
  Y = Y
)

library(glmmTMB)
fit <- glmmTMB(Y ~ x1 + x2 + x3 + (1+x1|id),
               data = df, family = poisson, ziformula = ~ 1 + x1 + (1+x1|id))#,
               #control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))
fit
VarCorr(fit)$cond$id


# Quick simulation --------------------------------------------------------
simData <- function(){
  n <- 200
  mi <- 5
  N <- n * mi
  # x1, x2, x3 (time cont bin)
  x1 <- 0:(mi - 1)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, .5)  
  # coefficients
  beta0 <- 0.1; beta1 <- 0.05; beta2 <- -0.1; beta3 <- 0.1
  # REs
  D <- matrix(c(0.5^2, 0, 0, 0.2^2),2,2)
  REs <- MASS::mvrnorm(n, c(0,0), D)
  b0 <- REs[,1]; b1 <- REs[,2]
  
  Yi <- list()
  for(i in 1:n){
    lambdai <- beta0+b0[i] + x1 * (beta1 + b1[i]) + x2[i] * beta2 + x3[i] * beta3
    Yi[[i]] <- rpois(mi, lambda = exp(lambdai))
  }
  Y <- do.call(c, Yi)
  
  df <- data.frame(
    id = rep(1:n, each = mi),
    x1 = rep(x1, n),
    x2 = rep(x2, each = mi),
    x3 = rep(x3, each = mi),
    b0 = rep(b0, each = mi),
    b1 = rep(b1, each = mi),
    Y = Y
  )
  
  df
}

fitData <- function(){
  df <- simData()
  glmfit <- glmer(Y ~ x1 + x2 + x3 + (1 + x1|id), data = df, family = poisson)
  c(fixef(glmfit), diag(VarCorr(glmfit)$id))
}

fits <- replicate(100, fitData(), simplify = F)
fits.bind <- as.data.frame(do.call(rbind, fits))
names(fits.bind) <- c(paste0("beta_", 0:3), "D11", "D22")

library(tidyverse)
gather(fits.bind) %>% 
  ggplot(aes(x = value)) + 
  geom_density() + 
  facet_wrap(~key, scales = "free") + 
  theme_bw()
# Verified!