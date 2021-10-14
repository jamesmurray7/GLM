#' ####
#' Simulating simple Poisson structure with random effects
#' and fitting using lme4::glmer...
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

library(lme4)
glmfit <- glmer(Y ~ x1 + x2 + x3 + (1 + x1|id), data = df, family = poisson)
VarCorr(glmfit)$id


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