#' #####
#' Testing functions for GLMM negbin test.
#' ####
rm(list=ls())
source('./simData.R')
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
sourceCpp('nb.cpp')

set.seed(123)
test <- simData()
fit <- glmmTMB(Y ~ time + cont + bin + (1+time|id), test, dispformula = ~1, family = nbinom2)
fit

n <- 250
Y <- X <- Z <- K <- list()
for(i in 1:n){
  i.dat <- test[test$id == i, ]
  Y[[i]] <- i.dat$Y
  X[[i]] <- model.matrix(~time+cont+bin, i.dat)
  Z[[i]] <- model.matrix(~time, i.dat)
  K[[i]] <- cbind(unique(i.dat$cont), unique(i.dat$bin))
}

true.ll <- logLik(fit)
beta <- fixef(fit)$cond
theta <- sigma(fit)
D <- matrix(VarCorr(fit)$cond$id,2,2)
b <- as.matrix(ranef(fit)$cond$id)
b <- lapply(1:n, function(i) b[i,])

nb.ll <- function(b, X, Y, Z, beta, theta){
  mu <- exp(X %*% beta + Z %*% b)
  crossprod(Y, log(mu)) - crossprod(theta + Y, log(mu + theta)) + sum(lgamma(theta + Y) - lgamma(Y + 1)) - lgamma(theta) + theta * log(theta)
}

D.ll <- function(b, D){
  -length(b)/2*log(2*pi) - 1/2 * log(det(D)) -1/2 * crossprod(b, solve(D) %*% b)
}

Dll <- ll <- c()
for(i in 1:n){ 
  ll[i] <- nb.ll(b[[i]], X[[i]], Y[[i]], Z[[i]], beta, theta)
  Dll[i] <- D.ll(b[[i]], D)
}
sum(ll) + sum(Dll)
true.ll
