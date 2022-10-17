#' Same as "first test.R" but with simulated data, mapply steps and so on.

rm(list=ls())
library(Rcpp)
library(glmmTMB)
library(RcppArmadillo)
sourceCpp('funs.cpp')
source('simData.R')

# Poisson body simulation
data <- simData_joint()$data
n <- 250
X <- Z <- Y <- vector('list', n)
for(i in 1:n){
  X[[i]] <- model.matrix(~time+cont+bin, data[data$id == i,])
  Z[[i]] <- model.matrix(~1, data[data$id == i,])
  Y[[i]] <- matrix(data[data$id == i, 'Y'], nc = 1)
}

inits <- glmmTMB(Y~time+cont+bin + (1|id), data = data, family = genpois, dispformula=~1)
beta <- fixef(inits)$cond
b <- ranef(inits)$cond$id
b <- lapply(1:n, function(i) b[i,])
phi <- 0
D <- VarCorr(inits)$cond$id

ll_beta <- function(beta, X, Z, b, Y, phi){
  ll_genpois(exp(X %*% beta + Z %*% b), phi, Y)
}
ll_beta(beta, X = X[[1]], Z = Z[[1]], b = b[[1]], Y = Y[[1]], phi = phi)

# Test for the update to \beta...
pracma_beta <- mapply(function(X, Z, b, Y){
  a <- pracma::grad(ll_beta, beta, X = X, Z = Z, b = b, Y = Y, phi = phi)
  b <- pracma::hessian(ll_beta, beta, X = X, Z = Z, b = b, Y = Y, phi = phi)
  list(grad = a, hess = b)
}, X = X, Z = Z, b = b, Y = Y, SIMPLIFY = F)

analy_beta <- mapply(function(X, Z, b, Y){
  mu <- exp(X %*% beta + Z %*% b)
  grad <- numeric(ncol(X))
  for(p in 1:ncol(X)){
    xp <- X[, p, drop = F]
    # for(j in 1:length(xp)){
    # this <- xp[j] + (Y[j] - 1) * (xp[j] * mu[j])/(mu[j] + phi * Y[j]) - xp[j] * mu[j] / (phi + 1)
    # This also works ...
    # print()
    grad[p] <- sum(xp + (Y - 1) * (xp * mu)/(mu + phi * Y) - xp * mu / (phi + 1))
  }
  H <- matrix(0, ncol(X), ncol(X))
  for(j in 1:nrow(X)){
    xj <- X[j,,drop=T]
    a <- (phi * (Y[j] - 1.) * Y[j]/((mu[j] + phi * Y[j])^2) - mu[j]/(1+phi)) * xj %*% t(xj)
    H <- H + a
  }
  list(grad=grad,hess=H)
}, X = X, Z = Z, b = b, Y = Y, SIMPLIFY = F)

beta - solve(Reduce('+', lapply(analy_beta, el, 2)), Reduce('+', lapply(analy_beta, el, 1)))


# \phi -- > Scalar dispersion only (I think).
ll_phi <- function(phi, X, Z, b, Y, beta){
  ll_genpois(exp(X %*% beta + Z %*% b), phi, Y)
}
pracma_phi <- mapply(function(X, Z, b, Y){
  a <- pracma::grad(ll_phi, phi, X = X, Z = Z, b = b, Y = Y, beta = beta)
  b <- pracma::hessian(ll_phi, phi, X = X, Z = Z, b = b, Y = Y, beta = beta)
  list(grad = a, hess = b)
}, X = X, Z = Z, b = b, Y = Y, SIMPLIFY = F)

analy_phi <- mapply(function(X,Z,b,Y){
  mu <- exp(X %*% beta + Z %*% b)
  grad <- Hess <- 0
  for(j in 1:length(mu)){
    grad <- grad + (
      (Y[j]-1) * Y[j] / (mu[j] + phi * Y[j]) - Y[j]/(1+phi) + (mu[j] - Y[j])/((phi + 1)^2)
    )
    Hess <- Hess + (
      2 * (Y[j] - mu[j]) / ((phi + 1)^3) + Y[j]/((1+phi)^2) - 
        Y[j]^2 * (Y[j]-1) / ((mu[j] + phi * Y[j])^2)
    )
  }
  list(grad = grad, hess = Hess)
}, X = X, Z = Z, b = b, Y = Y, SIMPLIFY = F)

phi - sum(sapply(analy_phi, el, 1))/sum(sapply(analy_phi, el, 2))

# Raneffs \b
ll_b <- function(b, X, Z, Y, beta, phi){
  ll_genpois(exp(X %*% beta + Z %*% b), phi, Y)
}

pracma_b <- mapply(function(X, Z, b, Y){
  a <- pracma::grad(ll_b, b, X = X, Z = Z, Y = Y, beta = beta, phi = phi)
  b <- pracma::hessian(ll_b, b, X = X, Z = Z, Y = Y, beta = beta, phi = phi)
  list(grad = a, hess = b)
}, X = X, Z = Z, b = b, Y = Y, SIMPLIFY = F)

analy_b <- mapply(function(Z,X,b,Y){
  mu <- exp(X %*% beta + Z %*% b)
  grad <- numeric(ncol(Z))
  for(q in 1:ncol(Z)){
    Zq <- Z[, q, drop = F]
    grad[q] <- sum(Zq + (Y - 1) * (Zq * mu)/(mu + phi * Y) - Zq * mu / (phi + 1))
  }
  H <- matrix(0, ncol(Z), ncol(Z))
  for(j in 1:nrow(Z)){
    zj <- Z[j,,drop=T]
    H <- H + (phi * (Y[j] - 1) * Y[j] / ((mu[j] + phi * Y[j])^2) - mu[j] / (phi + 1) ) * zj %*% t(zj)
  }
  list(grad = grad, hess = H)
}, X = X, Z = Z, b = b, Y = Y, SIMPLIFY = F)

#' ###
#' end
#' ###

joint_density_ddb(b[[1]], X[[1]], Y[[1]], Z[[1]], beta, phi, D, 
                  c(0), matrix(0,2,1), c(0), matrix(0,2,1), 1, c(0,0), 1, .6, c(-.2))
pracma::grad(joint_density, b[[1]], X = X[[1]], Y = Y[[1]], Z = Z[[1]], beta = beta, phi = phi, D = D,
             S = c(0), SS = matrix(0,2,1), Fi = c(0), Fu = matrix(0,2,1), l0i = 1, haz = c(0,0), Delta = 1, gamma = .6, zeta = c(-.2))


