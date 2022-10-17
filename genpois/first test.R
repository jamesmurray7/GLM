library(Rcpp)
library(RcppArmadillo)
sourceCpp('funs.cpp')

# Poisson body simulation
beta <- c(2, -.1, .1, -.2)
X <- cbind(1, 1:10, rnorm(1), rbinom(1,1,.5))
Z <- cbind(1, 1:10)
b <- rnorm(2, sd = .2)
Y <- rpois(10, exp(X %*% beta + Z %*% b))
phi <- .3

ll <- function(beta, X, Z, b, Y, phi){
  ll_genpois(exp(X %*% beta + Z %*% b), phi, Y)
}
ll(beta, X = X, Z = Z, b = b, Y = Y, phi = rep(phi, 10))

# Gradient test
pracma::grad(ll, beta, X = X, Z = Z, b = b, Y = Y, phi = rep(phi, 10))
target <- pracma::hessian(ll, beta, X = X, Z = Z, b = b, Y = Y, phi = rep(phi, 10))

out <- numeric(ncol(X))
for(p in 1:ncol(X)){
  xp <- X[, p, drop = F]
  # for(j in 1:length(xp)){
    # this <- xp[j] + (Y[j] - 1) * (xp[j] * mu[j])/(mu[j] + phi * Y[j]) - xp[j] * mu[j] / (phi + 1)
    # This also works ...
    # print()
  out[p] <- sum(xp + (Y - 1) * (xp * mu)/(mu + phi * Y) - xp * mu / (phi + 1))
  # }
}
out


# Hessian
H <- matrix(0, ncol(X), ncol(X))
for(p in 1:ncol(X)){
  # xp <- X[,p,drop=F]
  for(j in 1:nrow(X)){
    xj <- X[j,,drop=T]
    one <- ((phi * Y[j]^2 - Y[j]) / ((mu[j] + phi * Y[j])^2)) * xj %*% t(xj)
    two <- (mu[j]/(phi + 1)) * xj %*% t(xj)
    H <- H + (one - two)
  }
}

H <- matrix(0, ncol(X), ncol(X))
for(p in 1:ncol(X)){
  # xp <- X[,p,drop=F]
  for(j in 1:nrow(X)){
    xj <- X[j,,drop=T]
    one <- ((phi * (Y[j] - 1) * Y[j]) / ((mu[j] + phi * Y[j])^2) - mu[j]/(phi+1)) * xj %*% t(xj)
    H <- H + one
  }
}


# \phi -- > Scalar dispersion only (I think).
ll_phi <- function(phi, X, Z, b, Y, beta){
  ll_genpois(exp(X %*% beta + Z %*% b), phi, Y)
}
pracma::grad(ll_phi, phi, X = X, Z = Z, b = b, Y = Y, beta = beta)
pracma::hessian(ll_phi, phi, X = X, Z = Z, b = b, Y = Y, beta = beta)

grad <- Hess <- numeric(1); 
for(j in 1:length(mu)){
  grad <- grad + (
    (Y[j]-1) * Y[j] / (mu[j] + phi * Y[j]) - Y[j]/(1+phi) + (mu[j] - Y[j])/((phi + 1)^2)
  )
  Hess <- Hess + (
    2 * (Y[j] - mu[j]) / ((phi + 1)^3) + Y[j]/((1+phi)^2) - 
      Y[j]^2 * (Y[j]-1) / ((mu[j] + phi * Y[j])^2)
  )
}

# Raneffs \b
ll_b <- function(b, X, Z, Y, beta, phi){
  ll_genpois(exp(X %*% beta + Z %*% b), phi, Y)
}
pracma::grad(ll_b, b, X = X, Z = Z, Y = Y, beta = beta, phi = phi)
pracma::hessian(ll_b, b, X = X, Z = Z, Y = Y, beta = beta, phi = phi)

grad.b <- numeric(ncol(Z))
for(q in 1:ncol(Z)){
  Zq <- Z[, q, drop = F]
  # for(j in 1:length(xp)){
  # this <- xp[j] + (Y[j] - 1) * (xp[j] * mu[j])/(mu[j] + phi * Y[j]) - xp[j] * mu[j] / (phi + 1)
  # This also works ...
  # print()
  grad.b[q] <- sum(Zq + (Y - 1) * (Zq * mu)/(mu + phi * Y) - Zq * mu / (phi + 1))
}


Hess.b <- matrix(0, ncol(Z), ncol(Z))
for(j in 1:nrow(Z)){
  zj <- Z[j,,drop=T]
  Hess.b <- Hess.b + (phi * (Y[j] - 1) * Y[j] / ((mu[j] + phi * Y[j])^2) - mu[j] / (phi + 1) ) * zj %*% t(zj)
}








