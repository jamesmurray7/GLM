#' ####
#' early tests for negbinom ll
#' ###

rm(list=ls())

theta <- 1.1
X <- cbind(1, 0:9, rnorm(1), rbinom(1, 1, 0.5))
Z <- X[,c(1,2)]
beta <- c(1, 0.05, 0.5, -0.4)
D <- diag(c(0.5^2, 0.2^2))
b <- MASS::mvrnorm(1, c(0, 0), D)

mu <- exp(X %*% beta + Z %*% b)
mu

Y <- MASS::rnegbin(mu, theta = theta)


MASS::glm.nb(Y~X[,2]+X[,3]+X[,4])


# Simulating using MASS::rnegbin ------------------------------------------

void_fn <- function(){
  n <- 150
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  theta <- 1.5
  Y <- c()
  for(j in 1:n){
    Y[j] <- MASS::rnegbin(1, mu = exp(1 * beta[1] + x1[j] * beta[2] + x2[j] * beta[3]), theta = 1.5)
  }
  fit <- MASS::glm.nb(Y~x1+x2)
  return(list(fit = fit, beta_ests = coef(fit), response = Y, data = cbind(1, x1, x2)))
}

fits <- replicate(100, void_fn(), simplify = F)  
colMeans(do.call(rbind, lapply(fits, '[[', 2)))     ## This seems to be good for beta = c(1, 0.05, 0.5).
do.call(c, lapply(lapply(fits, '[[', 1), logLik))   ## About -350

# Get something to match
logLik(fits[[1]][[1]])
beta <- fits[[1]][[2]]
X <- fits[[1]][[4]]
Y <- fits[[1]][[3]]
mu <- exp(X %*% beta)
theta <- fits[[1]][[1]]$theta
sum(lgamma(theta + Y) - lgamma(theta) - lfactorial(Y) + theta * log(theta) - (theta + Y) * log(mu + theta) ) + Y %*% log(mu)


# Derivatives -------------------------------------------------------------
# \beta ---
# First derivative for beta
crossprod(X, Y) - crossprod(X, (Y + theta) * exp(X %*% beta) / (exp(X %*% beta) + theta))

# Check using numDeriv::grad
testfn <- function(beta, X, Y, theta){
  sum(lgamma(theta + Y) - lgamma(theta) - lfactorial(Y) + theta * log(theta) - (theta + Y) * log(exp(X %*% beta) + theta) ) + Y %*% log(exp(X %*% beta))
}

testfn(beta, X, Y, theta) # same ll

# Compare
numDeriv::grad(testfn, beta, X = X, Y = Y, theta = theta) 
crossprod(X, Y) - crossprod(X, (Y + theta) * exp(X %*% beta) / (exp(X %*% beta) + theta))

# Second derivative for \beta
sdbeta <- function(beta, X, Y, theta){
  t0 <- exp(X %*% beta)
  t1 <- (Y + theta) * t0
  t2 <- t0 + theta
  -1 * (
    crossprod(diag(x = c(t1 / t2)) %*% X, X) - crossprod(diag(x = c(t1 * t0  / (t2 * t2))) %*% X, X)
  )
}

# Then these are approximately the same!
sdbeta(beta, X, Y, theta)
numDeriv::hessian(testfn, beta, X = X, Y = Y, theta = theta) 

microbenchmark::microbenchmark(
  `manual` = {sdbeta(beta, X, Y, theta)},
  `numDeriv` = {numDeriv::hessian(testfn, beta, X = X, Y = Y, theta = theta) }
)


# \theta ---
# First derivative of theta
sum(
  digamma(theta+Y)/gamma(theta+Y)- digamma(theta)/gamma(theta) + 1 + log(theta) - ((theta + Y)/(exp(X %*% beta) + theta) + log(exp(X %*% beta + theta) + theta))
)

sum(digamma(theta+Y)/gamma(theta+Y)) - digamma(theta)/gamma(theta) + 1 + log(theta) - sum(((theta + Y)/(exp(X %*% beta) + theta) + log(exp(X %*% beta + theta) + theta)))


test_theta_fn <- function(theta, X, Y, beta){
  sum(lgamma(theta + Y) - lgamma(theta) - lfactorial(Y) + theta * log(theta) - (theta + Y) * log(exp(X %*% beta) + theta) ) + Y %*% log(exp(X %*% beta))
}

numDeriv::grad(test_theta_fn, theta, X = X, Y = Y, beta = beta)
