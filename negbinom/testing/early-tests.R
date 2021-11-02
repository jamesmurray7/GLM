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

#' ######
#' \beta
#' ######

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

microbenchmark::microbenchmark( # And numderiv appx 10x slower 
  `manual` = {sdbeta(beta, X, Y, theta)},
  `numDeriv` = {numDeriv::hessian(testfn, beta, X = X, Y = Y, theta = theta) }
)


#' ######
#' \theta
#' ######

# First derivative of theta
theta.testfn <- function(theta, X, Y, beta){
  sum(lgamma(theta + Y) - lgamma(theta) - lfactorial(Y) + theta * log(theta) - (theta + Y) * log(exp(X %*% beta) + theta) ) + Y %*% log(exp(X %*% beta))
}

#' 01/11/21 -- Can't work out how to d/dtheta the above. digamma(x)/gamma(x) but where to put summations!
numDeriv::grad(theta.testfn, theta, X = X, Y = Y, beta = beta)
theta + solve(-numDeriv::hessian(theta.testfn, theta, X = X, Y = Y, beta = beta),numDeriv::grad(theta.testfn, theta, X = X, Y = Y, beta = beta))

#'               ~~~~~~~~~~~~~~~~~~~~~~~~~~~                    
#' ##############################################################
#' Repeating for GLMM version
#' ##############################################################

rm(list=ls())
# GLMM --------------------------------------------------------------------
theta <- 1.1
X <- cbind(1, 0:9, rnorm(1), rbinom(1, 1, 0.5))
Z <- X[,c(1,2)]
beta <- c(1, 0.05, 0.5, -0.4)
D <- diag(c(0.5^2, 0.2^2))
b <- MASS::mvrnorm(1, c(0, 0), D)

mu <- exp(X %*% beta + Z %*% b)
Y <- MASS::rnegbin(mu, theta = theta)

void_fn <- function(){
  n <- 150
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  b <- MASS::mvrnorm(n, c(0, 0), D)
  theta <- 1.15
  Y <- list()
  for(j in 1:n){
    X <- cbind(1, 0:9, x1[j], x2[j]); Z <- X[, 1:2]
    mu <- exp(X %*% beta + Z %*% b[j,])
    temp <- temp2 <- numeric(10)
    for(k in seq_along(mu)){
      temp[k] <- MASS::rnegbin(1, mu = mu[k], theta = theta)
    }
    temp2 <- MASS::rnegbin(10, mu = mu, theta = theta)
    Y[[j]] <- data.frame(id = j, X, temp, temp2)
  }
  do.call(rbind, Y)
}
test <- void_fn()


library(glmmTMB)
# Looks as though the temp2 simulation approach (above) better //
# Note that if we set dispersion to be time-varying, we'd need to add ~1+time to dispformula argument
fit <- glmmTMB(temp~X2+X3+X4+(1+X2|id), test, family = nbinom2)
fit2 <- glmmTMB(temp2~X2+X3+X4+(1+X2|id), test, family = nbinom2)

# Is this any different under family = nbinom1? (specifically 'linear')?
fit2.lin <- glmmTMB(temp2~X2+X3+X4+(1+X2|id), test, family = nbinom1)
summary(fit2)
summary(fit2.lin) # Okay this is awful

# Can we speed it up a little without making the model not converge?
fit2.control <- glmmTMB(temp2~X2+X3+X4+(1+X2|id), test, family = nbinom2,
                        control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))

microbenchmark::microbenchmark( # Inclusion of control parameter seems to shave about one second off.
  `with control` = {fit2.control <- glmmTMB(temp2~X2+X3+X4+(1+X2|id), test, family = nbinom2,
                                           control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))},
  `vanilla` = {fit2.van <- glmmTMB(temp2~X2+X3+X4+(1+X2|id), test, family = nbinom2)},
  times = 10
)

#'               ~~~~~~~~~~~~~~~~~~~~~~~~~~~                    
#' ##############################################################
#' Testing simData.R
#' ##############################################################
rm(list=ls())
library(survival)
# simData test ------------------------------------------------------------
source('Documents/GLMM/negbinom/simData.R')
source('Documents/GLMM/negbinom/inits/inits.R')
source('Documents/GLMM/negbinom/inits/MVLME.R')
source('Documents/GLMM/DataFunctions/longFns.R')

args(simData)
beta <- rbind(c(0.75, -0.05, -0.2, 0.2),
              c(1, 0.05,  0.5, -0.80))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5

gamma <- c(-0.5, 0.8)
eta <- c(-0.3, 0.5)

data <- simData(250, 15, beta, D, gamma, eta, theta = c(-5, 0.2)) # appx. 50%
# Check all ids actually here, if not just rerun line above
length(unique(data$id))
ph <- coxph(Surv(survtime, status) ~ cont + bin, data = dplyr::distinct(data, id, cont, bin, survtime, status))
# Looks alright.

X <- getXi(data, 2); Y <- getYi(data, 2); Z <- getZi(data, 2)
Xk <- splitXks(data, 2); Zk <- splitZks(data, 2); mi <- getmi(data, 2); Yk <- splitYks(Y, mi, 2)

inits.long <- Longit.inits(2, data)
mvlme.step <- mvlme(data, Y, X, Z, Yk, Xk, Zk, mi, inits.long, 2, 4, verbose = T)


# Derivatives wrt b -------------------------------------------------------
b0 <- as.matrix(Ranefs(inits.long)[,1:4])
beta=c(0.75, -0.05, -0.2, 0.2)
testfn <- function(b, X, Y, Z, beta, theta){
  sum(lgamma(theta + Y) - lgamma(theta) - lfactorial(Y) + theta * log(theta) - (theta + Y) * log(exp(X %*% beta + Z %*% b) + theta)) + Y %*% log(exp(X %*% beta + Z %*% b))
}

testfn(b = b0[1,1:2], X = X, Y = Y, Z = Z, beta=c(0.75, -0.05, -0.2, 0.2), theta=inits.long$theta.init[1])

b_ll_grad <- function(b, X, Y, Z, beta, theta){
  crossprod(Z, Y) - crossprod(Z, (theta + Y) * exp(X %*% beta + Z %*% b) / (exp(X %*% beta + Z %*% b) + theta)) - solve(D[1:2,1:2]) %*% b
}

b_ll_grad(b, X, Y, Z, beta, theta)  # And these match, which is good!
numDeriv::grad(testfn, b, X = X, Y=Y, Z=Z, beta=beta, theta=theta)

# And second derivative
b_ll_inf <- function(b, X, Y, Z, beta, theta){
  t0 <- exp(X %*% beta + Z %*% b)
  t1 <- (Y + theta) * t0
  t2 <- t0 + theta
  lhs <- crossprod(diag(x=c(t1/t2)) %*% Z, Z)
  rhs <- crossprod(diag(x=c(t1*t0/(t2*t2))) %*% Z, Z)
  -1 * (lhs-rhs)
}

b_ll_inf(b, X, Y, Z, beta, theta)  # And these match too!
numDeriv::hessian(testfn, b, X = X, Y=Y, Z=Z, beta=beta, theta=theta)


# Now check c++ implementation here
library(Rcpp)
library(RcppArmadillo)
sourceCpp('~/Documents/GLMM/negbinom/ll.cpp')
nb_ll(b, Y, lfactorial(Y), X, Z, D[1:2, 1:2], beta, theta, 1)
nb_grad(b, Y, lfactorial(Y), X, Z, D[1:2, 1:2], beta, theta, 1)
nb_hess(b, X, Z, Y, D[1:2,1:2], beta, theta)
