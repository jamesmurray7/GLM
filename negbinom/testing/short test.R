beta <- c(1, 0.05, 0.5, -0.4)
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

fit <- void_fn()
logLik(fit$fit)

# Get something to match
beta <- fit[[2]]
X <- fit[[4]]
Y <- fit[[3]]
mu <- exp(X %*% beta)
theta <- fit[[1]]$theta
sum(lgamma(theta + Y) - lgamma(theta) - lfactorial(Y) + theta * log(theta) - (theta + Y) * log(mu + theta) ) + Y %*% log(mu)
sum(dnbinom(Y, 1, mu = exp(X %*% beta)*log(theta), log = T)) # A little bit different?

# what are phis in glmmA?
sum(GLMMadaptive::negative.binomial()$log_dens(y = Y, eta = X %*% beta, mu_fun = exp, phis = theta, eta_zi = 0))
sum(GLMMadaptive::negative.binomial()$log_dens(y = Y, eta = X %*% beta, mu_fun = exp, phis = log(theta), eta_zi = 0))
# So GLMMadaptive is parameterised as phi = log(theta) -> exp(phi) = theta.

# Score on eta
crossprod(X, GLMMadaptive::negative.binomial()$score_eta(Y, exp(X %*% beta), log(theta), 0))
crossprod(X, theta * (Y-exp(eta))/(theta+exp(eta)))

# second derivative
scoreBeta <- function(beta, Y, X, theta){
  eta <- X %*% beta
  crossprod(X, theta * (Y-exp(eta))/(theta+exp(eta)))
}

Hb <- GLMMadaptive:::fd_vec(beta, scoreBeta, Y = Y, X = X, theta = theta)
beta-solve(Hb,scoreBeta(beta, Y, X, theta))

# Score on theta
eta <- X %*% beta
St <- digamma(theta + Y) - digamma(theta) + log(theta) + 1 - ((theta+Y)/(exp(eta) + theta) + log(exp(eta) + theta))

scoreTheta <- function(theta, Y, X, beta){
  eta <- X %*% beta
  St <- digamma(theta + Y) - digamma(theta) + log(theta) + 1 - ((theta+Y)/(exp(eta) + theta) + log(exp(eta) + theta))
  sum(St) 
}

scoreTheta(theta, Y, X, beta)
sum(GLMMadaptive::negative.binomial()$score_phis(Y, exp(eta), log(theta), 0)) # GLMMadaptive is post-multiplied by theta(?)

Ht <- GLMMadaptive:::fd(theta, scoreTheta, Y = Y, X = X, beta = beta)
theta-scoreTheta(theta, Y, X, beta)/Ht

    