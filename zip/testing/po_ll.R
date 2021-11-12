Y <- rpois(150, 8)
X <- rnorm(150, 4, 2)
glmfit <- glm(Y~X, family = 'poisson')

Xm <- cbind(1,X)
po.logLik <- function(beta){
  ll <- -sum(exp(Xm %*% beta)) + sum(Y %*% (Xm %*% beta)) - sum(lfactorial(Y))
  -ll
}

beta <- nlm(po.logLik, c(0,0))$estimate
beta

library(numDeriv)
po.grad <- function(beta){
  t(-crossprod(Xm, exp(Xm %*% beta))) + Y %*% Xm
}

crossprod(-diag(as.numeric(exp(Xm %*% beta))) %*% Xm, Xm)

sderiv <- function(beta){
  crossprod(-diag(as.numeric(exp(Xm %*% beta))) %*% Xm, Xm)
} #check this with jacobian(po.grad!)

# jacobian(po.grad, beta); sderiv(beta)
 
# While loop NR
beta <- c(0.001,0.001)
iter <- 1
diff <- 100; tol <- 1e-5
while(diff > tol){
  beta.new <- beta + solve(-sderiv(beta), c(po.grad(beta)))
  iter <- iter + 1
  diff <- max(abs(beta-beta.new)/(abs(beta)))
  message("\nIteration ", iter, " Relative difference ", round(diff, 5),
          "\nNew beta = ", paste0(sapply(beta.new, round, 4), collapse =' '))
  beta <- beta.new
}
