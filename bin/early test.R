X <- replicate(1000, cbind(1, rnorm(1)), simplify = F)
X <- do.call(rbind, X)
beta <- c(1, 0.4)

eta <- X %*% beta
y <- rbinom(1000, 1, plogis(eta))

fit <- glm(y ~ X[,2], family = binomial)
true.ll <- logLik(fit)

fit.beta <- fit$coefficients
etahat <- X %*% fit.beta
mu <- plogis(etahat)
y %*% log(mu) + crossprod(1-y, log(1-mu)) ## The same ----
sum(dbinom(y, 1, plogis(etahat),T))

# GLMM --------------------------------------------------------------------
n <- 250
ntms <- 10
beta <- c(1, 0.10, 0.33, -0.5) # fixed effects coefficients
D <- matrix(c(0.5, 0, 0, 0.1), 2, 2)
b <- MASS::mvrnorm(n, mu=c(0, 0), Sigma = matrix(c(0.5, 0, 0, 0.1), 2, 2))

df <- data.frame(
  id = rep(1:n, each = ntms),
  time = rep(0:(ntms-1), n),
  cont = rep(rnorm(n), each = ntms),
  bin = rep(rbinom(n, 1, 0.5), each = ntms)
)

X <- model.matrix(~time+cont+bin, df)
Z <- model.matrix(~time, df)
# Linear predictor
eta <- X %*% beta + rowSums(Z * b[df$id, ])
df$Y <- rbinom(n * ntms, 1, plogis(eta))

fit <- glmmTMB::glmmTMB(Y ~ time + cont + bin + (1 + time|id),
                        data = df, family = binomial)

true.ll <- logLik(fit)
fit.beta <- glmmTMB::fixef(fit)$cond
fit.b <- as.matrix(glmmTMB::ranef(fit)$cond$id)
etahat <- X %*% fit.beta + rowSums(Z * fit.b[df$id,])
mu <- plogis(etahat)

Y %*% log(mu) + crossprod(1-Y, log(1-mu))
sum(dbinom(Y, 1, mu, TRUE))

# is it missing f(b)?
D <- as.matrix(glmmTMB::VarCorr(fit)$cond$id)
Dll <- -1/2 * log(2*pi) - 1/2 * log(det(D)) - 1/2 * sum(diag(fit.b %*% solve(D) %*% t(fit.b)))

true.ll
sum(dbinom(Y, 1, mu, TRUE)) + Dll # closer
