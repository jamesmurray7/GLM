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
y %*% log(mu) + crossprod(1-y, log(1-mu))


# GLMM --------------------------------------------------------------------

X <- replicate(100, cbind(1, 0:4, rnorm(1), rbinom(1, 1, .5)), simplify = F); X <- do.call(rbind, X)
Z <- replicate(100, cbind(1, 0:4), simplify = F); Z <- do.call(rbind, Z)

df <- data.frame(
  id = rep(seq_len(100), each = 5),
  X
)

beta <- c(1, 0.10, 0.33, -0.5) # fixed effects coefficients
D11 <- 0.48 # variance of random intercepts
D22 <- 0.1 # variance of random slopes

b <- cbind(rnorm(100, sd = sqrt(D11)), rnorm(100, sd = sqrt(D22)))
# linear predictor
eta <- X %*% beta + rowSums(Z * b[df$id, ])

y <- rbinom(500, 1, plogis(eta))

fit <- glmmTMB::glmmTMB(y ~ X2 + X3 + X4 + (1 + X2|id),
                        data = df, family = binomial)

true.ll <- logLik(fit)
fit.beta <- glmmTMB::fixef(fit)$cond
fit.b <- as.matrix(glmmTMB::ranef(fit)$cond$id)
etahat <- X %*% fit.beta + rowSums(Z * fit.b[df$id,])
mu <- plogis(etahat)
mu <- exp(etahat)/(1+exp(etahat))

y %*% log(mu) + crossprod(1-y, log(1-mu))

ll.fun <- function(b, y, X, beta, Z){
  eta <- X %*% beta + Z %*% b
  mu <- plogis(eta)
  ll <- y %*% log(mu) + crossprod(1 - y, log(1 - mu))
  -ll
}

ll.fun(fit.b[1,], y[1:5], X[df$id==1, ], fit.beta, Z[df$id == 1, ])
nlm(ll.fun, fit.b[1,], y = y[1:5], X = X[df$id==1, ], beta = fit.beta, Z = Z[df$id == 1, ])

ucminf::ucminf(c(fit.b[1,]), ll.fun, NULL, y = y[1:5], X = X[df$id==1, ], beta = fit.beta, Z = Z[df$id == 1, ])
