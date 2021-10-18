aa <- test[test$id == 1, ]
# Data objects
Y <- c(aa$Y.1, aa$Y.2)
X <- cbind(1, aa$time, aa$cont, aa$bin)
Z <- cbind(1, aa$time)
XX <- as.matrix(Matrix::bdiag(X, X))
ZZ <- as.matrix(Matrix::bdiag(Z, Z))
# b and beta
b <- MASS::mvrnorm(mu = rep(0, 4), Sigma = D)
betac <- c(beta[1, ], beta[2, ])

# Survival stuff
eta <- c(0,1)
gamma <- c(-1, 1)
source('~/Documents/PhD/Bernhardt/DataFunctions/survFns.R')
ph <- coxph(Surv(surv.time, status) ~ cont + bin, 
            dplyr::distinct(test, id, cont, bin, surv.time, status))
test$survtime <- test$surv.time
sv <- surv.mod(ph, test)
Fi <- sv$Fi; Di <- sv$Di
Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u
K <- as.matrix(cbind(aa[1, ]$cont, aa[1, ]$bin))


# Likelihood --------------------------------------------------------------

# Poisson LL for subject i
sum(-lfactorial(Y)) + sum(-exp(XX %*% betac + ZZ %*% b)) + Y %*% (XX %*% betac + ZZ %*% b)


# RE LL for subject i
-4/2 * log(2 * pi) - 1/2 * log(det(D)) -1/2 * crossprod(b, solve(D) %*% b) # 1x1

# Survival LL for subject i
Di[1] * log(l0i[1]) + Di[1] * (K %*% eta + c(Fi[1, ], Fi[1, ]) %*% (rep(gamma, each = 2) * b)) - 
  l0u[[1]] %*% (exp(K %*% eta) %x% exp(cbind(Fu[[1]], Fu[[1]]) %*% (rep(gamma, each = 2) * b)))

# C++
sourceCpp('../GLM/Poisson/ll.cpp')
# args(ll) -->
# b, Y, lfactY, X, Z, D, K, Delta, l0i, Fi, l0u, Fu, beta, eta, gr, rvFi, nK, q
ll(b, Y, lfactorial(Y), XX, ZZ, D, K, Di[1], l0i[1], Fi[1,], l0u[[1]], Fu[[1]], 
   betac, eta, rep(gamma, each = 2), c(Fi[1, ], Fi[1, ]), 2, 4)


# Gradient terms ----------------------------------------------------------
-crossprod(ZZ, exp(XX %*% betac + ZZ %*% b)) + crossprod(ZZ, Y) - solve(D) %*% b + 
Di[1] * c(Fi[1, ], Fi[1, ]) * rep(gamma, each = 2)
crossprod(cbind(Fu[[1]], Fu[[1]]),
          l0u[[1]] * (exp(K %*% eta) %x% exp(cbind(Fu[[1]], Fu[[1]]) %*% (rep(gamma, each = 2) * b)))
) * rep(gamma, each = 2)

# args(gradll) -->
# b, Y, X, Z, D, K, Delta, l0i, l0u, Fu, g, beta, eta, gr, rvFi, nK
gradll(b, Y, XX, ZZ, D, K, Di[1], l0i[1], l0u[[1]], Fu[[1]], betac,
       eta, rep(gamma, each = 2), c(Fi[1, ], Fi[1, ]), 2)
