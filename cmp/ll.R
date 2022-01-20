#' #######
#' Functions for the CMP log-likelihood
#' (as I don't think I can implement it all into one fn!)
#' #######

suppressWarnings(sourceCpp('ll.cpp'))
source('_Functions.R')

ll_cmp <- function(X, Y, Z, G, 
                   b, beta, delta, D, summax = 100){
  eta <- X %*% beta + Z %*% b
  mu <- exp(eta)
  nu <- exp(G %*% delta)
  lambda <- mapply(getlambda, mu, nu, summax = summax)
  lY <- lfactorial(Y)
  ll_cmpC(lambda, nu, summax, Y, lY)
}

test <- simData_joint()
data <- test$data
n <- length(unique(data$id))
X <- Z <- G <- Y <- list()
for(i in 1:n){
  X[[i]] <- model.matrix(~ time + cont + bin, data[data$id == i, ])
  Z[[i]] <- G[[i]] <- model.matrix(~ time, data[data$id == i, ])
  Y[[i]] <- data[data$id == i, 'Y']
}

.ll <- function(b, X, Y, Z, G, beta, delta, D, summax = 100){
  -1 * (ll_cmp(X, Y, Z, G, b, beta, delta, D, summax) + ll_b(b, D))
}

.bfits <- mapply(function(b, X, Y, Z, G){
  u <- ucminf::ucminf(b, .ll, NULL,
                      X, Y, Z, G, beta, delta, D, summax = 10,
                      control = list(xtol = 1e-4, grtol = 1e-6))
  list(u$par, u$invhessian.lt)
}, b = bl, X = X, Y = Y, Z = Z, G = G, SIMPLIFY = F)

b.hat <- lapply(.bfits, el, 1)
Sigmai <- lapply(lapply(.bfits, el, 2), vech2mat, 2)

tau.long <- mapply(function(Z, S){
  unname(sqrt(diag(tcrossprod(Z %*% S, Z))))
}, Z = Z, S = Sigmai, SIMPLIFY = F)

