#' #######
#' Functions for the CMP log-likelihood
#' (as I don't think I can implement it all into one fn!)
#' #######

library(Rcpp)
library(RcppArmadillo)

suppressWarnings(sourceCpp('ll.cpp'))
source('_Functions.R')
source('simData.R')

ll_cmp <- function(X, Y, Z, G, 
                   b, beta, delta, D, summax = 100){
  eta <- X %*% beta + Z %*% b
  mu <- exp(eta)
  nu <- exp(G %*% delta)
  lambda <- mapply(getlambda, mu, nu, summax = summax)
  lY <- lfactorial(Y)
  ll_cmpC(lambda, nu, summax, Y, lY)
}

.ll <- function(b, X, Y, Z, G, beta, delta, D, summax = 100){
  -1 * (ll_cmp(X, Y, Z, G, b, beta, delta, D, summax) + ll_b(b, D))
}

Score_b <- function(b, X, Y, Z, G, beta, delta, D, summax = 100){ # This incorrect currently
  eta <- X %*% beta + Z %*% b
  mu <- exp(eta)
  nu <- exp(G %*% delta)
  lambda <- mapply(getlambda, mu, nu, summax = summax)
  V <- V_mu_lambda(mu, lambda, nu, summax = summax)
  # S <- Smu(mu, lambda, nu, Y, summax = summax)
  Sb <- crossprod(mu*(Y-mu)/V, Z)
  -(Sb + t(S_ll_b(b, D)))
}

