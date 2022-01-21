#' #######
#' Functions for the CMP log-likelihood
#' (as I don't think I can implement it all into one fn!)
#' #######

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


# Data and data matrices
test <- simData_joint()
data <- test$data
n <- length(unique(data$id))
X <- Z <- G <- Y <- list()
for(i in 1:n){
  X[[i]] <- model.matrix(~ time + cont + bin, data[data$id == i, ])
  Z[[i]] <- G[[i]] <- model.matrix(~ time, data[data$id == i, ])
  Y[[i]] <- data[data$id == i, 'Y']
}

# Initial conditions (Using negbinom2 until can think of something better!)
fit <- glmmTMB::glmmTMB(Y ~ time + cont + bin + (1 + time|id), 
                        dispformula = ~ time,
                        data, 
                        family = glmmTMB::nbinom2)
beta <- glmmTMB::fixef(fit)$cond
delta <- glmmTMB::fixef(fit)$disp
b <- glmmTMB::ranef(fit)$cond$id
bl <- lapply(1:n, function(i) matrix(b[i,], nr = 1))
D <- matrix(glmmTMB::VarCorr(fit)$cond$id,2,2)

.ll <- function(b, X, Y, Z, G, beta, delta, D, summax = 100){
  -1 * (ll_cmp(X, Y, Z, G, b, beta, delta, D, summax) + ll_b(b, D))
}

Score_b <- function(X, Y, Z, G, 
                    b, beta, delta, D, summax = 100){
  eta <- X %*% beta + Z %*% b
  mu <- exp(eta)
  nu <- exp(G %*% delta)
  lambda <- mapply(getlambda, mu, nu, summax = summax)
  S <- Smu(mu, lambda, nu, Y, summax = summax)
  -(crossprod(Z, S) + S_ll_b(b, D))
}


# b.hat ----
.bfits <- mapply(function(b, X, Y, Z, G){
  u <- ucminf::ucminf(b, .ll, Score_b,
                      X = X, Y = Y, Z = Z, G = G, beta = beta, delta = delta, D = D, summax = 10,
                      control = list(xtol = 1e-4, grtol = 1e-6))
  list(u$par, u$invhessian.lt)
}, b = bl, X = X, Y = Y, Z = Z, G = G, SIMPLIFY = F)

b.hat <- lapply(.bfits2, el, 1)
Sigmai <- lapply(lapply(.bfits2, el, 2), vech2mat, 2)


# D update ----
Drhs <- mapply(function(b, S){
  S + tcrossprod(b)
}, b = b.hat, S = Sigmai, SIMPLIFY = F)

Reduce('+', Drhs)/n

# beta update ----

tau.long <- mapply(function(S, Z){
  unname(sqrt(diag(tcrossprod(Z %*% S, Z))))
}, S = Sigmai, Z = Z, SIMPLIFY = F)

.ngh <- 3
gh <- statmod::gauss.quad.prob(.ngh, 'normal')
w <- gh$w; v <- gh$n

.Sbeta <- function(beta, b, X, Y, Z, G, delta, summax = 10){
  eta <- X %*% beta + Z %*% b
  mu <- exp(eta)
  nu <- exp(G %*% delta)
  lambda <- mapply(getlambda, mu, nu, summax = summax)
  S <- Smu(mu, lambda, nu, Y, summax = summax)
  crossprod(X, S)
}

Sbeta <- mapply(function(b, X, Y, Z, G){
  out <- vector('list', 2)
  out[[1]] <- .Sbeta(beta, b, X, Y, Z, G, delta)
  out[[2]] <- GLMMadaptive:::fd_vec(beta, .Sbeta, b = b, X = X, Y = Y, Z = Z, G = G, delta = delta)
  out
}, b = b.hat, X = X, Y = Y, Z = Z, G = G, SIMPLIFY = F)

Sb <- lapply(Sbeta, el, 1)
Hb <- lapply(Sbeta, el, 2)
beta - solve(Reduce('+', Hb), Reduce('+', Sb))
