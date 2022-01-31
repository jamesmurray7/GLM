#' #######
#' testing.R
#' Run-through of E- and M- steps 
#' #######

source('ll.R')

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
delta <- glmmTMB::fixef(fit)$disp # delta <- c(1, 0)
b <- glmmTMB::ranef(fit)$cond$id
bl <- lapply(1:n, function(i) as.numeric(b[i,]))
D <- matrix(glmmTMB::VarCorr(fit)$cond$id,2,2)

# b.hat ----
.bfits <- mapply(function(b, X, Y, Z, G){
  u <- ucminf::ucminf(b, .ll, Score_b,  
                      X = X, Y = Y, Z = Z, G = G, beta = beta, delta = delta, D = D, summax = 10)
  list(u$par, u$invhessian.lt)
}, b = bl, X = X, Y = Y, Z = Z, G = G, SIMPLIFY = F)

# Posterior mode of b ~ N(bhat, Sigma)
b.hat <- lapply(.bfits, el, 1)
# And posterior variance 
Sigmai <- mapply(function(b, X, Y, Z, G){
  solve(pracma::hessian(.ll, b, X = X, Y = Y, Z = Z, G = G, beta = beta, delta = delta, D = D, summax = 10))
}, b = b.hat, X = X, Y = Y, Z = Z, G = G, SIMPLIFY = F)

Sigmai <- lapply(Sigmai, function(x){
  if(det(x) <= 0) return(as.matrix(Matrix::nearPD(x)$mat))
  else return(x)
})

# objects for mu and nu
mu <- mapply(function(X, Z, b){
  exp(X %*% beta + Z %*% b)
}, X = X, Z = Z, b = b.hat, SIMPLIFY = F)

nu <- mapply(function(G){
  exp(G %*% delta)
}, G = G, SIMPLIFY = F)

# lambda, from uniroot on mu and nu
lambda <- mapply(getlambda, mu, nu, summax = 50, SIMPLIFY = F)
log.lambda <- lapply(lambda, log)
# Z, normalising constant using lambda and nu
logZ_ <- mapply(logZ_c, log.lambda, nu, summax = 10, SIMPLIFY = F)
Z_ <- lapply(logZ_, exp)

lfactY <- mapply(E_logfactY, lambda, nu, logZ_, summax = 10, SIMPLIFY = F)
YlfactY <- mapply(E_YlogfactY, lambda, nu, logZ_, summax = 10, SIMPLIFY = F)
means <- mapply(E_means, lambda, nu, logZ_, summax = 10, SIMPLIFY = F)
V <- mapply(V_mu_lambda, mu, lambda, nu, summax = 10, SIMPLIFY = F)

# D update ----
Drhs <- mapply(function(b, S){
  S + tcrossprod(b)
}, b = b.hat, S = Sigmai, SIMPLIFY = F)

Reduce('+', Drhs)/n


# beta update -------------------------------------------------------------
# No quadrature
Sbeta <- mapply(function(Y, mu, V, X){
  crossprod(((Y - mu) / V) * mu, X)
}, Y = Y, mu = mu, V = V, X = X)

W1 <- mapply(function(mu, X, V){
  nr <- nrow(X); nc <- ncol(X)
  out <- matrix(0, nc, nc)
  for(i in 1:nr){
    out <- out + mu[i,]/V[i,] * tcrossprod(X[i,])
  }
  out
},  mu = mu, X = X, V = V, SIMPLIFY = F)

beta + solve(Reduce('+', W1), rowSums(Sbeta))

# Quadrature
.ngh <- 9
gh <- statmod:::gauss.quad.prob(.ngh, 'normal')
w <- gh$w;v <- gh$n
tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), Z = Z, S = Sigmai, SIMPLIFY = F)
qSbeta <- mapply(function(Y, X, Z, b, nu, tau){
  lhs <- numeric(length(Y))
  for(l in 1:.ngh){
    mu <- exp(X %*% beta + Z %*% b + tau * v[l])
    lam <- getlambda(mu, nu, 10)
    V <- V_mu_lambda(mu, lam, nu, 10)
    lhs <- lhs + w[l] * mu * (Y - mu) / V
  }
  crossprod(lhs, X)  
}, Y = Y, X = X, Z = Z, b = b.hat, nu = nu, tau = tau, SIMPLIFY = F)

# Functionise around \beta for forward differencing
Score_beta_quad <- function(beta, Y, X, Z, b, nu, tau){
  lhs <- numeric(length(Y))
  for(l in 1:.ngh){
    mu <- exp(X %*% beta + Z %*% b + tau * v[l])
    lam <- getlambda(mu, nu, 10)
    V <- V_mu_lambda(mu, lam, nu, 10)
    lhs <- lhs + w[l] * mu * (Y - mu) / V
  }
  crossprod(lhs, X) 
}

Sb <- mapply(Score_beta_quad, beta = replicate(n, beta, simplify = F),
             Y = Y, X = X, Z = Z, b = b.hat, nu = nu, tau = tau, SIMPLIFY = F) # score_beta

H_beta_quad <- mapply(function(Y, X, Z, b, nu, tau){
  GLMMadaptive:::fd_vec(beta, Score_beta_quad, Y = Y, X = X, Z = Z, b = b, nu = nu, tau = tau)
}, Y = Y, X = X, Z = Z, b = b.hat, nu = nu, tau = tau, SIMPLIFY = F)

Hb <- Reduce('+', H_beta_quad)

beta-solve(Hb, colSums(do.call(rbind, Sb)))

#' ========================================================================
#' ========================================================================
#' ========================================================================
#' ========================================================================

# delta update ------------------------------------------------------------
# No quadrature
A <- mapply(function(lY, mu, YlY) YlY - mu * lY, lY = lfactY, mu = mu, YlY = YlfactY, SIMPLIFY = F)
B <- lfactY
C <- mapply(Var_logfactY, lambda, nu, logZ_, summax = 10, SIMPLIFY = F)

Snu <- mapply(function(A, Y, mu, V, B){
  A * (Y - mu) / V - (lgamma(Y + 1) - B)
}, A = A, Y = Y, mu = mu, V = V, B = B, SIMPLIFY = F)

sdel <- mapply(function(Sn, nu, G){
  crossprod(Sn * nu, G)
}, Sn = Snu, nu = nu, G = G)

# or e.g crossprod(crossprod(((A[[1]]/V[[1]]+C[[1]])%*%crossprod(nu[[1]])), G[[1]]))?
W2 <- mapply(function(A, V, C, nu, G){
  lhs <- (-(A)^2/V + C) * nu^2 # m x 1
  crossprod(crossprod(lhs, G))
}, A = A, V = V, C = C, nu = nu, G = G, SIMPLIFY = F)

delta + solve(Reduce('+', W2), rowSums(sdel))

# Using Quadrature
Score_nu <- function(A, Y, mu, V, B) A * (Y - mu) / V - (lgamma(Y + 1) - B)

Score_delta <- function(delta, Y, X, Z, b, G, tau){
  nu <- exp(G %*% delta)
  lhs <- numeric(length(Y))
  for(l in 1:.ngh){
    mu <- exp(X %*% beta + Z %*% b + v[l] * tau)
    lambda <- getlambda(mu, nu, 10)
    log.lambda <- log(lambda)
    logZ_ <- logZ_c(log.lambda, nu, 10)
    lfY <- E_logfactY(lambda, nu, logZ_, 10)
    A <- E_YlogfactY(lambda, nu, logZ_, 10) - lfY * mu
    B <- lfY
    V <- V_mu_lambda(mu, lambda, nu, 10)
    lhs <- lhs + w[l] * Score_nu(A, Y, mu, V, B)
  }
  crossprod(lhs * nu, G)
}

Sd <- mapply(Score_delta, replicate(n, delta, F),
             Y = Y, X = X, Z = Z, b = b.hat, G = G, tau = tau)

Hd <- mapply(function(Y, X, Z, b, G, tau){
  GLMMadaptive:::fd_vec(
    delta, Score_delta, Y = Y, X = X, Z = Z, b = b, G = G, tau = tau
  )
}, Y = Y, X = X, Z = Z, b = b.hat, G = G, tau = tau, SIMPLIFY = F)


delta - solve(Reduce('+', Hd), c(rowSums(Sd)))
