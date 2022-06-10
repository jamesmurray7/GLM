#' #######
#' testing.R
#' Run-through of E- and M- steps 
#' #######
rm(list=ls())
library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
source('ll.R')
source('survFns.R')
sourceCpp('cmp.cpp')

# Data and data matrices
test <- simData_joint()
data <- test$data; survdata <- test$surv.data
n <- length(unique(data$id))
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)
sv <- surv.mod(ph, survdata, NULL)
S <- X <- Z <- G <- Y <- list()
for(i in 1:n){
  X[[i]] <- model.matrix(~ time + cont + bin, data[data$id == i, ])
  Z[[i]] <- G[[i]] <- model.matrix(~ time, data[data$id == i, ])
  Y[[i]] <- data[data$id == i, 'Y']
  S[[i]] <- t(unname(cbind(data[data$id==i, 'cont'], data[data$id==i, 'bin']))[1, ])
}
Fu <- sv$Fu; Fi <- sv$Fi; l0u <- sv$l0u; l0i <- as.list(sv$l0i); Delta <- as.list(sv$Di)
Fi <- lapply(1:n, function(i) Fi[i,,drop = F])
SS <- lapply(1:n, function(i){
  x <- apply(S[[i]],2,rep,nrow(Fu[[i]]))
  if(any(class(x)=='numeric'))x <- t(x)
  x
})
gamma <- 0.4; zeta <- c(0.1, -0.2)
# Initial conditions 
fitP <- glmmTMB(Y ~ time + cont + bin + (1 + time|id),
               data, family = poisson)
beta <- glmmTMB::fixef(fitP)$cond
delta <- delta <- c(0,0) # glmmTMB::fixef(fit)$disp # delta <- c(1, 0) ## True value -0.1, 0.2
b <- glmmTMB::ranef(fitP)$cond$id
bl <- lapply(1:n, function(i) as.numeric(b[i,]))
D <- matrix(glmmTMB::VarCorr(fitP)$cond$id,2,2)

mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = bl, SIMPLIFY = F)
nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)

lambdas <- mapply(function(mu, nu) lambda_uniroot_wrap(1e-6, 1e3, mu, nu, 10), mu = mus, nu = nus, SIMPLIFY = F)
sapply(1:10, function(i) uniroot(mu_lambdaZ_eq2, interval=c(1e-6, 1e3), mus[[2]][i], nus[[2]][i], summax = 100)$root)
lY <- lapply(Y, lfactorial)

# bhat
b.hat <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
  optim(b, joint_density, joint_density_ddb,
        X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
        S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
        gamma = gamma, zeta = zeta, summax = 100, method = 'BFGS')$par
}, b = bl, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)

Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
  solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
        S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
        gamma = gamma, zeta = zeta, summax = 100, eps = 1e-3))
}, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)

apply(do.call(rbind, lapply(Sigma, diag)),2,function(x)any(x<0)) # Check, setting summax ~20 seems to alleviate issue here.
                                                                 # unsure as to why this is the case(!)
which(do.call(rbind, lapply(Sigma, diag))<0,arr.ind=T)           # For checking which id is going weird.
which(sapply(Sigma, det) < 0)                                    # Checking non- semi-pos-def'ness

# \beta -------------------------------------------------------------------

Sbeta <- function(beta, X, Y, Z, G, b, delta, summax){
  mu <- exp(X %*% beta + Z %*% b)
  nu <- exp(G %*% delta)
  lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
  V <- V_mu_lambda(mu, lambda, nu, summax)
  if(length(mu) > 1) lhs <- diag(c(mu)) else lhs <- diag(mu)
  crossprod(lhs %*% X, ((Y-mu) / V))
}

Sb <- mapply(function(X, Y, Z, G, b){
  Sbeta(beta, X, Y, Z, G, b, delta, summax = 100)
}, X = X, Y = Y, Z = Z, G = G, b = b.hat, SIMPLIFY = F)

Hb <- mapply(function(X, Y, Z, G, b){
  GLMMadaptive:::fd_vec(beta, Sbeta, X, Y, Z, G, b, delta, summax = 100)
}, X = X, Y = Y, Z = Z, G = G, b = b.hat, SIMPLIFY = F)

# Huang (2016) W1 asymp. variance of \beta
W1 <- mapply(function(X, Z, G, b){
  getW1(X, Z, G, b, beta, delta, summax = 100)
}, X = X, G = G, Z = Z, b = b.hat, SIMPLIFY = F)

Hb[[1]]
W1[[1]]
beta-solve(Reduce('+', Hb), Reduce('+', Sb))
beta-solve(Reduce('+', W1), Reduce('+', Sb))  # Looks good no matter what approach - likely a matter of benchmarking & choosing.

# \delta ------------------------------------------------------------------
Sdelta <- function(delta, X, Y, lY, Z, b, G, beta, summax){
  mu <- exp(X %*% beta + Z %*% b)
  nu <- exp(G %*% delta)
  lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
  V <- V_mu_lambda(mu, lambda, nu, summax)
  ABC <- calc.ABC(mu, nu, lambda, summax)
  Snu <- ABC$A * (Y-mu) / V - (lY-ABC$B)
  if(length(nu) > 1) lhs <- diag(c(nu)) else lhs <- diag(nu)
  crossprod(lhs %*% G, Snu)
}

Sd <- mapply(function(X, Y, lY, Z, G, b){
  Sdelta(delta, X, Y, lY, Z, b, G, beta, summax = 100)
}, X = X, Y = Y, lY=lY, Z = Z, G = G, b = b.hat, SIMPLIFY = F)

Hd <- mapply(function(X, Y, lY, Z, G, b){
  GLMMadaptive:::fd_vec(delta, Sdelta, X = X, Y = Y, lY = lY, Z = Z, b = b, G = G, beta = beta, summax = 100)
}, X = X, Y = Y, lY=lY, Z = Z, G = G, b = b.hat, SIMPLIFY = F)

delta - solve(Reduce('+', Hd), Reduce('+', Sd))


# Testing for W2 from Huang paper ...
# Get all ABC vectors
mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b.hat, SIMPLIFY = F)
nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)
lambdas <- mapply(function(mu, nu) lambda_uniroot_wrap(1e-6, 1e3, mu, nu, 10), mu = mus, nu = nus, SIMPLIFY = F)

ABCs <- mapply(calc.ABC, mu = mus, nu = nus, lambda = lambdas, summax = 100, SIMPLIFY = F)
Vs <- mapply(V_mu_lambda, mu = mus, lambda = lambdas, nu = nus, summax = 100, SIMPLIFY = F)

Sd2 <- mapply(function(ABC, Y, mu, V, nu, G){
  crossprod(((ABC$A * (Y - mu) / V - lgamma(Y + 1) + ABC$B) * nu), G)
}, ABC = ABCs, Y = Y, mu = mus, V = Vs, nu = nus, G = G) ## VERY similar to Sd as calculated above! A matter of benchmarking ...

Hd2 <- mapply(function(ABC, nu, V, G){
  W2 <- matrix(0, ncol(G), ncol(G))
  for(j in 1:nrow(G)){
    W2 <- W2 + ((-(ABC$A[j])^2/V[j]+ABC$C[j])*nu[j]^2)*G[j,] %*% t(G[j,])
  }
  -W2
}, ABC = ABCs, nu = nus, V = Vs, G = G, SIMPLIFY = F)

Hd[[1]]
Hd2[[1]]


delta - solve(Reduce('+', Hd), Reduce('+', Sd))
delta - solve(Reduce('+', Hd2), rowSums(Sd2))

GH <- statmod::gauss.quad.prob(3,'normal')
w <- GH$w;v <- GH$n
tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), Z = Z, S = Sigma, SIMPLIFY = F)
ABCs2 <- mapply(calc2.ABC, mus, nus, lambdas, 100, tau, w, v, SIMPLIFY = F)
calc2.ABC(mus[[1]], nus[[1]], lambdas[[1]], 100, tau[[1]], w, v)

Sd2q <- mapply(function(ABC, Y, mu, V, nu, G, tau){
  lhs <- numeric(length(mu))
  for(l in 1:length(w)) lhs <- lhs + w[l] * ABC$A * (Y - mu * exp(tau * v[l])) / V
  crossprod(((lhs - lgamma(Y + 1) + ABC$B) * nu), G)
}, ABC = ABCs, Y = Y, mu = mus, V = Vs, nu = nus, G = G, tau= tau, SIMPLIFY = F) ## VERY similar to Sd as calculated above! A matter of benchmarking ...

#'------------------------------------------
#'  OLD VERSION BELOW ----------------------
#'------------------------------------------
#'------------------------------------------


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
