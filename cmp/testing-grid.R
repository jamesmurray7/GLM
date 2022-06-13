rm(list=ls())
library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
source('ll.R')
source('survFns.R')
sourceCpp('grid-test.cpp')

# Load parameter matrices
load('data.nosync/lambda.RData')
load('data.nosync/V.RData')
load('data.nosync/logZ.RData')

# Data and data matrices
test <- simData_joint(n = 250, delta = c(-0.5, 0.1), ntms = 15, theta = c(-3, .25))
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
lY <- lapply(Y, lfactorial)
mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = bl, SIMPLIFY = F)
nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)

# See plot of min and max values
plot(do.call(rbind, lapply(mus, range)), pch = 20,xlab='min',ylab='max'); abline(h=10,v=10,lty=5,col='red')

# Optim for b.hat
# bhat
b.hat <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
  optim(b, joint_density, joint_density_ddb,
        X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
        S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
        gamma = gamma, zeta = zeta, 
        lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, method = 'BFGS')$par
}, b = bl, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)

Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
  solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
                          S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                          gamma = gamma, zeta = zeta, lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, eps = 1e-3))
}, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)

# Check \Sigma for any non-pos-definite-ness 
which(sapply(Sigma, det) < 0)
which(apply(sapply(Sigma, diag), 1, function(x) any(x<0)), arr.ind = T)


# Update for \beta --------------------------------------------------------
mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b.hat, SIMPLIFY = F)
mus2 <- lapply(mus, mu_fix)
lambdas <- mapply(function(mu, nu){
  m <- (mu/0.01) - 0.01; n <- (nu/0.01) - 0.01
  getlambda(m, n, lambda.mat)
}, mu = mus2, nu = nus, SIMPLIFY = F)

Vs <- mapply(function(mu, nu){
  m <- (mu/0.01) - 0.01; n <- (nu/0.01) - 0.01
  getV(m, n, V.mat)
}, mu = mus2, nu = nus, SIMPLIFY = F)

logZ <- mapply(function(mu, nu){
  m <- (mu/0.01) - 0.01; n <- (nu/0.01) - 0.01
  getV(m, n, logZ.mat)
}, mu = mus2, nu = nus, SIMPLIFY = F)


Sb <- mapply(Sbeta, X, Y, mus2, nus, lambdas, Vs, SIMPLIFY = F)
Hb <- mapply(getW1, X, mus, nus, lambdas, Vs, SIMPLIFY = F)
beta - solve(Reduce('+', Hb), Reduce('+', Sb))

# Update for \delta -------------------------------------------------------
E.lfactorialY <- function(lambda, nu, Z, summax){ # mu, nu, vectors
  out <- matrix(0, nr = length(lambda), nc = summax)
  for(j in 1:summax){
    out[, j] <- lgamma(j) * exp((j-1) * log(lambda) - nu * lgamma(j) - Z)
  }
  apply(out, 1, sum)
}
E.YlfactorialY <- function(lambda, nu, Z, summax){
  out <- matrix(0, nr = length(lambda), nc = summax)
  for(j in 1:summax){
    out[, j] <- exp(
      log(j - 1) + log(lgamma(j)) + (j - 1) * log(lambda) - nu * lgamma(j) - Z
    )
  }
  apply(out, 1, sum)
}
V.lfactorialY <- function(lambda, nu, Z, summax, B){
  out <- matrix(0, nr = length(lambda), nc = summax)
  for(j in 1:summax){
    out[, j] <- lgamma(j)^2 * exp((j-1)*log(lambda) - nu * lgamma(j) - Z)
  }
  apply(out, 1, sum) - B^2
}

calc.ABC <- function(mu, nu, lambda, Z, summax){ # NB: Z is log(Z)...
  # lambda <- lambda_uniroot_wrap(1e-6, 1e3, mu, nu, summax)
  B <- E.lfactorialY(lambda, nu, Z, summax)
  A <- E.YlfactorialY(lambda, nu, Z, summax) - mu * B
  C <- V.lfactorialY(lambda, nu, Z, summax, B) # c is potentially needed in W2 matrix creation, remove if not!
  list(A = A, B = B, C = C)
}

ABC <- mapply(calc.ABC, mus2, nus, lambdas, logZ, 100, SIMPLIFY=F)

Sd <- mapply(function(ABC, Y, mu, V, nu, G){
  crossprod(((ABC$A * (Y - mu) / V - lgamma(Y + 1) + ABC$B) * nu), G)
}, ABC = ABC, Y = Y, mu = mus2, V = Vs, nu = nus, G = G, SIMPLIFY = F)
Hd <- mapply(getW2, ABC, Vs, nus, G, SIMPLIFY = F)

delta-solve(Reduce('+', Hd), c(Reduce('+', Sd)))
