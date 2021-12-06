#' #####
#' Testing functions for joint model negbin test.
#' ####
rm(list=ls())
source('./simData.R')
source('./inits.R')
source('survFns.R')
library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
sourceCpp('nb.cpp')

set.seed(123)
testsim <- simData_joint()
test <- testsim$data
survdata <- testsim$surv.data
fit <- glmmTMB(Y ~ time + cont + bin + (1+time|id), test, dispformula = ~1, family = nbinom2)
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)


n <- 250
Y <- X <- Z <- K <- list()
for(i in 1:n){
  i.dat <- test[test$id == i, ]
  Y[[i]] <- i.dat$Y
  X[[i]] <- model.matrix(~time+cont+bin, i.dat)
  Z[[i]] <- model.matrix(~time, i.dat)
  K[[i]] <- cbind(unique(i.dat$cont), unique(i.dat$bin))
}

inits.long <- Longit.inits(test)
b <- Ranefs(inits.long)
beta <- inits.long$beta.init
theta <- inits.long$theta.init
D <- inits.long$D.init
inits.surv <- TimeVarCox(test, b)
b <- lapply(1:n, function(i) b[i, ])

gamma <- inits.surv$inits[3]
eta <- inits.surv$inits[1:2]

#svmod
sv <- surv.mod(ph, test, inits.surv$l0.init)
Fi <- lapply(1:n, function(i) sv$Fi[i, ])
Fu <- sv$Fu
KK <- sapply(1:n, function(i){
  x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
  if('numeric'%in%class(x)) x <- t(as.matrix(x))
  x
})
l0i <- as.list(sv$l0i); l0u <- sv$l0u
Delta <- as.list(sv$Di)

# Estimates for b and subsequently for Sigmai
b.hat <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, haz){
  ucminf::ucminf(b, joint_density, joint_density_ddb,
                 X = X, Y = Y, Z = Z, beta = beta, theta = theta, D = D,
                 Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu,
                 haz = haz, gamma = gamma, eta = eta)$par
}, b = b, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)

b.hat2 <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, haz){
  ucminf::ucminf(b, joint2_density, joint_density_ddb,
                 X = X, Y = Y, Z = Z, beta = beta, theta = theta, D = D,
                 Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu,
                 haz = haz, gamma = gamma, eta = eta)$par
}, b = b, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)

Sigmai <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, haz){
  solve(joint_density_sdb(b, X, Y, Z, beta, theta, D, Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta, eps = 1e-4))
}, b = b.hat, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)

# Update to D
Drhs <- mapply(function(b, S){
  S + tcrossprod(b)
}, S = Sigmai, b = b.hat, SIMPLIFY = F)

# beta
Sb <- mapply(function(X, Y, Z, b){
  Sbeta(beta, X, Y, Z, b, theta)
}, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)

Hb <- mapply(function(X, Y, Z, b){
  Hbeta(beta, X, Y, Z, b, theta, 1e-4)
}, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)

# theta
St <- mapply(function(X, Y, Z, b){
  Stheta(theta, beta, X, Y, Z, b)
}, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)

Ht <- mapply(function(X, Y, Z, b){
  Htheta(theta, beta, X, Y, Z, b, 1e-4)
}, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)

# gamma and eta
aa <- statmod::gauss.quad.prob(9, 'normal')
w <- aa$w; v <- aa$n
Sge <- mapply(function(b, Delta, Fi, K, KK, Fu, l0u, S){
  Sgammaeta(c(gamma, eta), b, Delta, Fi, K, KK, Fu, l0u, S, w, v)
}, b = b.hat, Delta = Delta, Fi = Fi, K = K, KK = KK, Fu = Fu, l0u = l0u, S = Sigmai)

Hge <- mapply(function(b, Delta, Fi, K, KK, Fu, l0u, S){
  Hgammaeta(c(gamma, eta), b, Delta, Fi, K, KK, Fu, l0u, S, w, v, 1e-4)
}, b = b.hat, Delta = Delta, Fi = Fi, K = K, KK = KK, Fu = Fu, l0u = l0u, S = Sigmai, SIMPLIFY = F)

# updates
D.new <- Reduce('+', Drhs)/250
beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
theta.new <- theta-sum(do.call(c, St))/sum(do.call(c, Ht))
theta2.new <- theta-sum(do.call(c, St2))/sum(do.call(c, Ht2))
gamma.eta.new <- c(gamma, eta) - solve(Reduce('+', Hge), rowSums(Sge))
lambda <- lambdaUpdate(sv$surv.times, sv$ft, gamma, eta, K, Sigmai, b.hat, w, v)

# Baseline hazard
l0 <- sv$nev/rowSums(lambda)
l0.new <- sv$nev/rowSums(lambda)
l0u.new <- lapply(l0u, function(x){
  ll <- length(x); l0.new[1:ll]
})
l0i.new <- c()
l0i.new[which(unlist(Delta) == 0)] <- 0 
l0i.new[which(unlist(Delta) == 1)] <- l0.new[match(survdata[which(unlist(Delta)==1), 'survtime'], sv$ft)]
l0i.new <- as.list(l0i.new)

