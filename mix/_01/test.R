#' #####
#' Testing functions for joint model poisson test.
#' ####
rm(list=ls())
source('./simData.R')
source('./inits.R')
source('./survFns.R')
library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
sourceCpp('mix.cpp')

set.seed(123)
testsim <- simData_joint()
test <- testsim$data
survdata <- testsim$surv.data
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)

n <- 250
m <- Y <- X <- Z <- K <- list()
for(i in 1:n){
  i.dat <- test[test$id == i, ]
  m[[i]] <- nrow(i.dat)
  Y[[i]] <- cbind(i.dat$Y.1, i.dat$Y.2, i.dat$Y.3)
  X[[i]] <- model.matrix(~time+cont+bin, i.dat)
  Z[[i]] <- model.matrix(~time, i.dat)
  K[[i]] <- cbind(unique(i.dat$cont), unique(i.dat$bin))
}

inits.long <- Longit.inits(test)
b <- Ranefs(inits.long)
beta <- inits.long$beta.init
var.e <- inits.long$var.e.init
V <- lapply(m, function(x) diag(x = var.e, nr = x, nc = x))
theta <- inits.long$theta.init
D <- inits.long$D.init
inits.surv <- TimeVarCox(test, b)
b <- lapply(1:n, function(i) b[i, ])
rm(inits.long) # large object

gamma <- inits.surv$inits[3:5]
eta <- inits.surv$inits[1:2]

#svmod
sv <- surv.mod(ph, test, inits.surv$l0.init)
Fi <- lapply(1:n, function(i) do.call(c, replicate(3, sv$Fi[i, ], simplify = F)))
Fu <- sv$Fu
KK <- sapply(1:n, function(i){
  x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
  if('numeric'%in%class(x)) x <- t(as.matrix(x))
  x
})
l0i <- as.list(sv$l0i); l0u <- sv$l0u
Delta <- as.list(sv$Di)

# Estimates for b and subsequently for Sigmai
b.hat <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, haz){
  ucminf::ucminf(b, joint_density, joint_density_ddb,
                 X = X, Z = Z, beta = beta, V = V, D = D,
                 Y_1 = Y[, 1], Y_2 = Y[, 2], Y_3 = Y[, 3], F, 0.0,
                 Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu,
                 haz = haz, gamma = rep(gamma, each = 2), eta = eta)$par
}, b = b, X = X, Y = Y, Z = Z, V = V, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
KK = KK, Fu = Fu, haz = l0u, SIMPLIFY = F)
# RE-related objects, bmat and bsplit
bmat <- lapply(b.hat, matrix, nc = 2, byr = T)
bsplit <- lapply(b.hat, function(y) lapply(split(seq(6), rep(seq(3), each = 2)), function(x) y[x]))

Sigmai <- mapply(function(b, X, Z, V, Y, Delta, K, Fi, l0i, KK, Fu, l0u){
  solve(joint_density_sdb(b = b, X = X, Z = Z, beta = beta, V = V, D = D,
                          Y_1 = Y[,1], Y_2 = Y[,2], Y_3 = Y[,3], F, 0.0, 
                          Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu, 
                          haz = l0u, gamma = rep(gamma, each = 2), eta = eta, eps = 1e-3))
}, b = b.hat, X = X, Z = Z, V = V, Y = Y, Delta = Delta, K = K, Fi = Fi,
   l0i= l0i, KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)
# Split out into constituent block-diag pieces.
SigmaiSplit <- lapply(Sigmai, function(y) lapply(split(seq(6), rep(seq(3), each = 2)), function(x) y[x,x]))

# check diags
diags <- do.call(rbind, lapply(Sigmai, diag))
checkzerodiag <- apply(diags, 2, function(x) sum(x <= 0))
if(any(checkzerodiag > 0)) warning('Negative diagonal elements, sdb misspecified')

# Update to D
Drhs <- mapply(function(b, S){
  S + tcrossprod(b)
}, S = Sigmai, b = b.hat, SIMPLIFY = F)

# beta
Sb <- mapply(function(X, Y, Z, b, V){
  Sbeta(beta, X, Y[,1], Y[,2], Y[,3], Z, b, V, F, 0)
}, X = X, Y = Y, Z = Z, b = b.hat, V = V, SIMPLIFY = F)

Hb <- mapply(function(X, Y, Z, b, V){
  Hbeta(beta, X, Y[,1], Y[,2], Y[,3], Z, b, V, F, 0, 1e-4)
}, X = X, Y = Y, Z = Z, b = b.hat, V = V, SIMPLIFY = F)

# var.e (Gaussian submodel)
tau <- mapply(function(Z, S){
  unname(sqrt(diag(tcrossprod(Z %*% S[[1]], Z))))
}, Z = Z, S = SigmaiSplit, SIMPLIFY = F)

var.e.update <- mapply(function(b, Y, X, Z, tau){
  out <- numeric(2)
  for(i in 1:9){
    out[1] <- out[1] + (w[i] * crossprod(Y[,1] - X %*% beta[1:4] - Z %*% b[1:2] - tau * v[i]))
  }
  out[2] <- length(Y[,1])
  out
}, b = b.hat, Y = Y, X = X, Z = Z, tau = tau)

var.e.update2 <- mapply(function(Y, X, Z, b, tau){
  vare_update(Y[,1], X, Z, beta[1:4], b[[1]], tau, w, v)
}, Y = Y, X = X, Z = Z, b = bsplit, tau = tau)

# gamma and eta
aa <- statmod::gauss.quad.prob(9, 'normal')
w <- aa$w; v <- aa$n
Sge <- mapply(function(b, S, K, KK, Fu, Fi, l0u, Delta){
  Sgammaeta(c(gamma, eta), b, S, K, KK, Fu, Fi[1:2], l0u, Delta, w, v, 1e-4)
}, b = bmat, S = SigmaiSplit, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)

Hge <- mapply(function(b, S, K, KK, Fu, Fi, l0u, Delta){
  Hgammaeta(c(gamma, eta), b, S, K, KK, Fu, Fi[1:2], l0u, Delta, w, v, 1e-4)
}, b = bmat, S = SigmaiSplit, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)

# updates
# Longitudinal parts
D.new <- Reduce('+', Drhs)/n
beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
# Residual variance for var.e
var.e.update <- rowSums(var.e.update2)
var.e.new <- var.e.update[1]/var.e.update[2]
gamma.eta.new <- c(gamma, eta) - solve(Reduce('+', Hge), rowSums(Sge))
lambda <- lambdaUpdate(sv$surv.times, sv$ft, gamma, eta, K, SigmaiSplit, bsplit, w, v)

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

