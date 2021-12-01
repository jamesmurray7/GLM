setwd('~/Documents/GLMM/bin/')
library(survival)
library(Rcpp)
library(RcppArmadillo)

source('_Functions.R')
source('survFns.R')
source('inits.R')
sourceCpp('bin.cpp')

# Simdata
dd <- simData_joint()
ph <- coxph(Surv(survtime, status) ~ cont + bin, dd$surv.data)
n <- nrow(dd$surv.data)
inits.long <- Longit.inits(dd$data)
beta <- inits.long$beta.init
D <- inits.long$D.init
b <- Ranefs(inits.long)
inits.surv <- TimeVarCox(dd$data, b)
b <- lapply(1:n, function(i) b[i, ])
gamma <- inits.surv$inits[3]
eta <- inits.surv$inits[1:2]

# Get model matrices
K <- Y <- X <- Z <- list()
for(i in 1:n){
  i.dat <- dd$data[dd$data$id == i, ]
  Y[[i]] <- i.dat$Y
  X[[i]] <- model.matrix(~time + cont + bin, i.dat)
  Z[[i]] <- model.matrix(~time, i.dat)
  K[[i]] <- unname(cbind(unique(i.dat$cont), unique(i.dat$bin)))
}

# Survival objects
sv <- surv.mod(ph, dd$data, inits.surv$l0.init)
Fi <- sv$Fi; Fu <- sv$Fu;
Delta <- as.list(sv$Di)
l0i <- as.list(sv$l0i)
l0u <- sv$l0u
Fi <- lapply(1:n, function(i) Fi[i, ])
KK <- sapply(1:n, function(i){
  x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
  if('numeric' %in% class(x)) x <- t(as.matrix(x))
  x
})


# Check ll
ll <- c()
for(i in 1:n){
  ll[i] <- joint_densitytemp(b[[i]], X[[i]], Y[[i]], Z[[i]], beta, D)
}

b.hat <- mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, l0u){ 
  ucminf::ucminf(b, joint_density, joint_density_ddb,
                 X = X, Y = Y, Z = Z, beta = beta, D = D, Delta = Delta, K = K,
                 Fi = Fi, l0i = l0i, KK = KK, Fu = Fu, haz = l0u, gamma = gamma, eta = eta)$par
}, b = b, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)

Sigmai <-  mapply(function(b, X, Y, Z, Delta, K, Fi, l0i, KK, Fu, l0u){ 
  solve(-joint_density_sdb(b, X = X, Y = Y, Z = Z, beta = beta, D = D, Delta = Delta, K = K,
                    Fi = Fi, l0i = l0i, KK = KK, Fu = Fu, haz = l0u, gamma = gamma, eta = eta, eps = 1e-4))
}, b = b.hat, X = X, Y = Y, Z = Z, Delta = Delta, K = K, Fi = Fi, l0i = l0i, KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)


# D
Drhs <- mapply(function(S, b){
  S + tcrossprod(b)
}, S = Sigmai, b = b.hat, SIMPLIFY = F)



# \beta

# qfun <- function(beta, X, Z, Y, b){
#   eta <- X %*% beta + Z %*% b;
#   log_density(eta, Y)
# }

beta.update <- mapply(function(X, Y, Z, b){
  out <- list()
  out[[1]] <- Sbeta(beta, X, Y, Z, b)
  out[[2]] <- Hbeta(beta, X, Y, Z, b, eps = 1e-4)
#  out[[3]] <- pracma::grad(qfun, beta, X = X, Y = Y, Z = Z, b = b)
#  out[[4]] <- pracma::hessian(qfun, beta, X = X, Y = Y, Z = Z, b = b)
  out
}, X = X, Y = Y, Z = Z, b = b.hat, SIMPLIFY = F)

beta.new <- beta - solve(Reduce('+', lapply(beta.update, '[[', 2)),
                         rowSums(do.call(cbind, lapply(beta.update, '[[', 1))))

# \gamma and \eta
a <- statmod::gauss.quad.prob(9, 'normal')
w <- a$w; v <- a$n

Sge <- mapply(function(b, Delta, Fi, K, KK, Fu, l0u, S){
  Sgammaeta(c(gamma, eta), b, Delta, Fi, K, KK, Fu, l0u, S, w, v)
}, b = b.hat, Delta = Delta, Fi = Fi, K = K, KK = KK, Fu = Fu, l0u = l0u, S = Sigmai)

Hge <- mapply(function(b, Delta, Fi, K, KK, Fu, l0u, S){
  Hgammaeta(c(gamma, eta), b, Delta, Fi, K, KK, Fu, l0u, S, w, v, 1e-4)
}, b = b.hat, Delta = Delta, Fi = Fi, K = K, KK = KK, Fu = Fu, l0u = l0u, S = Sigmai, SIMPLIFY = F)

gammaeta.new <- c(gamma, eta) - solve(Reduce('+', Hge), rowSums(Sge))

lambda <- lambdaUpdate(sv$surv.times, sv$ft, gamma, eta, K, Sigmai, b.hat, w, v)
l0 <- sv$nev/rowSums(lambda)
