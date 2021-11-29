rm(list=ls())
setwd('~/Documents/GLMM/zip/zhu')
library(glmmTMB)
source('./_Functions.R')
source('./inits.R')
source('../survFnsInt.R')
sourceCpp('./zip.cpp')

data <- simData_zip_joint(n = 250, theta = c(-3, 0.25)); n <- 250
ph <- coxph(Surv(survtime, status) ~ bin, data = data$surv.data)
inits.long <- Longit.inits(data$data)
beta <- inits.long$beta.init
alpha <- inits.long$alpha.init
D <- inits.long$D.init
b <- Ranefs(inits.long)
inits.surv <- TimeVarCox(data$data, b)
eta <- inits.surv$inits[1]
gamma <- inits.surv$inits[2:3]

sv <- surv.mod(ph, data$data, inits.surv$l0.init)

b <- lapply(1:n, function(i) b[i,])
indzi <- zi.ind <- 2

X <- Y <- Z <- Xz <- Zz <- K <- list()
for(i in 1:n){
  i.dat <- data$data[data$data$id == i, ]
  X[[i]] <- model.matrix(~time + bin, i.dat)
  Xz[[i]] <- model.matrix(~time + bin, i.dat)
  Z[[i]] <- Zz[[i]] <- model.matrix(~1, i.dat)
  Y[[i]] <- i.dat$Y
  K[[i]] <- as.matrix(unname(unique(i.dat[,'bin'])))
  row.names(K[[i]]) <- NULL
}      

sv <- surv.mod(ph, data$data, inits.surv$l0.init)
surv.times <- sv$surv.times; ft <- sv$ft; nev <- sv$nev
Fu <- sv$Fu
Delta <- sv$Di
l0i <- sv$l0i
l0u <- sv$l0u
KK <- sapply(1:n, function(i){
  x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
  if('numeric' %in% class(x)) x <- t(as.matrix(x))
  x
})

aa <- statmod::gauss.quad.prob(9, 'normal')
w <- aa$w; v <- aa$n

b.hat <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta){
  ucminf::ucminf(b, joint_density, joint_density_ddb,
                 Y, X, Z, Xz, Zz, beta, alpha, D, indzi, gamma, K, eta, l0u, KK, Fu, l0i, Delta)$par
}, b = b, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), SIMPLIFY = F)

Sigmai <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta){
  solve(joint_density_sdb(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, gamma, K, eta, l0u, KK, Fu, l0i, Delta, eps = 1e-4))
}, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), SIMPLIFY = F)
                
Drhs <- mapply(function(b, S){
  S + tcrossprod(b)
}, b = b.hat, S = Sigmai, SIMPLIFY = F)

# Necessary steps for update to (beta, alpha)
# Sba <- mapply(function(b, Y, X, Z, Xz, Zz, S){
#   S2betaalpha(c(beta, alpha), b, Y, X, Z, Xz, Zz, length(beta), length(alpha),
#              indzi, S, w, v)
# }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, S = Sigmai, SIMPLIFY = F)
# 
# Hba <- mapply(function(b, Y, X, Z, Xz, Zz, S){
#   Hbetaalpha(c(beta, alpha), b, Y, X, Z, Xz, Zz, length(beta), length(alpha),
#               indzi, S, w, v, eps = 1e-4)
# }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, S = Sigmai, SIMPLIFY = F)

ba <- mapply(function(b, Y, X, Z, Xz, Zz){
  out <- list()
  out[[1]] <- pracma::grad(b_logdensity2, c(beta, alpha), b = b, Y=Y,X=X,Z=Z,Xz=Xz,Zz=Zz,
                           beta_length = length(beta), alpha_length = length(alpha), D=D, indzi=2)
  out[[2]] <- pracma::hessian(b_logdensity2, c(beta, alpha), b = b, Y=Y,X=X,Z=Z,Xz=Xz,Zz=Zz,
                              beta_length = length(beta), alpha_length = length(alpha), D=D, indzi=2)
  out
}, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, SIMPLIFY = F)

# Necessary steps for update to (gamma, eta)
tau <- mapply(function(Fu, S){
  diag(tcrossprod(Fu %*% S, Fu))
}, Fu = Fu, S = Sigmai, SIMPLIFY = F)

Sge <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta, tau){
  Sgammaeta(c(gamma,eta), b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, K, 
            l0u, KK, Fu, l0i, Delta, tau, w, v)
},  b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), tau = tau, SIMPLIFY = F)

Hge <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta, tau){
  Hgammaeta(c(gamma,eta), b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, K, 
            l0u, KK, Fu, l0i, Delta, tau, w, v,1e-4)
},  b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), tau = tau, SIMPLIFY = F)

# Updates
# D
D.new <- Reduce('+', Drhs)/n
# \beta and \alpha
Sba <- colSums(do.call(rbind, lapply(ba, '[[', 1)))
Hba <- Reduce('+', lapply(ba, '[[', 2))
# Sba <- rowSums(do.call(cbind, Sba))
# Hba <- Reduce('+', Hba)
beta.alpha.new <- c(beta,alpha) - solve(Hba, Sba)
# \gamma and \eta
Sge <- rowSums(do.call(cbind, Sge))
Hge <- Reduce('+', Hge)

gamma.eta.new <- c(gamma,eta)-solve(Hge,Sge)

  
  

