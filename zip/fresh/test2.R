rm(list=ls())
library(glmmTMB)
library(survival)
source('./_Functions.R')
source('../inits/inits.R')
source('../survFnsInt.R')
library(Rcpp)
library(RcppArmadillo)
sourceCpp('../zip.cpp')

beta <- c(1.5, 0.05, 0.33, 0.50)
alpha <- c(-0.5, 0.25)
D <- diag(c(.5^2, .15^2))
gamma <- .5
n <- 250
ntms <- 15
data <- simData_zip_joint(n, ntms, beta, alpha, D,
                         theta = c(-4.5, 0.2), gamma = gamma)

# inits
inits.long <- Longit.inits(data$data)
b <- Ranefs(inits.long)
inits.surv <- TimeVarCox(data$data, b)
beta <- inits.long$beta.init
alpha <- inits.long$alpha.init
b <- lapply(1:n, function(i) b[i,])
indzi <- zi.ind <- 2
D <- inits.long$D.init
gamma <- inits.surv$inits[3]
eta <- inits.surv$inits[1:2]

X <- Y <- Z <- Xz <- Zz <- K <- list()
for(i in 1:n){
  i.dat <- data$data[data$dat$id == i, ]
  X[[i]] <- model.matrix(~time + cont + bin, i.dat)
  Xz[[i]] <- model.matrix(~time, i.dat)
  Z[[i]] <- Zz[[i]] <- model.matrix(~1, i.dat)
  Y[[i]] <- i.dat$Y
  K[[i]] <- as.matrix(unname(unique(i.dat[,c('cont', 'bin')])))
  row.names(K[[i]]) <- NULL
}          

ph <- coxph(Surv(survtime, status) ~ cont + bin, distinct(data$data, id, survtime, status, cont, bin))
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

# Checking densities are similar +
joint_density(b[[1]], Y[[1]], X[[1]], Z[[1]], Xz[[1]], Zz[[1]], beta, alpha, D, 2, 
              gamma, K[[1]], eta, l0u[[1]], KK[[1]], Fu[[1]], l0i[1], Delta[1])
b_logdensity(b[[1]], Y[[1]], X[[1]], Z[[1]], Xz[[1]], Zz[[1]], beta, alpha, D, 2)
# + first maximisations for b are also similar.
ucminf::ucminf(b[[1]], b_logdensity, b_score, Y[[1]], X[[1]], Z[[1]], Xz[[1]], Zz[[1]], 
               beta, alpha, D, indzi)$par
ucminf::ucminf(b[[1]], joint_density, NULL, Y[[1]], X[[1]], Z[[1]], Xz[[1]], Zz[[1]], 
               beta, alpha, D, indzi, gamma, K[[1]], eta, l0u[[1]], KK[[1]], Fu[[1]], l0i[1], Delta[1], hessian = 2)

ucminf::ucminf(b[[1]], joint_density, joint_density_ddb, Y[[1]], X[[1]], Z[[1]], Xz[[1]], Zz[[1]], 
               beta, alpha, D, indzi, gamma, K[[1]], eta, l0u[[1]], KK[[1]], Fu[[1]], l0i[1], Delta[1])$par

solve(joint_density_sdb(b[[1]],Y[[1]], X[[1]], Z[[1]], Xz[[1]], Zz[[1]], 
                  beta, alpha, D, indzi, gamma, K[[1]], eta, l0u[[1]], KK[[1]], Fu[[1]], l0i[1], Delta[1], eps=1e-4))

# df/dgamma
Delta[1] * sum(b[[1]]) - crossprod(l0u[[1]] * exp(KK[[1]] %*% eta + gamma * Fu[[1]] %*% b[[1]]), Fu[[1]] %*% b[[1]])
# df/deta
Delta[1] * t(K[[1]]) - crossprod(KK[[1]], l0u[[1]] * exp(KK[[1]] %*% eta + Fu[[1]] %*% (gamma * b[[1]])))

gh <- statmod::gauss.quad.prob(3, 'normal')
w <- gh$w; v <- gh$n

EMupdate <- function(b, Y, X, Z, Xz, Zz, 
                     beta, D, alpha, indzi,
                     gamma, K, eta, l0u, KK, Fu, l0i, Delta, w, v){
  
  # E-step
  b.hat <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta){
    ucminf::ucminf(b, joint_density, joint_density_ddb, 
                   Y, X, Z, Xz, Zz, beta, alpha, D, indzi, gamma, K, eta, l0u, KK, Fu, l0i, Delta)$par
  }, b = b, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
     KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), SIMPLIFY = F)
  
  Sigmai <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta){
    solve(joint_density_sdb(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, gamma, K, eta, l0u, KK, Fu, l0i, Delta, eps = 1e-4))
  }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
     KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), SIMPLIFY = F)
  S <- lapply(Sigmai, function(y) lapply(split(1:2, c(1,2)), function(x) as.matrix(y[x,x])))

  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, b = b.hat, S = Sigmai, SIMPLIFY = F)
  
  ba <- mapply(function(b, Y, X, Z, Xz, Zz, S){
    beta_alpha_update(beta, alpha, b, Y, X, Z, Xz, Zz, D, S, indzi, w, v, eps = 1e-4)
  }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, S = S, SIMPLIFY = F)
  
  tau <- mapply(function(Fu, S){
    sqrt(diag(tcrossprod(Fu %*% S, Fu)))
  }, Fu = Fu, S = Sigmai, SIMPLIFY = F)
  
  gammaeta <- mapply(function(b, Y, X, Z, Xz, Zz, K, l0u, KK, Fu, l0i, Delta, tau){
    gamma_eta_update(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, gamma, K, eta, 
                     l0u, KK, Fu, l0i, Delta, tau, w, v)
  },  b = b.hat, Y = Y, X = X, Z = Z, Xz = Xz, Zz = Zz, K = K, l0u = l0u, 
      KK = KK, Fu = Fu, l0i = as.list(l0i), Delta = as.list(Delta), tau = tau, SIMPLIFY = F)
  
  # M-step
  # ZIP and random effects part.
  Sbeta <- Reduce('+', lapply(ba, '[[', 1))
  Hbeta <- Reduce('+', lapply(ba, '[[', 3))
  Salpha <- Reduce('+', lapply(ba, '[[', 2))
  Halpha <- Reduce('+', lapply(ba, '[[', 4))
  
  D.new <- Reduce('+', Drhs)/n
  beta.new <- beta - solve(Hbeta, Sbeta)
  alpha.new <- alpha - solve(Halpha, Salpha)
  
  # Survival part
  Sgamma <- sum(do.call(c, lapply(gammaeta, function(x) x$Sgamma)))
  Seta <- rowSums(do.call(cbind, lapply(gammaeta, function(x) x$Seta)))
  Hgamma <- sum(do.call(c, lapply(gammaeta, function(x) x$Hgamma)))
  Heta <- Reduce('+', lapply(gammaeta, function(x) x$Heta))
  Hgammaeta <- rowSums(do.call(cbind, lapply(gammaeta, function(x) x$Hgammaeta)))
  
  # Form score vector and information matrix
  Sgammaeta <- c(Sgamma, Seta)
  Imat <- as.matrix(Matrix::bdiag(Hgamma, Heta))
  Imat[2:3, 1] <- Imat[1, 2:3] <- Hgammaeta
  gammaeta.new <- c(gamma, eta) - solve(Imat, Sgammaeta)
  
  # Update for baseline hazard
  
  lambda <- lambdaUpdate(surv.times, ft, gamma, eta, K, Sigmai, b.hat, n, w, v)
  l0.new <- sv$nev/rowSums(lambda)
  l0.new <- nev/rowSums(lambda)
  l0u.new <- lapply(l0u, function(x){
    ll <- length(x); l0.new[1:ll]
  })
  l0i.new <- c()
  l0i.new[which(Delta == 0)] <- 0 
  l0i.new[which(Delta == 1)] <- l0.new[match(data$surv.data[which(Delta==1), 'survtime'], ft)]
  
  return(list(
    D.new = D.new, beta.new = beta.new, alpha.new = alpha.new, b.hat = b.hat,
    gamma.new = gammaeta.new[1], eta.new = gammaeta.new[2:3],
    l0.new = l0.new, l0u.new = l0u.new, l0i.new = l0i.new
  ))
}

# update <- EMupdate(b, Y, X, Z, Xz, Zz, beta, D, alpha, 2, gamma, K, eta, l0u, KK, Fu, l0i, Delta, w, v)

vech <- function(x) x[lower.tri(x, diag = T)]
params <- c(vech(D), beta, alpha, gamma, eta)
diff <- 100; tol <- 1e-2; iter <- 0
while(diff > tol){
  update <- EMupdate(b, Y, X, Z, Xz, Zz, beta, D, alpha, 2, gamma, K, eta, l0u, KK, Fu, l0i, Delta, w, v)
  params.new <- c(vech(update$D.new), update$beta.new, update$alpha.new, update$gamma.new, update$eta.new)
  diff <- max(
    abs(params.new-params)/(abs(params)+1e-3)
  )
  message('\nIteration ', iter + 1)
  print(sapply(params.new, round, 4))
  message('Maximum relative difference: ', round(diff, 5))
  b <- update$b.hat
  D <- update$D.new
  beta <- update$beta.new
  alpha <- update$alpha.new
  gamma <- update$gamma.new
  eta <- update$eta.new
  l0 <- update$l0.new; l0u <- update$l0u.new; l0i <- update$l0i.new
  params <- c(vech(D), beta, alpha, gamma, eta)
  iter <- iter + 1
}
params
params.start
