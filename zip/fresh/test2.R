rm(list=ls())
library(glmmTMB)
library(survival)
source('./_Functions.R')
source('../inits/inits.R')
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
               beta, alpha, D, indzi, gamma, K[[1]], eta, l0u[[1]], KK[[1]], Fu[[1]], l0i[1], Delta[1])$par

ucminf::ucminf(b[[1]], joint_density, joint_density_ddb, Y[[1]], X[[1]], Z[[1]], Xz[[1]], Zz[[1]], 
               beta, alpha, D, indzi, gamma, K[[1]], eta, l0u[[1]], KK[[1]], Fu[[1]], l0i[1], Delta[1])$par

# df/dgamma
Delta[1] * sum(b[[1]]) - crossprod(l0u[[1]] * exp(KK[[1]] %*% eta + gamma * Fu[[1]] %*% b[[1]]), Fu[[1]] %*% b[[1]])

# df/deta
Delta[1] * t(K[[1]]) - crossprod(KK[[1]], l0u[[1]] * exp(KK[[1]] %*% eta + Fu[[1]] %*% (gamma * b[[1]])))

# TO DO:: 
# score for b from JOINT density
# also parameter updates for gamma and eta.

a <- sqrt(diag(gamma^2 * tcrossprod(Fu[[1]] %*% Sigmai, Fu[[1]])))
b <- gamma * sqrt(diag(tcrossprod(Fu[[1]] %*% Sigmai, Fu[[1]])))
EMupdateCpp <- function(b, Y, X, Z, Xz, Zz, 
                     beta, D, alpha, gh, indzi = 2){
  
  b.hat <- mapply(function(b, Y, X, Z, Xz, Zz){
    ucminf::ucminf(b, b_logdensity, b_score, Y, X, Z, Xz, Zz, 
                   beta, alpha, D, indzi)$par
  }, b = b, Y = Y, X = X, Z = Z, Xz = Xzi, Zz = Zzi, SIMPLIFY = F)
  
  Sigmai <- mapply(function(b, Y, X, Z, Xz, Zz){
    solve(fd_b(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, 1e-4))
  }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xzi, Zz = Zzi, SIMPLIFY = F)
  
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, b = b.hat, S = Sigmai, SIMPLIFY = F)
  
  ba <- mapply(function(b, Y, X, Z, Xz, Zz, Sigmai){
    beta_alpha_new(b, Y, X, Z, Xz, Zz, beta, D, alpha, Sigmai, gh)
  }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xzi, Zz = Zzi, Sigmai = Sigmai, SIMPLIFY = F)
  
  # Updates
  Sbeta <- Reduce('+', lapply(ba, '[[', 1))
  Hbeta <- Reduce('+', lapply(ba, '[[', 3))
  Salpha <- Reduce('+', lapply(ba, '[[', 2))
  Halpha <- Reduce('+', lapply(ba, '[[', 4))
  
  D.new <- Reduce('+', Drhs)/n
  beta.new <- beta - solve(Hbeta, Sbeta)
  alpha.new <- alpha - solve(Halpha, Salpha)
  
  return(list(
    D.new = D.new, beta.new = beta.new, alpha.new = alpha.new, b.hat = b.hat
  ))
}
cd <- GLMMadaptive:::cd_vec
diff <- 100; iter <- 0; tol <- 5e-3
beta <- fixef(testfit)$cond
alpha <- fixef(testfit)$zi
update <- EMupdateCpp(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, gh = 9)
params.new <- c(vech(update$D.new), update$beta.new, update$alpha)
diff <- max(
  abs(params.new-params)/(abs(params)+1e-3)
)
message('\nIteration ', iter + 1)
message('Maximum relative difference: ', round(diff, 5))
b <- update$b.hat
D <- update$D.new
beta <- update$beta.new
alpha <- update$alpha.new
params <- c(vech(D), beta, alpha)
iter <- iter + 1

b <- as.matrix(cbind(ranef(testfit)$cond$id, ranef(testfit)$zi$id))
b <- lapply(1:n, function(i) b[i,])
indzi <- zi.ind <- 2
D <- diag(sapply(1:2, function(x) VarCorr(testfit)[[x]]$id))



vech <- function(x) x[lower.tri(x, diag = T)]
params <- c(vech(D), beta, alpha)
diff <- 100
while(diff > tol){
  update <- EMupdateCpp(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, gh = 9)
  params.new <- c(vech(update$D.new), update$beta.new, update$alpha)
  diff <- max(
    abs(params.new-params)/(abs(params)+1e-3)
  )
  message('\nIteration ', iter + 1)
  message('Maximum relative difference: ', round(diff, 5))
  b <- update$b.hat
  D <- update$D.new
  beta <- update$beta.new
  alpha <- update$alpha.new
  params <- c(vech(D), beta, alpha)
  iter <- iter + 1
}
params
params.start
