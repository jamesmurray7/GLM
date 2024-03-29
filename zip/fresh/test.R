rm(list=ls())
library(glmmTMB)
source('./_Functions.R')

beta <- c(4, 0.5, 1, 3)
alpha <- c(-2, 0.25)
n <- 250; ntms <- 15
D <- diag(c(.5^2, .15^2))
test <- simData_zip(n, ntms, beta, alpha, D)
testfit <- glmmTMB::glmmTMB(Y ~ time  +cont  +bin + (1|id),
                            ziformula = ~ time + (1|id), data = test, family = poisson)
# pscl::zeroinfl(y~time+cont+bin | time,data=test) # Could we do something with this and speed up a little?
beta <- fixef(testfit)$cond
alpha <- fixef(testfit)$zi
b <- as.matrix(cbind(ranef(testfit)$cond$id, ranef(testfit)$zi$id))
b <- lapply(1:n, function(i) b[i,])
indzi <- zi.ind <- 2
D <- diag(sapply(1:2, function(x) VarCorr(testfit)[[x]]$id))

X <- Y <- Z <- Xzi <- Zzi <- list()
for(i in 1:n){
  i.dat <- test[test$id == i, ]
  X[[i]] <- model.matrix(~time + cont + bin, i.dat)
  Xzi[[i]] <- model.matrix(~time, i.dat)
  Z[[i]] <- Zzi[[i]] <- model.matrix(~1, i.dat)
  Y[[i]] <- i.dat$Y
}          

zi <- zip()

b.hat <- mapply(function(b, Y, X, Z, Xzi, Zzi){
  ucminf::ucminf(b, .b()$logfb, .b()$score.logfb, Y, X, Z, Xzi, Zzi, 
                 beta, D, alpha, zi$ldens, indzi, zi$Seta, zi$Setazi)$par
}, b = b, Y = Y, X = X, Z = Z, Xzi = Xzi, Zzi = Zzi, SIMPLIFY = F)

Sigmai <- mapply(function(b, Y, X, Z, Xzi, Zzi){
  solve(.b()$H(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, zi$ldens, indzi, zi$Seta, zi$Setazi))
}, b = b.hat, Y = Y, X = X, Z = Z, Xzi = Xzi, Zzi = Zzi, SIMPLIFY = F)
                
Drhs <- mapply(function(b, S){
  S + tcrossprod(b)
}, b = b.hat, S = Sigmai, SIMPLIFY = F)

ba <- mapply(function(b, Y, X, Z, Xzi, Zzi, Sigmai){
  beta_alpha(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, zi$ldens, indzi, zi$Seta, zi$Setazi, 9, Sigmai)
}, b = b.hat, Y = Y, X = X, Z = Z, Xzi = Xzi, Zzi = Zzi, Sigmai = Sigmai, SIMPLIFY = F)

# Updates

Sbeta <- Reduce('+', lapply(ba, '[[', 1))
Hbeta <- Reduce('+', lapply(ba, '[[', 3))
Salpha <- Reduce('+', lapply(ba, '[[', 2))
Halpha <- Reduce('+', lapply(ba, '[[', 4))

D.new <- Reduce('+', Drhs)/n
beta.new <- beta - solve(Hbeta, Sbeta)
alpha.new <- alpha - solve(Halpha, Salpha)

# Let's try and functionise this

EMstep <- function(b, Y, X, Z, Xzi, Zzi, 
                   beta, D, alpha, gh){
  b.hat <- mapply(function(b, Y, X, Z, Xzi, Zzi){
    ucminf::ucminf(b, .b()$logfb, .b()$score.logfb, Y, X, Z, Xzi, Zzi, 
                   beta, D, alpha, zi$ldens, indzi, zi$Seta, zi$Setazi)$par
  }, b = b, Y = Y, X = X, Z = Z, Xzi = Xzi, Zzi = Zzi, SIMPLIFY = F)
  
  Sigmai <- mapply(function(b, Y, X, Z, Xzi, Zzi){
    solve(.b()$H(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, zi$ldens, indzi, zi$Seta, zi$Setazi))
  }, b = b.hat, Y = Y, X = X, Z = Z, Xzi = Xzi, Zzi = Zzi, SIMPLIFY = F)
  
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, b = b.hat, S = Sigmai, SIMPLIFY = F)
  
  ba <- mapply(function(b, Y, X, Z, Xzi, Zzi, Sigmai){
    beta_alpha(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, zi$ldens, indzi, zi$Seta, zi$Setazi, gh, Sigmai)
  }, b = b.hat, Y = Y, X = X, Z = Z, Xzi = Xzi, Zzi = Zzi, Sigmai = Sigmai, SIMPLIFY = F)
  
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

update <- EMstep(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, 9)
# unpack
params <- c(vech(D), beta, alpha)
params.new <- c(vech(update$D.new), update$beta.new, update$alpha)
diff <- max(
  abs(params.new-params)/(abs(params)+1e-3)
)
message('Maximum relative difference: ', round(diff, 5))
  
# let's try full-run while looping ----------------------------------------
diff <- 100; iter <- 0; tol <- 5e-3
beta <- fixef(testfit)$cond
alpha <- fixef(testfit)$zi
b <- as.matrix(cbind(ranef(testfit)$cond$id, ranef(testfit)$zi$id))
b <- lapply(1:n, function(i) b[i,])
indzi <- zi.ind <- 2
D <- diag(sapply(1:2, function(x) VarCorr(testfit)[[x]]$id))

X <- Y <- Z <- Xzi <- Zzi <- list()
for(i in 1:n){
  i.dat <- test[test$id == i, ]
  X[[i]] <- model.matrix(~time + cont + bin, i.dat)
  Xzi[[i]] <- model.matrix(~time, i.dat)
  Z[[i]] <- Zzi[[i]] <- model.matrix(~1, i.dat)
  Y[[i]] <- i.dat$y
}          

zi <- zip()

vech <- function(x) x[lower.tri(x, diag = T)]
params <- c(vech(D), beta, alpha)

while(diff > tol){
  update <- EMstep(b, Y, X, Z, Xzi, Zzi, beta, D, alpha, 3)
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

# C++ versions ------------------------------------------------------------
library(Rcpp)
library(RcppArmadillo)
sourceCpp('./wip.cpp')

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
