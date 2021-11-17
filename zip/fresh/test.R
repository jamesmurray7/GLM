rm(list=ls())
library(glmmTMB)
source('~/Documents/GLMM/zip/fresh/_Functions.R')

beta <- c(-.5, .1, 0.5, 1)
alpha <- c(-5, 0.25)
n <- 250; ntms <- 8
D <- list(D1 = diag(c(.5^2, .25^2)), D2 = diag(c(.2^2, .05^2)))
D <- diag(c(.5^2, .15^2))
test <- simData_zip(n, ntms, beta, alpha, D)
testfit <- glmmTMB::glmmTMB(y ~ time  +cont  +bin + (1|id),
                            ziformula = ~ time + (1|id), data = test, family = poisson)

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
