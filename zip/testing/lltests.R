tr <- function(x) sum(diag(x))
test <- simData()
fit <- glmmTMB(y ~ time + cont + bin + (1 + time|id), 
               data = test, 
               family = poisson, 
               ziformula = ~ time + (1 + time|id))
ll.true <- logLik(fit)

X <- model.matrix(~time+cont+bin, test)
Y <- test$y
Z <- model.matrix(~time, test)
Xzi <- model.matrix(~time, test)
Zzi <- model.matrix(~time, test)

alpha <- fixef(fit)$zi
beta <- fixef(fit)$cond
b <- as.matrix(ranef(fit)$cond$id)
bzi <- as.matrix(ranef(fit)$zi$id)

Dl <- as.matrix(VarCorr(fit)$cond$id)
Dz <- as.matrix(VarCorr(fit)$zi$id)

WZ <- Xzi %*% alpha + rowSums(Zzi * bzi[test$id, , drop=F])
XZ <- X %*% beta + rowSums(Z * b[test$id, , drop = F])
# Index for non-zero Y
ii <- Y>0
l1 <- sum(log(exp(WZ[!ii]) + exp(-exp(XZ[!ii]))))
l2 <- crossprod(Y[ii], XZ[ii]) - sum(exp(XZ[ii])) - sum(lfactorial(Y[ii]))
l3 <- sum(log(1 + exp(WZ)))

b.ll <- -1/2 * log(2*pi) -1/2 * log(det(Dl)) - 1/2 * tr(b %*% solve(Dl) %*% t(b))
bz.ll <- -1/2 * log(2*pi) - 1/2 * log(det(Dz)) - 1/2 * tr(bzi %*% solve(Dz) %*% t(bzi))

# Adding on theRE ll's?
l1 + l2 - l3 + b.ll + bz.ll

ll.true # A little bit off(?)

l1+l2-l3
