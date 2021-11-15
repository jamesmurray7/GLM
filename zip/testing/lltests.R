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

getL <- function(data, i){
  i.dat <- data[data$id == i, ]
  b <- b[i,]; bzi <- bzi[i,]
  X <- model.matrix(~time + cont + bin, i.dat)
  Z <- Xzi <- Zzi <- model.matrix(~time, i.dat)
  Y <- i.dat$y
  ii <- Y == 0
  if(sum(ii) > 0){
    l1 <- sum(exp(Xzi[ii,] %*% alpha + Zzi[ii,] %*% bzi) + exp(-exp(X[ii,] %*% beta + Z[ii,] %*% b)))
  }else{
    l1 <- 0
  }
  ii <- Y > 0 
  if(sum(ii) > 0){
    l2 <- crossprod(Y[ii], X[ii,] %*% beta + Z[ii,] %*% b) - sum(exp(X[ii,] %*% beta + Z[ii,] %*% b)) - sum(lfactorial(Y[ii]))
  }else{
    l2 <- 0
  }
  l3 <- sum(log(1 + exp(Xzi %*% alpha + Zzi %*% bzi)))
  list(l1, l2, l3, l1 + l2 - l3)
}

aa <- c()
for(i in 1:250) aa[i] <- getL(test, i)[[4]]

getL1Derivs <- function(data, i){
  i.dat <- data[data$id == i, ]
  b <- b[i,]; bzi <- bzi[i,]
  X <- model.matrix(~time + cont + bin, i.dat)
  Z <- Xzi <- Zzi <- model.matrix(~time, i.dat)
  Y <- i.dat$y
  # L1
  ii <- Y == 0
  diag <- function(x) base::diag(x, nc = sum(ii))
  if(sum(ii) > 0){
    XZ <- X[ii,] %*% beta + Z[ii,] %*% b
    ZI <- Xzi[ii,] %*% alpha + Zzi[ii,] %*% bzi
    # Scores
    Sbeta1 <- -crossprod(X[ii,,drop=F],
                         exp(-exp(XZ)) * exp(XZ)/(exp(-exp(XZ)) + exp(ZI)))
    Salpha1 <- crossprod(Xzi[ii,,drop=F],
                         exp(ZI) / (exp(ZI) + exp(-exp(XZ))))
    Sb1 <- -crossprod(Z[ii,,drop=F],
                      exp(-exp(XZ)) * exp(XZ)/(exp(-exp(XZ)) + exp(ZI)))
    Sbzi1 <- crossprod(Zzi[ii,,drop=F],
                       exp(ZI) / (exp(ZI) + exp(-exp(XZ))))
    # Second Derivatives
    # I(beta) and I(b) #
    t0 <- exp(XZ); t1 <- exp(-t0); t2 <- t1 + exp(ZI); t3 <- t1 * t0
    Ibeta1 <- -1 * (crossprod(diag(x = c(t3 / t2)) %*% X[ii,,drop=F], X[ii,,drop=F]) - 
              crossprod(diag(x = c(t0 * t1 * t0 / t2)) %*% X[ii,,drop=F], X[ii,,drop=F]) + 
              crossprod(diag(x = c(t3 * t1 * t0 / (t2 * t2))) %*% X[ii,,drop=F], X[ii,,drop=F]) )
    Ib1 <- -1 * (crossprod(diag(x = c(t3 / t2)) %*% Z[ii,,drop=F], Z[ii,,drop=F]) - 
                 crossprod(diag(x = c(t0 * t1 * t0 / t2)) %*% Z[ii,,drop=F], Z[ii,,drop=F]) + 
                 crossprod(diag(x = c(t3 * t1 * t0 / (t2 * t2))) %*% Z[ii,,drop=F], Z[ii,,drop=F]) )
    t0 <- exp(ZI); t1 <- t0 + exp(-exp(XZ))
    # I(alpha) and I(bzi) #
    Ialpha1 <- crossprod(diag(x = c(t0 / t1)) %*% Xzi[ii,,drop=F], Xzi[ii,,drop=F]) - 
               crossprod(diag(x = c(t0 * t0 / (t1 * t1))) %*% Xzi[ii,,drop=F], Xzi[ii,,drop=F])
    Ibzi1 <- crossprod(diag(x = c(t0 / t1)) %*% Zzi[ii,,drop=F], Zzi[ii,,drop=F]) - 
              crossprod(diag(x = c(t0 * t0 / (t1 * t1))) %*% Zzi[ii,,drop=F], Zzi[ii,,drop=F])
    # Cross-terms #
    t0 <- exp(XZ); t1 <- exp(ZI); t2 <- exp(-t0); t3 <- t1 + t2
    Ibbzi1 <- crossprod(diag(x = c(t1 * t2 * t0 / (t3 * t3))) %*% Zzi[ii,,drop=F], Z[ii,,drop=F])
  }else{
    Sbeta1 <- rep(0, length(beta))
    Salpha1 <- rep(0, length(alpha))
    Sb1 <- rep(0, length(b)); Sbzi1 <- rep(0, length(bzi))
    Ibeta1 <- Ialpha1 <- Ib1 <- Ibzi1 <- Ibbzi1 <- 0
  }
  return(list(
    Sbeta1, Salpha1, Sb1, Sbzi1, Ibeta1, Ialpha1, Ib1, Ibzi1, Ibbzi1
  ))
}

getL1Derivs(test, 1)
