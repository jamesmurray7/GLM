library(Rcpp);library(RcppArmadillo)

setwd('~/Documents/GLMM/negbinom/')
source('simData.R')
sourceCpp('./ll.cpp')
source('./inits/inits.R')
source('./inits/MVLME.R')
source('../DataFunctions/longFns.R')
vech <- function(x) x[lower.tri(x, diag = T)]
args(simData)
beta <- rbind(c(0.75, -0.05, -0.2, 0.2),
              c(1, 0.05,  0.5, -0.80))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5

gamma <- c(-0.5, 0.8)
eta <- c(-0.3, 0.5)

nK <- 2; q <- 4

data <- simData(250, 15, beta, D, gamma, eta, theta = c(-5, 0.2))   # About 50%

X <- getXi(data, 2); Y <- getYi(data, 2); Z <- getZi(data, 2)
Xk <- splitXks(data, 2); Zk <- splitZks(data, 2); mi <- getmi(data, 2); Yk <- splitYks(Y, mi, 2)

inits.long <- Longit.inits(2, data)
mvlme.step <- mvlme(data, Y, X, Z, inits.long, 2, 4, verbose = T)

beta <- mvlme.step$beta; beta.inds <- list(1:4, 5:8)
theta <- mvlme.step$theta
b <- mvlme.step$b; b.inds <- list(1:2,3:4)
D <- mvlme.step$D

inits.long$beta.init

muk <- mapply(function(X, Z, b){
  out <- list()
  for(k in 1:nK){
    out[[k]] <- exp(X[[k]] %*% beta[beta.inds[[k]]] + Z[[k]] %*% b[b.inds[[k]]])
  }
  out
}, X=Xk,Z=Zk,b=b,SIMPLIFY = F)

Stheta <- mapply(function(Yk, muk){
  out <- list()
  for(k in 1:nK){
    out[[k]] <- sum(digamma(theta[k] + Yk[, k]) - (theta[k] + Yk[,k])/(muk[[k]]+theta[k])+log(muk[[k]]+theta[k])) - digamma(theta[k]) + 1 + log(theta[k])
  }
  out
}, Yk = Yk, muk = muk, SIMPLIFY = F)

Itheta <- mapply(function(Yk, muk){
  out <- list()
  for(k in 1:nK){
    out[[k]] <- sum(trigamma(theta[k] + Yk[, k]) - (Yk[, k] - muk[[k]])/((muk[[k]] + theta[k])^2) + (muk[[k]] + theta[k])^(-1)) - trigamma(theta[k]) + 1/theta[k]
  }
  out
}, Yk = Yk, muk = muk, SIMPLIFY = F)

SthetaSum <- sapply(1:nK, function(i){
  sum(do.call(c, lapply(Stheta, '[[', i)))
})

IthetaMatrix <- diag(
  sapply(1:nK, function(i){
    sum(do.call(c, lapply(Itheta, '[[', i)))
  })
)

theta + solve(IthetaMatrix, SthetaSum)

########################################
########################################
########################################
########################################

# Theta as a vector
theta.list <- lapply(1:n, function(i){
  out <- list()
  for(k in 1:nK) out[[k]] <- rep(theta[k], length(Yk[[i]][, k]))
  out
})

# Take id=1, k=1
t1 <- theta.list[[1]][[1]]
y1 <- Yk[[1]][,1]
m1 <- muk[[1]][[1]]

(t1+y1) %*% log(m1+t1)
crossprod((t1+y1) ,log(m1+t1))

theta_ll <- function(theta, Y, mu){
  sum(lgamma(theta + Y)) - sum(lgamma(theta)) + theta %*% log(theta) - crossprod(theta + Y, log(mu + theta))
}

theta_ll(theta.list[[1]][[1]], Yk[[1]][,1], muk[[1]][[1]])
numDeriv::grad(theta_ll, theta.list[[1]][[1]], Y = Yk[[1]][, 1], mu = muk[[1]][[1]])

Stheta <- function(theta, Y, mu){ 
  digamma(theta+Y) - digamma(theta) + 1 + log(theta) - (log(mu + theta) + (theta + Y)/(mu + theta))
}
# OK this matches with numDeriv above
Stheta(theta.list[[1]][[1]], Yk[[1]][, 1], muk[[1]][[1]])

H <- numDeriv::hessian(theta_ll, theta.list[[1]][[1]], Y = Yk[[1]][, 1], mu = muk[[1]][[1]])
max(H[lower.tri(H)]) # off-diagonal terms essentially zero, good for us. 

Itheta <- function(theta, Y, mu){
  ones <- rep(1, length(theta))
  term1 <- diag(x = c(trigamma(theta + Y))) - diag(x = c(trigamma(theta))) + diag(x = c(ones / theta))
  term2 <- diag(x = c(ones / (mu + theta)))
  term3 <- diag(x = c((theta + Y) / ((mu + theta) * (mu + theta))))
  term1 - (2 * term2 - term3)
}

myH <- Itheta(theta.list[[1]][[1]], Yk[[1]][, 1], muk[[1]][[1]])
all.equal(diag(myH), diag(H))

# Construct an update using this approach
# What's the longest profile? Need to make a matrix of this length for each K
LL <- sapply(1:nK, function(k) max(do.call(c, lapply(lapply(Xk, '[[', k), nrow))))

Sthetas <- Ithetas <- list()
for(i in 1:n){
  thisS <- thisI <- list()
  for(k in 1:nK){
    tempSk <- Stheta(theta.list[[i]][[k]], Yk[[i]][, k], muk[[i]][[k]])
    tempSk <- c(tempSk, rep(0, LL[k] - length(tempSk)))                         # Concatenate with a load of zeroes so we can rbind later
    thisS[[k]] <- tempSk
    tempIk <- -1 * Itheta(theta.list[[i]][[k]], Yk[[i]][, k], muk[[i]][[k]])
    tempIk <- diag(tempIk)
    tempIk <- diag(x = c(tempIk, rep(0, LL[k] - length(tempIk))))               # Concatenate with a load of zeroes so we can rbind later
    thisI[[k]] <- tempIk
  }
  Sthetas[[i]] <- thisS
  Ithetas[[i]] <- thisI
}

St <- lapply(1:nK, function(k) colSums(do.call(rbind, lapply(Sthetas, '[[', k))))
It <- lapply(1:nK, function(k) Reduce('+', lapply(Ithetas, '[[', k)))

rep(theta[1], LL[1]) + solve(It[[1]], St[[1]]) # looks more plausible than previous attempts
rep(theta[2], LL[2]) + solve(It[[2]], St[[2]]) # looks more plausible than previous attempts

mean(rep(theta[1], LL[1]) + solve(It[[1]], St[[1]]))
mean(rep(theta[2], LL[2]) + solve(It[[2]], St[[2]]))

# NB how does theta look from specifying time-varying in glmmTMB call?

glmmTMB(Y.1~time+cont+bin+(1+time|id), dispformula = ~ 1 + time,
        family = nbinom2, data = data, control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3))) -> fit

getME(fit, 'beta')
cbind(1, 1:15) %*% fixef(fit)$disp


# V = mu(1+mu/theta)
# => theta = mu/(V/mu-1)

X1 <- do.call(rbind, lapply(Xk, '[[', 1)); beta1 <- beta[beta.inds[[1]]]
Z1 <- do.call(rbind, lapply(Zk, '[[', 1)); b1 <- do.call(rbind, lapply(b, function(x) x[b.inds[[1]]]))

Zb <- mapply(function(Z, b){
  Z[[1]] %*% b[b.inds[[1]]]
}, Z = Zk, b=b, SIMPLIFY = F)

mu1 <- X1 %*% beta1 + do.call(rbind, Zb)

V1 <- var(do.call(rbind, lapply(Yk, '[[', 1)))

sum(mu1)/(sum(c(V1)/mu1)-1)
