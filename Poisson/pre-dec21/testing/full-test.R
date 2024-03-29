rm(list=ls())   
# Prerequisites -----------------------------------------------------------
vech <- function(x) x[lower.tri(x, diag = T)]
tr <- function(x) sum(diag(x))
an <- as.numeric
repCols <- function(X, n = 2) X[,rep(1:2, n)] 
repVec <- function(x, n = 2) rep(x, n)

library(survival)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('Poisson/ll.cpp')
sourceCpp('temp-gammaCalc.cpp')
source('Simulations/simData.R')
source('./DataFunctions/longFns.R')
source('./DataFunctions/survFns.R')

# Simulate some Data
args(simData)
beta <- rbind(c(0.75, -0.05, -0.2, 0.2),
              c(1, 0.05,  0.5, -0.80))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5

gamma <- c(-0.5, 0.8)
eta <- c(-0.3, 0.5)

data <- simData(250, 15, beta, D, gamma, eta, theta = c(-3, 0.2)) # appx. 50%
# Check all ids actually here, if not just rerun line above
length(unique(data$id))
ph <- coxph(Surv(survtime, status) ~ cont + bin, data = dplyr::distinct(data, id, cont, bin, survtime, status))

nK <- 2; q <- 2*nK
X <- getXi(data, 2)
Y <- getYi(data, 2)
Z <- getZi(data, 2)

#' Initial conditions -----
source('inits/inits.R')
inits.long <- Longit.inits(2, data)
inits.surv <- TimeVarCox(data, Ranefs(inits.long))

sv <- surv.mod(ph, data, inits.surv$l0.init)
# Extract all survival-related objects
ft <- sv$ft; nev <- sv$nev
surv.ids <- sv$surv.ids; surv.times <- sv$surv.times
Di <- sv$Di
l0 <- sv$l0; l0i <- sv$l0i; l0u <- sv$l0u
Fi <- sv$Fi; Fu <- sv$Fu
Fi.list <- lapply(1:nrow(Fi), function(i) Fi[i, ])
rvFi.list <- lapply(1:nrow(Fi), function(i) do.call(c, replicate(nK, Fi[i,], simplify = F)))
K <- getKi(data, ph); n <- 250
Krep <- sapply(1:n, function(x){  # For updates to \eta
  x <- apply(K[[x]], 2, rep, nrow(Fu[[x]]))
  if("numeric" %in% class(x)) x <- t(as.matrix(x))
  x
}) 

#' MVLME Step ----
source('./inits/MVLME.R')
mvlme.fit <- mvlme(data, Y, X, Z, inits.long, nK = 2, q = 4,  
                   mvlme.tol = 5e-3, verbose = F)

beta <- mvlme.fit$beta
D <- mvlme.fit$D
b <- mvlme.fit$b # currently a list
gamma <- inits.surv$inits[3:4]; gr <- rep(gamma, each = nK)
eta <- inits.surv$inits[1:2]

#' Quadrature ----
gh <- statmod::gauss.quad.prob(9, 'normal')
w <- gh$w; v <- gh$n

#' Start EM ----
EMtime <- proc.time()[3]

#' #################
#' E-step =========#
#' #################
b.hat <- mapply(function(b, X, Y, lfactY, Z, K, Delta, l0i, Fi, l0u, Fu, rvFi){
  ucminf::ucminf(
    b, ll, gradll, Y, lfactY, X, Z, D, K, Delta, l0i, Fi, l0u, Fu, 
    beta, eta, gr, rvFi, nK, q
  )$par}, 
  b = b, X = X, Y = Y, lfactY = lapply(Y, lfactorial), Z = Z, K = K, Delta = as.list(Di),
  l0i = as.list(l0i), Fi = Fi.list, l0u = l0u, Fu = Fu, rvFi = rvFi.list, SIMPLIFY = F)

Sigmai <- mapply(function(b, Y, lfactY, X, Z, K, l0u, Fu){
  solve(-1 * hessll(b, Y, lfactY, X, Z, D, K, l0u, Fu, beta, eta, gr, nK))
  }, b = b.hat, Y = Y, lfactY = lapply(Y, lfactorial), X = X, Z = Z, K = K,
  l0u = l0u, Fu = Fu, SIMPLIFY = F)

#' Steps to update Longitudinal and RE parameters --------
Sbeta <- mapply(function(Y, X, Z, b){      # Worth noting that Simplify = F is slightly faster
  crossprod(X, Y) - crossprod(X, exp(X %*% beta + Z %*% b))
}, b = b.hat, Y = Y, X = X, Z = Z, SIMPLIFY = F)

# Test with quadrature
tau.long <- mapply(function(S, Z){
  sqrt(diag(tcrossprod(Z %*% S, Z)))
}, S = Sigmai, Z = Z, SIMPLIFY = F)

Sbeta2 <- mapply(function(Y, X, Z, b, tau){
  rhs <- numeric(length(beta))
  for(l in 1:9) rhs <- rhs + w[l] * crossprod(X, exp(X %*% beta + Z %*% b + tau * v[l]))
  crossprod(X, Y) - rhs
}, Y = Y, X = X, Z = Z, b = b.hat, tau = tau.long, SIMPLIFY = F)

Ibeta <- mapply(function(X, Z, b){
  crossprod(-diag(c(exp(X %*% beta + Z %*% b))) %*% X, X)
}, X = X, Z = Z, b = b.hat, SIMPLIFY = F)

Ibeta2 <- mapply(function(X, Z, b, tau){
  out <- matrix(0, length(beta), length(beta))
  for(l in 1:9) out <- out + -w[l] * crossprod(diag(c(exp(X %*% beta + Z %*% b + tau * v[l]))) %*% X, X)
  out
}, X = X, Z = Z, b = b.hat, tau = tau.long, SIMPLIFY = F)


#' Steps to update D -----------
D.new <- mapply(function(S, b){
  S + tcrossprod(b)
}, S = Sigmai, b = b.hat, SIMPLIFY = F)

#' Steps to update survival parameters ------------
inds <- split(seq(nK * 2), rep(1:nK, each = 2))
b.hat.split <- lapply(b.hat, function(y) lapply(inds, function(x) y[x]))
# Define mu.surv
mu.surv <- mapply(function(K, Fu, b){
  rhs <- 0
  for(k in 1:nK) rhs <- rhs + gamma[k] * b[inds[[k]]]
  exp(K %*% eta + Fu %*% rhs)
}, K = Krep, Fu = Fu, b = b.hat, SIMPLIFY = F)

# Define tau objects
S <- lapply(Sigmai, function(y) lapply(inds, function(x) y[x,x]))

tau <- mapply(function(Fu, S){
  out <- numeric(nrow(Fu))
  for(k in 1:nK) out <- out + diag(gamma[k]^2 * tcrossprod(Fu %*% S[[k]], Fu))
  out
}, Fu = Fu, S = S, SIMPLIFY = F)

tau.tilde <- mapply(function(Fu, S){
  mat <- matrix(0, nr = nK, nc = nrow(Fu))
  for(k in 1:nK) mat[k, ] <- diag(tcrossprod(Fu %*% S[[k]], Fu))
  mat
}, Fu = Fu, S = S, SIMPLIFY = F)

tau.surv <- lapply(tau, sqrt)
tau2.surv <- lapply(tau, function(x){
  x <- x^(-0.5)
  if(any(is.infinite(x))) x[which(is.infinite(x))] <- 0   # Avoid NaN
  x
})

#' S(\gamma) ----
gh.nodes <- 9
Sgamma <- mapply(function(Delta, tau.surv, mu.surv, l0u, Fu, Fi, b){
  t(Sgammacalc(gamma, Delta, tau.surv, mu.surv, l0u, Fu, Fi, w, v, b, nK, gh.nodes))
}, Delta = as.list(Di), tau.surv = tau.surv, mu.surv = mu.surv, l0u = l0u,
   Fu = Fu, Fi = Fi.list, b = b.hat.split, SIMPLIFY = F)

#' I(\gamma)
Igamma <- mapply(function(tau.tilde, tau.surv, tau2.surv, mu.surv, Fu, l0u, b){
  gamma2Calc(gamma, tau.tilde, tau.surv, tau2.surv, mu.surv, w, v, Fu, l0u, b, nK, gh.nodes)
}, tau.tilde = tau.tilde, tau.surv = tau.surv, tau2.surv = tau2.surv,
   mu.surv = mu.surv, Fu = Fu, l0u = l0u, b = b.hat.split, SIMPLIFY = F)
  
#' S(\eta) ----
Seta <- mapply(function(K, KK, Delta, l0u, mu.surv, tau.surv){
 cen <- Delta %*% K
 rhs <- c(0, 0)
 for(l in 1:gh.nodes) rhs <- rhs + w[l] * t(KK) %*% (l0u * (mu.surv * exp(tau.surv * v[l])))
 cen-t(rhs)
}, K = K, KK = Krep, Delta = as.list(Di), l0u = l0u, mu.surv = mu.surv, tau.surv = tau.surv, SIMPLIFY = F)

#' I(\eta) ----
Ieta <- mapply(function(K, KK, tau.surv, mu.surv, l0u){
  Ietacalc(2, K, KK, tau.surv, mu.surv, l0u, w, v, gh.nodes)
}, K = K, KK = Krep, tau.surv = tau.surv, mu.surv = mu.surv, l0u = l0u, SIMPLIFY = F)

#' Second derivates of \gamma and \eta ('cross-terms') -----

Igammaeta <- list() # for some reason I can't get mapply() to work with this 
for(i in 1:n){
  Igammaeta[[i]] <- Igammaetacalc(2, Krep[[i]], tau.surv[[i]], tau.tilde[[i]], mu.surv[[i]], l0u[[i]], Fu[[i]], 
                                  b.hat.split[[i]], gamma, w, v, nK, gh.nodes)
}


#' #################
#' M-step =========#
#' #################

#' \beta -----
beta.new <- beta + solve(-Reduce('+', Ibeta), 
                         rowSums(do.call(cbind, Sbeta)))

beta.new2 <- beta + solve(-Reduce('+', Ibeta2), 
                          rowSums(do.call(cbind, Sbeta2)))

#' D -----
D.new <- Reduce('+', D.new)/n

#' The baseline hazard, \lambda_0
lambda <- lambdaUpdate(surv.times, ft, gamma, eta, K, S,
                       b.hat.split, n, w, v, gh.nodes)

l0.new <- nev/rowSums(lambda)
l0u.new <- lapply(l0u, function(x){
  ll <- length(x); l0.new[1:ll]
})
l0i.new <- c()
l0i.new[which(Di == 0)] <- 0 
l0i.new[which(Di == 1)] <- l0.new[match(Fi[which(Di==1), 2], ft)]


#' (\gamma, \eta) ----
# Set up score vector and information matrix
Sge <- c(colSums(do.call(rbind, Sgamma)), colSums(do.call(rbind, Seta)))
Imat <- as.matrix(Matrix::bdiag(Reduce('+', Igamma),
                                Reduce('+', Ieta)))
# Fill in off-block diagonal with Igammaeta
eta.inds <- (nK+1):ncol(Imat)
for(k in 1:nK){
  Imat[k, eta.inds] <- Imat[eta.inds, k] <- rowSums(do.call(cbind, lapply(Igammaeta, '[[', k)))
}

gamma.eta.new <- c(gamma, eta) + solve(Imat, Sge)

#' And then update + print etc. ----

#' ####
#' Room for misc. testing 
#' ####
sourceCpp('temp-gammaCalc.cpp')

gammaeta_ll(c(gamma, eta), Di[1], K[[1]], Krep[[1]], gr, rvFi.list[[1]], b.hat[[1]], l0u[[1]],
            Fu[[1]], b.hat.split[[1]], S[[1]], nK, w, v, gh.nodes)

imat <- list()
for (i in 1:n){
  imat[[i]] <- -1 * numDeriv::hessian(gammaeta_ll, c(gamma, eta), method="Richardson",method.args=list(),
                                      Di[i], K[[i]], Krep[[i]], gr, rvFi.list[[i]], b.hat[[i]], l0u[[i]],
                                      Fu[[i]], b.hat.split[[i]], S[[i]], nK, w, v, gh.nodes)
}

Reduce('+', imat)  #  Way to test all



