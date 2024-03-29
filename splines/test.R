#' ####
#' Testing spline model fit
#' ####

# setwd('~/Documents/GLMM/splines/')
source('simData.R')
source('survFns.R')

# Stick with quadratic for means of testing 
vech <- function(x) x[lower.tri(x, diag = T)]
beta <- rbind(c(1, -0.2, 0.01, 0.33, -0.50),
              c(0, -0.5, 0.05, -0.33, 0.50),
              c(3, 0.1, -0.05, 0.5, 0.1))

D <- as.matrix(Matrix::bdiag(replicate(3, diag(c(0.5^2, .2^2, .05^2)), simplify = F)))

data <- poly.simData(n = 250, ntms = 10, beta = beta, D = D, theta0 = -4.5)

# plot to see if it looks decently quadratic...
# library(tidyverse)
data$dat %>%
  select(id, time, Y.1:Y.3) %>%
  pivot_longer(Y.1:Y.3) %>%
  ggplot(aes(x = time, y = value, group = id)) +
  geom_line(alpha = .33) +
  facet_wrap(~name, scales = 'free') + 
  geom_smooth(method = 'lm', formula = y ~ poly(x, 2, raw = T), se = T, lwd = .5,
              mapping = aes(group = NULL))

# Fit quickly using nlme as usual...
nlme::lme(fixed = Y.1 ~ time + I(time^2) + cont + bin,
          random = ~ time + I(time^2) | id,
          data = data$dat,
          method = "ML",
          control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

nlme::lme(fixed = Y.1 ~ bs(time, degree = 2) + cont + bin,
          random = ~ bs(time, degree = 2) | id,
          data = data$dat,
          method = "ML",
          control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

# Data matrices -----------------------------------------------------------
dat <- data$dat
sda <- data$survdat
degree <- 6
library(survival)
ph <- coxph(Surv(survtime, status) ~ cont + bin, sda)

basis <- getsurvbasis(dat, degree)

n <- nrow(sda)
# Get data matrices
nK <- 3
Z <- getZfromsurvbasis(basis, dat)
m <- Y <- Ymat <- X <- K <- list()
for(i in 1:n){
  i.dat <- dat[dat$id == i, ]
  m[[i]] <- rep(nrow(i.dat), 3)
  Y[[i]] <- c(i.dat$Y.1, i.dat$Y.2, i.dat$Y.3)
  Ymat[[i]] <- cbind(i.dat$Y.1, i.dat$Y.2, i.dat$Y.3)
  colnames(Ymat[[i]]) <- paste0('Y.', 1:nK)
  X[[i]] <- cbind(Z[[i]], i.dat$cont, i.dat$bin)
  colnames(X[[i]])[(degree + 2):ncol(X[[i]])] <- c('cont', 'bin')
  K[[i]] <- unname(cbind(unique(i.dat$cont), unique(i.dat$bin)))
}

# Block matrices
Xblock <- lapply(X, function(x) as.matrix(Matrix::bdiag(replicate(nK, x, simplify = F))))
Zblock <- lapply(Z, function(x) as.matrix(Matrix::bdiag(replicate(nK, x, simplify = F))))
XtX <- lapply(Xblock, crossprod)

df.basis <- as.data.frame(cbind(id = do.call(c, lapply(1:n, function(i) rep(i, m[[i]][1]))), do.call(rbind, Ymat), do.call(rbind, X)))

lme(Y.1 ~ basis1 + basis2 + basis3 + basis4 + basis5 + basis6 + cont + bin,
    random = ~ basis1 + basis2 + basis3 + basis4 + basis5 + basis6|id,
    data = df.basis,
    method = 'ML',
    control = nlme::lmeControl(opt = "optim", msTol = 1e-3)) -> manual.basislme

lme(Y.1 ~ bs(time, degree = 3) + cont + bin,
    random = ~ bs(time, degree = 3) | id,
    data = dat,
    method = 'ML',
    control = nlme::lmeControl(opt = "optim", msTol = 1e-3)) -> auto.basislme


source('inits.R')
inits.long <- Longit.inits(nK, df.basis, degree = 3, basis = 'survival')
b <- Ranefs(inits.long)
beta <- inits.long$beta.init
var.e <- inits.long$var.e.init
V <- lapply(m, function(iii) {
  diag(x = rep(var.e, iii), ncol = sum(iii))
})
D <- inits.long$D.init
inits.surv <- TimeVarCox2(dat, b, basis, nK = 3)
b <- lapply(1:n, function(i) b[i, ])

# Survival objects
source('survFns.R')
sv <- surv.mod(ph, dat, inits.surv$l0.init)
Delta <- as.list(sv$Di)
l0i <- as.list(sv$l0i)
l0u <- sv$l0u
Fi <- lapply(1:n, function(i) do.call(c, replicate(nK, sv$Fi[i, ], simplify = F)))
Fu <- sv$Fu
KK <- sapply(1:n, function(i){
  x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
  if('numeric'%in%class(x)) x <- t(as.matrix(x))
  x
})
gamma <- inits.surv$inits[3:length(inits.surv$inits)]
eta <- inits.surv$inits[1:2]

# Quadrature //
gh <- 3
aa <- statmod::gauss.quad.prob(gh, 'normal')
w <- aa$w; v <- aa$n

vD <- vech(D);
names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
params <- c(vD, beta, var.e, gamma, eta)

# Split beta -> do this at end of every M-step?
beta.split <- lapply(1:nK, function(k){
  ii <- grepl(paste0('^beta', k, '_'), names(beta))
  beta[ii]
})
# or get inds once
beta.inds <- lapply(1:nK, function(k){
  which(grepl(paste0('^beta', k, '_'), names(beta)))
}) # this clearly pref. option.


library(Rcpp)
library(RcppArmadillo)
sourceCpp('quad.cpp')

#' ###
#' E-step
#' ###

b.hat <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, l0u){
  ucminf::ucminf(b, joint_density, joint_density_db,
                 X = X, Y = Y, Z = Z, beta = beta, V = V, D = D,
                 K = K, KK = KK, Fi = Fi, Fu = Fu, l0u = l0u, l0i = l0i,
                 Delta = Delta, gamma = rep(gamma, each = (degree + 1)), eta = eta)$par
}, b = b, X = Xblock, Y = Y, Z = Zblock, V = V, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)

# Split out into response-oriented objects.
bmat <- lapply(b.hat, matrix, nc = (degree + 1), byr = T)
bsplit <- lapply(b.hat, function(y) lapply(split(seq(nK * (degree + 1)), rep(seq(nK), each = (degree + 1))), function(x) y[x]))

# Posterior variance
Sigmai <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, l0u){
 solve(joint_density_sdb(b=b, X = X, Y = Y, Z = Z, beta = beta, V = V, D = D,
                   K = K, KK = KK, Fi = Fi, Fu = Fu, l0u = l0u, l0i = l0i,
                   Delta = Delta, gamma = rep(gamma, each = (degree + 1)), eta = eta))
}, b = b.hat, X = Xblock, Y = Y, Z = Zblock, V = V, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)

do.call(rbind, lapply(Sigmai, diag)) # check for any negative variances

# Split into nK constituent sub-matrices
S <- lapply(Sigmai, function(y) lapply(split(seq(nK * (degree + 1)), rep(seq(nK), each = (degree + 1))), function(x) y[x, x]))

#' Step to update D -----
Drhs <- mapply(function(S, b) S + tcrossprod(b), S  = Sigmai, b = b.hat, SIMPLIFY = F)

#' Step to update beta ------
tau.long <- mapply(function(S, Z) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigmai, Z = Zblock, SIMPLIFY = F)

beta.rhs <- mapply(function(b, X, Y, Z, tau){
  mu <- Z %*% b
  rhs <- numeric(length(tau))
  for(l in 1:gh) rhs <- rhs + w[l] * tau * v[l]
  crossprod(X, Y - mu - rhs)
}, b = b.hat, X = Xblock, Y = Y, Z = Zblock, tau = tau.long, SIMPLIFY = F)

#' Step to update sigma^2 ------
tau.longK <- mapply(function(S, Z){
  out <- vector('list', nK)
  for(k in 1:nK) out[[k]] <- sqrt(diag(tcrossprod(Z %*% S[[k]], Z)))   # exploiting set-up of Z being the same
  out
}, S = S, Z = Z, SIMPLIFY = F)

mu.longK <- mapply(function(X, Z, b){
  out <- list()
  for(k in 1:nK) out[[k]] <- X %*% beta[beta.inds[[k]]] + Z %*% b[[k]]
  out
}, X = X, Z = Z, b = bsplit , SIMPLIFY = F)

Esigma <- mapply(function(Y, mu, tau){
  temp <- matrix(NA, nr = gh, nc = nK)
  for(k in 1:nK){
    for(l in 1:gh){
      temp[l, k] <- w[l] * crossprod(Y[, k] - mu[[k]] - tau[[k]] * v[l])
    }
  }
  colSums(temp)
}, Y = Ymat, mu = mu.longK, tau = tau.longK)

#' Step to update (gamma, eta)
Sge <- mapply(function(bmat, S, K, KK, Fu, Fi, l0u, Delta){
  Sgammaeta(c(gamma, eta), bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi[1:(degree + 1)], haz = l0u,
            Delta = Delta, w = w, v = v, eps = 1e-4)
}, bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)

Hge <- mapply(function(bmat, S, K, KK, Fu, Fi, l0u, Delta){
  Hgammaeta(c(gamma, eta), bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi[1:(degree + 1)], haz = l0u,
            Delta = Delta, w = w, v = v, eps = 1e-4)
}, bmat = bmat, S = S, K = K, KK = KK, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)

#' ###
#' M-step
#' ###

D.update <- Reduce('+', Drhs) / n
beta.update <- c(solve(Reduce('+', XtX)) %*% Reduce('+', beta.rhs))
var.e.update <- rowSums(Esigma) / colSums(do.call(rbind, m))
gammaeta.update <- c(gamma, eta) - solve(Reduce('+', Hge), rowSums(Sge))

# baseline hazard

# create a matrix which has basis for each possible failure time
basis.ft <- cbind(1, sv$basisdf[match(sv$ft, sv$basisdf[, 1]), -1])

lambda <- lambdaUpdate(sv$surv.times, sv$ft, basis.ft, gamma, eta, K, S, bsplit, w, v)
# Baseline hazard objects
l0.new <- sv$nev/rowSums(lambda)
l0u.new <- lapply(l0u, function(x){
  ll <- length(x); l0.new[1:ll]
})
l0i.new <- c()
l0i.new[which(unlist(Delta) == 0)] <- 0 
l0i.new[which(unlist(Delta) == 1)] <- l0.new[match(sda[which(unlist(Delta)==1), 'survtime'], sv$ft)]
l0i.new <- as.list(l0i.new)



