#' Unique run-through/ideas testing.

# rm(list=ls())
library(survival)
library(Rcpp)
library(RcppArmadillo)
library(glmmTMB)
source('_Functions.R')
source('simData.R')
source('../inits.R')
sourceCpp('test.cpp') # CHANGE TO HYBRID.CPP WHEN FINALISED.
vech <- function(x) x[lower.tri(x, T)]

# Model spec.
long.formula <- Y ~ time + cont + bin + (1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

# Simulate some data, high-ish deltas
data <- simData_joint2(n = 250, delta = c(.8, 0), 
               ntms = 10, theta = c(-2, .1), fup = 3,
               beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
               D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data

#' ##############################################################
#' ##############################################################
#' ##############################################################

#' Parsing formula objects ----
formulas <- parseFormula(long.formula)
surv <- parseCoxph(surv.formula, data)
n <- surv$n

#' Initial conditions ----
inits.long <- Longit.inits(long.formula, disp.formula, data)
inits.surv <- TimeVarCox(data, inits.long$b, surv$ph, formulas)

# Longitudinal parameters
beta <- inits.long$beta.init
D <- inits.long$D.init
b <- lapply(1:n, function(i) inits.long$b[i, ])
delta <- inits.long$delta.init

# Survival parameters
zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
names(zeta) <- paste0('zeta_', names(zeta))
gamma <- inits.surv$inits[grepl('gamma', names(inits.surv$inits))]

#' Data objects ----
sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
dmats <- createDataMatrices(data, formulas, disp.formula)
X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
lY <- lapply(Y, lfactorial)
G <- dmats$G                             # Dispersion data matrix
m <- sapply(Y, length)
# survival
Fi <- sv$Fi; Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u; Delta <- surv$Delta 
l0 <- sv$l0
S <- sv$S; SS <- sv$SS

#' Parameter vector and list ----
Omega <- list(D=D, beta = beta, delta = delta, gamma = gamma, zeta = zeta)
params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
            beta, delta, gamma, zeta)

# Initial mus and nus ----
summax <- max(sapply(Y, max)) + 10
mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b)
nus <- mapply(function(G) exp(G %*% delta), G = G)

# Create grids -----
# Define grid dimensions.
mm <- do.call(c, mus); mu.min <- min(mm); mu.max <- max(mm)
nn <- do.call(c, nus); nu.vec <- as.vector(unique(nn))
all.mus <- generate_mus(max(0, mu.min - 5), mu.max + 5)

lambda.mat <- structure(gen_lambda_mat(max(0, mu.min - 5), mu.max + 5, nu.vec, summax),
                        dimnames = list(as.character(all.mus),
                                        as.character(nu.vec))
                        )

logZ.mat <- structure(gen_logZ_mat(max(0, mu.min - 5), mu.max + 5, nu.vec, lambda.mat, summax),
                      dimnames = list(as.character(all.mus),
                                      as.character(nu.vec))
)

V.mat <- structure(gen_V_mat(max(0, mu.min - 5), mu.max + 5, nu.vec, lambda.mat, logZ.mat, summax),
                   dimnames = list(as.character(all.mus),
                                   as.character(nu.vec))
)

#' Find b.hat and Sigma
b.hat <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
  optim(b, joint_density, joint_density_ddb,
        X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
        S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
        gamma = gamma, zeta = zeta, 
        lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, 
        all_mus = all.mus, all_nus = nu.vec,
        summax = summax, method = 'BFGS')$par
}, b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)

Sigma <- mapply(function(b, X, Y, lY, Z, G, S, SS, Fi, Fu, l0i, l0u, Delta){
  solve(joint_density_sdb(b = b, X = X, Y = Y, lY = lY, Z = Z, G = G, beta = beta, delta = delta, D = D,
                          S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, haz = l0u, Delta = Delta,
                          gamma = gamma, zeta = zeta, lambdamat = lambda.mat, Vmat = V.mat, logZmat = logZ.mat, 
                          all_mus = all.mus, all_nus = nu.vec, summax = summax, eps = .001))
}, b = b.hat, X = X, Y = Y, lY = lY, Z = Z, G = G, S = S, SS = SS, Fi = Fi, Fu = Fu,
l0i = l0i, l0u = l0u, Delta = Delta, SIMPLIFY = F)

any(unlist(Sigma) <= 0)

# mus, nus, tau
mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b.hat, SIMPLIFY = F)
nus <- mapply(function(G) exp(G %*% delta), G = G, SIMPLIFY = F)
tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), S = Sigma, Z = Z, SIMPLIFY = F)

# Define NEXT iter grid dimensions.
mu.min.old <- mu.min; mu.max.old <- mu.max; all.mus.old <- all.mus
mm <- do.call(c, mus); mu.min <- min(mm); mu.max <- max(mm)
nn <- do.call(c, nus); nu.vec <- as.vector(unique(nn))
all.mus <- generate_mus(max(0, mu.min - 5), mu.max + 5)

# Could add few lines to dictate whether matrices need to be recalculated?
subseteq <- ifelse(mu.max < mu.max.old && mu.min > mu.min.old, T, F)

# Indices for lookup given new b.hat
m <- mapply(function(a) get_indices(a, all.mus), a = mus, SIMPLIFY = F)
n <- mapply(function(a) get_indices(a, nu.vec), a = nus, SIMPLIFY = F)

# lambdas, Vs for \beta update.
lambdas <- mapply(function(m, n) mat_lookup(m, n, lambda.mat), m = m, n = n, SIMPLIFY = F)
Vs <- mapply(function(m, n) mat_lookup(m, n, V.mat), m = m, n = n, SIMPLIFY = F)


## E-step
gh <- statmod::gauss.quad.prob(3, 'normal')
w <- gh$w; v <- gh$n

# D
D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
# \beta
Sb <- mapply(Sbeta, X, Y, mus, nus, lambdas, Vs, SIMPLIFY = F)
Hb <- mapply(getW1, X, mus, nus, lambdas, Vs, SIMPLIFY = F)
# \delta
Sd <- mapply(function(G, b, X, Z, Y, lY, tau){
  Sdelta_cdiff(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps=.Machine$double.eps^(1/3))
}, G = G, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)
Hd <- mapply(function(G, b, X, Z, Y, lY, tau){
  Hdelta(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax, eps=.Machine$double.eps^(1/4))
}, G = G, b = b.hat, X = X, Z = Z, Y = Y, lY = lY, tau = tau, SIMPLIFY = F)
# \gamma
Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
  Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 1e-3)
}, b = b.hat, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)
Hgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
  Hgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 1e-4)
}, b = b.hat, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
## M-step
(D.new <- Reduce('+', D.update)/n)
(beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb)))
(delta.new <- delta - solve(Reduce('+', Hd), Reduce('+', Sd)))
(gammazeta.new <- c(gamma, zeta) - solve(Reduce('+', Hgz), rowSums(Sgz)))
lambda.update <- lambdaUpdate(sv$surv.times, sv$ft.mat, gamma, zeta, S, Sigma, b.hat, w, v)
l0.new <- sv$nev/rowSums(lambda.update)
l0u.new <- lapply(l0u, function(ll){
  l0.new[1:length(ll)]
})
l0i.new <- l0.new[match(sv$Ti, sv$ft)] 
l0i.new[is.na(l0i.new)] <- 0


# return...
list(
  #blah blah, 
  D = D.new, beta = beta.new, delta = delta.new,       # <Y>
  gamma = gammazeta.new[1], zeta = gammazeta.new[-1],  # Survival
  l0 = l0.new, l0u = l0u.new, l0i = as.list(l0i.new),  # ---""---
  b = b.hat, mus = mus, subseteq = subseteq
)
# Check if possible mu values is a subset of all.mus.old


