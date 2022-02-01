#' ####
#' Testing quadratic setup
#' ####

setwd('~/Documents/PhD/GLM/quad/')
source('simData.R')
 
vech <- function(x) x[lower.tri(x, diag = T)]
beta <- rbind(c(1, -0.2, 0.01, 0.33, -0.50),
              c(0, -0.5, 0.05, -0.33, 0.50),
              c(3, 0.1, -0.05, 0.5, 0.1))

D <- as.matrix(Matrix::bdiag(replicate(3, diag(c(0.5^2, .2^2, .05^2)), simplify = F)))

data <- quad.simData(n = 250, ntms = 10, beta = beta, D = D)

# plot to see if it looks decently quadratic...
data$dat %>%
  select(id, time, Y.1:Y.3) %>%
  pivot_longer(Y.1:Y.3) %>%
  ggplot(aes(x = time, y = value, group = id)) +
  geom_line(alpha = .33) +
  facet_wrap(~name, scales = 'free') + 
  geom_smooth(method = 'lm', formula = y ~ poly(x, 2, raw = T), se = F,
              mapping = aes(group = NULL))

# Fit quickly using nlme as usual...
nlme::lme(fixed = Y.1 ~ time + I(time^2) + cont + bin,
          random = ~ time + I(time^2) | id,
          data = data$dat,
          method = "ML",
          control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

# Data matrices -----------------------------------------------------------
dat <- data$dat
sda <- data$survdat
ph <- coxph(Surv(survtime, status) ~ cont + bin, sda)

n <- nrow(sda)
# Get data matrices
m <- Y <- X <- Z <- K <- list()
for(i in 1:n){
  i.dat <- dat[dat$id == i, ]
  m[[i]] <- rep(nrow(i.dat), 3)
  Y[[i]] <- c(i.dat$Y.1, i.dat$Y.2, i.dat$Y.3)
  X[[i]] <- structure(model.matrix(~time + I(time^2) + cont + bin, i.dat),
                      dimnames = list(as.character(1:nrow(i.dat)),
                                      c('(Intercept)', 'time', 'time^2', 'cont', 'bin')))
  Z[[i]] <- structure(model.matrix(~time + I(time^2), i.dat),
                      dimnames = list(as.character(1:nrow(i.dat)),
                                      c('(Intercept)', 'time', 'time^2')))
  K[[i]] <- unname(cbind(unique(i.dat$cont), unique(i.dat$bin)))
}

# Block matrices
Xblock <- lapply(X, function(x) as.matrix(Matrix::bdiag(replicate(3, x, simplify = F))))
Zblock <- lapply(Z, function(x) as.matrix(Matrix::bdiag(replicate(3, x, simplify = F))))

inits.long <- Longit.inits(3, dat)
b <- Ranefs(inits.long)
beta <- inits.long$beta.init
var.e <- inits.long$var.e.init
V <- lapply(m, function(iii) {
  diag(x = rep(var.e, iii), ncol = sum(iii))
})
D <- inits.long$D.init
inits.surv <- TimeVarCox(dat, b)
b <- lapply(1:n, function(i) b[i, -10])
rm(inits.long) # large object

# Survival objects
sv <- surv.mod(ph, dat, inits.surv$l0.init)
Delta <- as.list(sv$Di)
l0i <- as.list(sv$l0i)
l0u <- sv$l0u
Fi <- lapply(1:n, function(i) do.call(c, replicate(3, sv$Fi[i, ], simplify = F)))
Fu <- sv$Fu
KK <- sapply(1:n, function(i){
  x <- apply(K[[i]], 2, rep, nrow(Fu[[i]]))
  if('numeric'%in%class(x)) x <- t(as.matrix(x))
  x
})
gamma <- inits.surv$inits[3:5]
eta <- inits.surv$inits[1:2]

# Quadrature //
gh <- 3
aa <- statmod::gauss.quad.prob(gh, 'normal')
w <- aa$w; v <- aa$n

vD <- vech(D);
names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
params <- c(vD, beta, var.e, gamma, eta)

b.hat <- mapply(function(b, X, Y, Z, V, Delta, K, Fi, l0i, KK, Fu, l0u){
  ucminf::ucminf(b, joint_density, joint_density_db,
                 X = X, Y = Y, Z = Z, beta = beta, V = V, D = D,
                 K = K, KK = KK, Fi = Fi, Fu = Fu, l0u = l0u, l0i = l0i,
                 Delta = Delta, gamma = rep(gamma, each = 3), eta = eta)$par
}, b = b, X = Xblock, Y = Y, Z = Zblock, V = V, Delta = Delta, K = K, Fi = Fi, l0i = l0i,
KK = KK, Fu = Fu, l0u = l0u, SIMPLIFY = F)

bmat <- lapply(b.hat, matrix, nc = 2, byr = T)
bsplit <- lapply(b.hat, function(y) lapply(split(seq(6), rep(seq(3), each = 2)), function(x) y[x]))


