rm(list=ls())
source('EM.R')
set.seed(8)
test <- simData_joint2(n = 250, delta.int = c(0.3, -0.3),
                       ntms = 15, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0), 2, 2))
set.seed(8)
stest <- simData_joint2_safe(n = 250, delta.int = c(1.0, -1.0), delta.time = c(0.25, 0.25),
                       ntms = 15, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0), 2, 2), diff.tol = 5)
data <- stest$data

# check.disps(test) # NYI

control <- list(verbose=T, debug = T)
disp.control <- list(delta.method = 'optim', min.profile.length = 3,
                     truncated = T, max.val = 2.25)
optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
summax=3

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~time
update.deltas <- F

summax.fn <- NULL
min.summax <- 20
delta.update.quad <- T
beta.update.quad <- T
initialise.delta <- F

fit <- EM(long.formula, surv.formula, disp.formula, data,
          control = control,
          disp.control = disp.control,
          optim.control = optim.control,
          summax.fn = summax.fn, min.summax = min.summax,
          delta.update.quad = delta.update.quad,
          beta.update.quad = beta.update.quad,
          initialise.delta = initialise.delta)

my.summary(fit)

# Play around with different arguments...
fit2 <- EM(long.formula, surv.formula, data,
           control = control,
           disp.control = disp.control,
           optim.control = optim.control,
           summax.fn = summax.fn, min.summax = min.summax,
           delta.update.quad = F,
           beta.update.quad = F,
           initialise.delta = F)

my.summary(fit2)

fit3 <- EM(long.formula, surv.formula, data,
           control = control,
           disp.control = disp.control,
           optim.control = optim.control,
           summax.fn = summax.fn, min.summax = min.summax,
           delta.update.quad = T,
           beta.update.quad = F,
           initialise.delta = F)

# Sink + compare
sink('./three-competing-intonly.txt')
cat('\n---\nInitial conditions on delta, quadrature on both beta and delta.\n\n')
my.summary(fit)
cat('\n---\nNo initial conditions on delta, quadrature on neither beta nor delta.\n\n')
my.summary(fit2)
cat('\n---\nNo initial conditions on delta, quadrature on delta only.\n\n')
my.summary(fit3)
sink()


# Repeat same on intercept-and-slope --------------------------------------
rm(list=ls())
source('EM.R')
test <- simData_joint2(n = 250, delta = c(-1,1), 
                       ntms = 15, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0.05), 2, 2))
data <- test$data

check.disps(test)

control <- list(verbose=T, debug = T)
disp.control <- list(delta.method = 'optim', min.profile.length = 3,
                     truncated = T, max.val = 2.5)
optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
summax=3

long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ bin
update.deltas <- F

summax.fn <- function(y) max(y) + 10
min.summax <- 20
delta.update.quad <- T
beta.update.quad <- F
initialise.delta <- T

fit <- EM(long.formula, surv.formula, data,
          control = control,
          disp.control = disp.control,
          optim.control = optim.control,
          summax.fn = summax.fn, min.summax = min.summax,
          delta.update.quad = delta.update.quad,
          beta.update.quad = beta.update.quad,
          initialise.delta = initialise.delta)

my.summary(fit)

# Play around with different arguments...
fit2 <- EM(long.formula, surv.formula, data,
           control = control,
           disp.control = disp.control,
           optim.control = optim.control,
           summax.fn = summax.fn, min.summax = min.summax,
           delta.update.quad = F,
           beta.update.quad = F,
           initialise.delta = F)

my.summary(fit2)

fit3 <- EM(long.formula, surv.formula, data,
           control = control,
           disp.control = disp.control,
           optim.control = optim.control,
           summax.fn = summax.fn, min.summax = min.summax,
           delta.update.quad = T,
           beta.update.quad = F,
           initialise.delta = F)

# Sink + compare
sink('./three-competing-intslope.txt')
cat('\n---\nInitial conditions on delta, quadrature on both beta and delta.\n\n')
my.summary(fit)
cat('\n---\nNo initial conditions on delta, quadrature on neither beta nor delta.\n\n')
my.summary(fit2)
cat('\n---\nNo initial conditions on delta, quadrature on delta only.\n\n')
my.summary(fit3)
sink()
