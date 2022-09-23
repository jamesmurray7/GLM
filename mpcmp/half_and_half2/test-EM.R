rm(list=ls())
source('EM.R')
test <- simData_joint2(n = 250, delta = c(1, -1.), 
                       ntms = 15, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0), 2, 2))
data <- test$data

check.disps(test)

control <- list(verbose=T)
disp.control <- list(delta.method = 'optim', min.profile.length = 3,
                     max.val = 10)
optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
summax=3

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
update.deltas <- F

summax.fn <- function(y) max(y) + 10
min.summax <- 20
delta.update.interval <- 1L
update.deltas <- T

fit <- EM(long.formula, surv.formula, data,
          control = control,
          disp.control = disp.control,
          optim.control = optim.control,
          summax.fn = summax.fn, min.summax = min.summax,
          update.deltas = update.deltas, delta.update.interval = delta.update.interval)


# Intercept + slope -------------------------------------------------------
rm(list=ls())
source('EM.R')
test <- simData_joint2(n = 250, delta = c(-1,1), 
                       ntms = 15, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0.05), 2, 2))
data <- test$data

long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1
control <- list(verbose=T)
summax.fn <- function(y) max(y) + 10
min.summax <- 20
disp.control <- list(delta.method = 'bobyqa', min.profile.length = 4,
                     max.val = 10)
optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)

long.formula <- Y~time+cont+bin+(1 + time|id)
surv.formula <- Surv(survtime, status) ~ bin
update.deltas <- F

fit2 <- EM(long.formula, surv.formula, data,
           control = control,
           disp.control = disp.control,
           optim.control = optim.control,
           summax.fn = summax.fn, min.summax = min.summax)
