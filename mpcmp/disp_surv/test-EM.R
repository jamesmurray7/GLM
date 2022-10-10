rm(list=ls())
source('EM.R')
source('_zzz.R')
stest <- simData_joint2(n = 250, delta.int = c(1.0, -1.0), delta.time = c(0.2, -0.1),
                       ntms = 15, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0), 2, 2), diff.tol = 6,
                       gamma.disp = -0.3)
data <- stest$data
disp.formula <- ~time
# check.disps(test) # NYI

fit <- EM(long.formula, surv.formula, disp.formula, data,
          control = control,
          disp.control = disp.control,
          optim.control = optim.control,
          summax.fn = summax.fn, min.summax = min.summax,
          delta.update.quad = delta.update.quad,
          beta.update.quad = beta.update.quad,
          initialise.delta = initialise.delta)

my.summary(fit)

