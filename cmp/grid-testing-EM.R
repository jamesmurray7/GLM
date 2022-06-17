rm(list=ls())
N <- 1e4
source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin
disp.formula <- ~time

data <- simData_joint(n = 250, delta = c(-0.5, 0.1), ntms = 15, theta = c(-3, .25))
data <- data$data

test <- EM(long.formula, disp.formula, surv.formula, data, N = N,
           control = list(verbose = T, summax.override = T, tol = 1e-2, gh.nodes = 3, debug = T),
           summax=100)
plot.stepmat(test)
test$SE
test$coeffs$delta
