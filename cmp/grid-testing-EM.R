rm(list=ls())
source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin
disp.formula <- ~time

data <- simData_joint(n = 100, delta = c(-0.5, 0.1), ntms = 15, theta = c(-3, .25))
data <- data$data

test <- EM(long.formula, disp.formula, surv.formula, data, 
           control = list(verbose = T, summax.override = T, tol = 5e-2, gh.nodes = 25),
           summax=100)
plot.stepmat(test)
