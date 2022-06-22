rm(list=ls())
N <- 1e3; pete.flag <- T
source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin
disp.formula <- ~time

data <- simData_joint(n = 250, delta = c(-0.2, 0.1), ntms = 15, theta = c(-3, .25),
                      beta = c(.1,0,.2,-.1),D = matrix(c(0.25^2, 
                                                         0, 0, 0.05^2), 2, 2))
data <- data$data

dtest <- EM(long.formula, disp.formula, surv.formula, data, N = N,
           control = list(verbose = T, summax.override = T, tol = 1e-2, gh.nodes = 3, debug = T),
           summax=100)
plot.stepmat(dtest)
dtest$SE
dtest$coeffs$delta
