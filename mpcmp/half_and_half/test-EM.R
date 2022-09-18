rm(list=ls())
source('EM.R')
test <- simData_joint(n = 250, delta = c(.8, -0.3), 
                      ntms = 10, theta = c(-2, .1), fup = 3,
                      beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                      D = matrix(c(0.25, 0, 0, 0), 2, 2))
data <- test$data

check.disps(test)

control <- list(verbose=T)
disp.control <- list(delta.method = 'bobyqa', min.profile.length = 3)
optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
summax=3

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin

fit <- EM(long.formula, disp.formula, surv.formula, data,
          control = list(verbose = T))
fit2 <- EM(long.formula, disp.formula, surv.formula, data,
           control = list(verbose = T),
           disp.control = list(re.maximise = F)) # See difference with/out remaximisation


# Intercept + slope -------------------------------------------------------
rm(list=ls())
source('EM.R')
test <- simData_joint2(n = 250, delta = c(.2,0), 
                       ntms = 10, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0.05), 2, 2))
data <- test$data

long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

fit <- EM(long.formula, disp.formula, surv.formula, data,
          control = list(verbose = T))
