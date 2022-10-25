# genpois
rm(list=ls())
source('EM.R')
test <- simData_joint()
data <- test$data
source('_zzz.R')
disp.formula <- ~time

fit <- EM(long.formula, disp.formula, surv.formula, 
          data, control = control,
          genpois.inits = T)

# Intslope ----------------------------------------------------------------
rm(list=ls())
source('EM.R')
test <- simData_joint(D = matrix(c(.2, 0,0,.05), 2, 2), diff.tol = 99,
                      phi = c(0.5, -0.25),
                      beta = c(0.5, 0.05, 0.1, -0.1))
data <- test$data
source('_zzz.R')
long.formula <- Y~time+cont+bin+(1+time|id)
disp.formula <- ~time

fit <- EM(long.formula, disp.formula, surv.formula, 
          data, control = control,
          genpois.inits = T)

my.summary(fit)