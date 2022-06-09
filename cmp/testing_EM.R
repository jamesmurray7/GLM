rm(list=ls())
source('EM.R')
long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin
disp.formula <- ~time

data <- simData_joint(n = 150)
data <- data$data

test <- EM(long.formula, disp.formula, surv.formula, data, control = list(verbose = T))
plot(test$stepmat)
