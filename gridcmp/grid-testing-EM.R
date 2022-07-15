# Test one ----------------------------------------------------------------
rm(list=ls())
source('grid-EM-nogrids.R')

test <- simData_joint2(n = 250, delta = c(.8,.1), 
                      ntms = 10, theta = c(-2, .1), fup = 3,
                      beta = c(0.0, 0.5, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                      D = matrix(c(0.25, 0, 0, 0.00), 2, 2))

# Testing dispersion of outputted 'Y'.
summary(with(test$data, tapply(Y, id, var))/with(test$data, tapply(Y, id, mean)))

range(test$data$Y)
data <- test$data
rm(to.remove)
N <- 1e4; pete.flag <- F
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

fit <- EM(long.formula, disp.formula, surv.formula, data, N = N,
          control = list(verbose = T, auto.summax = T, tol = 1e-2, gh.nodes = 3, debug = T))


# Intercept and slope now? ------------------------------------------------
rm(list=ls())
source('grid-EM-nogrids.R')
test <- simData_joint2(n = 250, delta = c(.8,0), 
                       ntms = 10, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0.05), 2, 2))


fit <- glmmTMB(Y~time+cont+bin+(1+time|id),test$data, family = poisson)
fit
VarCorr(fit)$cond$id
fixef(fit)$cond
range(test$data$Y)
data <- test$data

long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

fit <- EM(long.formula, disp.formula, surv.formula, data, N = 1e4,
          control = list(verbose = T, auto.summax = T, tol = 1e-2, gh.nodes = 3, debug = T))
