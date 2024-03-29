# Test one ----------------------------------------------------------------
rm(list=ls())
source('EM.R')
test <- simData_joint2(n = 250, delta = c(0.3, 0), 
                      ntms = 10, theta = c(-2, .1), fup = 3,
                      beta = c(0.0, 0.5, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                      D = matrix(c(0.25, 0, 0, 0.00), 2, 2)) ## Doesn't appear to work particularlily well; instead can simply save data from non-hybrid dir.

# Testing dispersion of outputted 'Y'.
summary(with(test$data, tapply(Y, id, var))/with(test$data, tapply(Y, id, mean)))

range(test$data$Y)
data <- test$data
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

fit <- EM(long.formula, disp.formula, surv.formula, data, summax = 50, 
          control = list(verbose = T, auto.summax = F, tol = 1e-2, gh.nodes = 3, debug = T))

plot.stepmat(fit)


# Intercept with bigger beta values ---------------------------------------
rm(list=ls())
source('EM.R')
test <- simData_joint(n = 250, delta = c(.4,0), 
                       ntms = 10, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0.00), 2, 2)) ## Doesn't appear to work particularlily well; instead can simply save data from non-hybrid dir.
fit <- glmmTMB(Y~time+cont+bin+(1|id),test$data, family = poisson) 
fit
VarCorr(fit)$cond$id
fixef(fit)$cond
range(test$data$Y)
data <- test$data

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

fit <- EM(long.formula, disp.formula, surv.formula, data,
          control = list(verbose = T), grid.summax = 'same')
fit$elapsed
plot.stepmat(fit)
plot.grid.and.step(fit)

# Let first grids be created on different truncation amount.  
fit2 <- EM(long.formula, disp.formula, surv.formula, data, summax = 100,
            control = list(verbose = T), grid.summax = 'a')
fit2$elapsed
plot.stepmat(fit2)
plot.grid.and.step(fit2)