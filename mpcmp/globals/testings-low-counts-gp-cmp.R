rm(list=ls())
#' Case study: Underdispersed low counts...
# setwd('~/Documents/GLMM/mpcmp/globals')
source('EM.R')

test1 <- simData_joint2(delta = c(1.35, 0.0), D = matrix(c(.2, 0,0,0),2,2))$data
quantile(with(test1, tapply(Y, id, var))/with(test1, tapply(Y, id, mean)),na.rm=T)

# Increase counts a little...
test.beta <- c(0.5, 0.05, 0.1, -0.1)
test2 <- simData_joint2(delta = c(1.35, 0.0), beta = test.beta, D = matrix(c(.2, 0,0,0),2,2))$data
quantile(with(test2, tapply(Y, id, var))/with(test2, tapply(Y, id, mean)),na.rm=T)

# How do these profiles look?
par(mfrow=c(1,2))
hist(test1$Y, breaks = 5)  # roughly in [0,5]
hist(test2$Y, breaks = 10) # roughly in [0,10)
par(mfrow=c(1,1))

a <- glmmTMB(Y~ time+ cont + bin + (1|id), test1, poisson)

# 'Standard' // Genpois inits dont appear to work too well for lower counts...
fit1a <- EM(long.formula = Y~ time+ cont + bin + (1|id),
            disp.formula = ~1,
            surv.formula = Surv(survtime, status) ~ cont + bin,
            data = test1, control = list(verbose = T,
                                         min.profile.length=1), # note this 
            optim.control = list(Hessian = 'grad', eps = 1e-3), 
            individual.summax = T, min.summax = 30)

# Place all on lower summax -> can we speed up fit?
fit1b <- EM(long.formula = Y~ time+ cont + bin + (1|id),
            disp.formula = ~1,
            surv.formula = Surv(survtime, status) ~ cont + bin,
            data = test1, control = list(verbose = T),
            optim.control = list(Hessian = 'grad', eps = 1e-3), 
            individual.summax = T,
            summax.fn = max, min.summax = 10)

# And same for set test2
b <- glmmTMB(Y~time+cont+bin+(1|id), test2, genpois)

fit2a <- EM(long.formula = Y~ time+ cont + bin + (1|id),
            disp.formula = ~1,
            surv.formula = Surv(survtime, status) ~ cont + bin,
            genpois.inits = F,
            data = test2, control = list(verbose = T, min.profile.length = 1),
            optim.control = list(Hessian = 'grad', eps = 1e-3), 
            individual.summax = T)

# Start with genpois delta init
fit2b <- EM(long.formula = Y~ time+ cont + bin + (1|id),
            disp.formula = ~1,
            surv.formula = Surv(survtime, status) ~ cont + bin,
            data = test2, control = list(verbose = T, min.profile.length = 1),
            genpois.inits = T,
            optim.control = list(Hessian = 'grad', eps = 1e-3), 
            individual.summax = T)

# Didnt really improve; start with poisson and reduce summax as before...
fit2c <- EM(long.formula = Y~ time+ cont + bin + (1|id),
            disp.formula = ~1,
            surv.formula = Surv(survtime, status) ~ cont + bin,
            data = test2, control = list(verbose = T, min.profile.length = 1),
            genpois.inits = T,
            optim.control = list(Hessian = 'grad', eps = 1e-3), 
            individual.summax = T, summax.fn = max, min.summax = 10)

# Compare with GP1
setwd('../../genpois/')
source('EM.R')

# On lower counts  -- 
# 'Standard' // Genpois inits dont appear to work too well for lower counts...
fit1g <- EM(long.formula = Y ~ time+ cont + bin + (1|id),
            disp.formula = ~1,
            surv.formula = Surv(survtime, status) ~ cont + bin,
            data = test1, control = list(verbose = T, debug = T), genpois.inits = T) # Fails!


fit2g <- EM(long.formula = Y ~ time+ cont + bin + (1|id),
            disp.formula = ~1,
            surv.formula = Surv(survtime, status) ~ cont + bin,
            data = test2, control = list(verbose = T, debug = T), genpois.inits = T) # Runs decently well.

#my.summary(fit2g)
# Y (Generalised Poisson)
# Estimate    SE   2.5%  97.5% p-value
# beta_(Intercept)    0.520 0.043  0.436  0.604   0.000
# beta_time           0.048 0.006  0.035  0.060   0.000
# beta_cont           0.075 0.031  0.014  0.135   0.015
# beta_bin           -0.123 0.060 -0.241 -0.005   0.041
# phi_(Intercept)    -0.396 0.000 -0.396 -0.395   0.000
# 
# Event-time sub-model: 
#   Estimate    SE   2.5%  97.5% p-value
# zeta_cont    0.030 0.111 -0.188  0.248   0.789
# zeta_bin    -0.527 0.171 -0.863 -0.191   0.002
# gamma        0.246 0.268 -0.279  0.771   0.359
# 
# Elapsed times as follows (3 iterations):
#   startup time                EM time        post processing Total computation time 
# 3.111                  0.918                  0.300                  4.436 
