rm(list=ls())
source('EM.R')
source('_zzz.R')
stest <- simData_joint2_safe(n = 250, delta.int = c(0.0, -0.5), delta.time = c(0.2, -0.2),
                       ntms = 15, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0), 2, 2), diff.tol = 5)
data <- stest$data

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


# Impose increasing (in either direction) variation over time -------------
rm(list=ls())
source('EM.R')
source('_zzz.R')
stest <- simData_joint2_safe(n = 250, delta.int = c(0.0, 0.0), delta.time = c(-0.2, -0.3),
                             ntms = 15, theta = c(-2, .1), fup = 3,
                             beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                             D = matrix(c(0.25, 0, 0, 0), 2, 2), diff.tol = 5)
data <- stest$data

# check.disps(test) # NYI

fit <- EM(long.formula, surv.formula, disp.formula, data,
          control = control,
          disp.control = disp.control,
          optim.control = optim.control,
          summax.fn = summax.fn, min.summax = min.summax,
          delta.update.quad = delta.update.quad,
          beta.update.quad = beta.update.quad,
          initialise.delta = initialise.delta)

# Comparison w/ poisson / nbin
fit$logLik
a <- glmmTMB(Y ~ time + cont + bin + (1|id), data, poisson) # Not particularly compelling!
b <- glmmTMB(Y ~ time + cont + bin + (1|id), data, nbinom2, disp = ~time)



# Underdispersion ---------------------------------------------------------
rm(list=ls())
source('EM.R')
source('_zzz.R')
stest <- simData_joint2_safe(n = 250, delta.int = c(2.0, 2.0), delta.time = c(0, -0.4),
                             ntms = 15, theta = c(-2, .1), fup = 3,
                             beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                             D = matrix(c(0.25, 0, 0, 0), 2, 2), diff.tol = 5)
data <- stest$data

# check.disps(test) # NYI
disp.formula <- ~time

fit <- EM(long.formula, surv.formula, disp.formula, data,
          control = control,
          disp.control = disp.control,
          optim.control = optim.control,
          summax.fn = summax.fn, min.summax = min.summax,
          delta.update.quad = delta.update.quad,
          beta.update.quad = beta.update.quad,
          initialise.delta = initialise.delta)

# Comparison w/ poisson / nbin
a <- glmmTMB(Y ~ time + cont + bin + (1|id), data, poisson) # Not particularly compelling!
# Can't fit corresponing genpois nor nbin2 model...
b <- glmmTMB(Y ~ time + cont + bin + (1|id), data, genpois, disp = ~as.factor(id) * time)
fit$logLik

# Repeat same on intercept-and-slope --------------------------------------
rm(list=ls())
source('EM.R')
source('_zzz.R')
stest <- simData_joint2_safe(n = 250, delta.int = c(2.0, 2.0), delta.time = c(0, -0.4),
                             ntms = 15, theta = c(-2, .1), fup = 3,
                             beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                             D = matrix(c(0.25, 0, 0, 0.05), 2, 2), diff.tol = 5)
data <- stest$data

# check.disps(test) # NYI
long.formula <- Y ~ time + cont + bin + (1 + time| id)

fit <- EM(long.formula, surv.formula, disp.formula, data,
          control = control,
          disp.control = disp.control,
          optim.control = optim.control,
          summax.fn = summax.fn, min.summax = min.summax,
          delta.update.quad = delta.update.quad,
          beta.update.quad = beta.update.quad,
          initialise.delta = initialise.delta)

# Comparison w/ poisson / nbin
fit$logLik
a <- glmmTMB(Y ~ time + cont + bin + (1+time|id), data, poisson) # Not particularly compelling!
b <- glmmTMB(Y ~ time + cont + bin + (1+time|id), data, genpois, disp = ~time)
a
b
fit$logLik


# Intercept-only dispersion + compares ------------------------------------
rm(list=ls())
source('EM.R')
source('_zzz.R')
stest <- simData_joint2_safe(n = 250, delta.int = c(1.0, 1.0), delta.time = c(0.0, 0.0),
                             ntms = 15, theta = c(-2, .1), fup = 3,
                             beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                             D = matrix(c(0.25, 0, 0, 0.05), 2, 2), diff.tol = 5)
data <- stest$data

# check.disps(test) # NYI
long.formula <- Y ~ time + cont + bin + (1 + time| id)
disp.formula <- ~ 1

fit <- EM(long.formula, surv.formula, disp.formula, data,
          control = control,
          disp.control = disp.control,
          optim.control = optim.control,
          summax.fn = summax.fn, min.summax = min.summax,
          delta.update.quad = delta.update.quad,
          beta.update.quad = beta.update.quad,
          initialise.delta = initialise.delta)

fit$logLik
# just poisson
system.time(a <- glmmTMB(long.formula, data, poisson))
glmmTMB:::logLik.glmmTMB(a)
# Negative binomial
aa <- glmmTMB(long.formula, data, nbinom2, dispformula = ~ as.factor(id))
# global dispersion (not subj. spec)
(b <- glmmTMB(long.formula, data, genpois, dispformula = ~ 1)) 
glmmTMB:::logLik.glmmTMB(b)
# subj-specific dispersion 
system.time(cc <- glmmTMB(long.formula, data, genpois, dispformula = ~ as.factor(id) - 1))
glmmTMB:::logLik.glmmTMB(cc)
-2 * -cc$fit$objective + 2 * (length(fit$modelInfo$inds.met) + 7) # not all inds are met.

# Quick comparison between genpois deltas and mine...
my <- fit$delta.df
tm <- data.frame(id = 1:250, delta_tmb = cc$fit$parfull[names(cc$fit$parfull) == 'betad'])
tm <- tm[tm$id %in% my$`(Intercept)`$id,]
# Look to be broadly the same (just opposite sign).
my <- merge(my$`(Intercept)`, tm, 'id')

# What about CMP?
system.time(d <- glmmTMB(long.formula, data, compois, dispformula = ~ as.factor(id) - 1))
# Doesnt everfit! sad!
