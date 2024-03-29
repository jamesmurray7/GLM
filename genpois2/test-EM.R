# genpois
rm(list=ls())
source('EM.R')
test <- simData_joint(delta=-0.5)
data <- test$data
source('_zzz.R')

fit <- EM(long.formula, disp.formula, surv.formula, 
          data, control = control,
          genpois.inits = T)

sdat <- test$surv.data
svpart <- coxph(Surv(survtime, status) ~ bin, data = sdat)
coxph(Surv(survtime, status) ~ bin, data= sdat, init = c(fit$coeffs$zeta), 
      control=coxph.control(iter.max = 0))$loglik[2]
fit$logLik
(a <- glmmTMB(long.formula, data, family = genpois()))
adj.ll <- logLik(a) + coxph(Surv(survtime, status) ~ bin, data= sdat, init = c(fit$coeffs$zeta), 
                  control=coxph.control(iter.max = 0))$loglik[2]
-2 * c(adj.ll) + (attr(logLik(a), 'df') + 1) * 2


# Intslope ----------------------------------------------------------------
rm(list=ls())
source('EM.R')
test <- simData_joint(D = matrix(c(.25, 0,0,.05), 2, 2), diff.tol = 99,
                      beta = c(4, -.1, .1, -0.5), phi = .5)
data <- test$data
source('_zzz.R')
long.formula <- Y~time+cont+bin+(1+time|id)

fit <- EM(long.formula, disp.formula, surv.formula, 
          data, control = control,
          genpois.inits = T)

sdat <- test$surv.data
svpart <- coxph(Surv(survtime, status) ~ bin, data = sdat)
coxph(Surv(survtime, status) ~ bin, data= sdat, init = c(fit$coeffs$zeta), 
      control=coxph.control(iter.max = 0))$loglik[2]
fit$logLik
(a <- glmmTMB(long.formula, data, family = genpois()))
adj.ll <- logLik(a) + coxph(Surv(survtime, status) ~ bin, data= sdat, init = c(fit$coeffs$zeta), 
                            control=coxph.control(iter.max = 0))$loglik[2]
-2 * c(adj.ll) + (attr(logLik(a), 'df') + 1) * 2


# Mini simulation study, to be appraised in parse.R
rm(list=ls())
setwd('../genpois2')
source('EM.R')
data <- replicate(100, simData_joint()$data, simplify = F) # default: phi=-0.5
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1
genpois.inits <- T
control <- list(verbose = F, delta.update.quad = T, include.all = T)


pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  ff <- tryCatch(suppressMessages(EM(long.formula, disp.formula, surv.formula, 
                                     d, control = control,
                                     genpois.inits = T)),
                 error = function(e) NULL)
  fits[[i]] <- ff
  utils::setTxtProgressBar(pb, i)
}
save(fits, file = './gp1set.RData')


fixef(a)$disp

setwd('../mpcmp/globals/')
source('EM.R')
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1
genpois.inits <- T
delta.update.quad <- beta.update.quad <- T
optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
individual.summax <- T
summax.fn <- NULL
min.summax <- 20

pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  ff <- tryCatch(suppressMessages(EM(long.formula, disp.formula, surv.formula,
                                     d, genpois.inits = T, optim.control = optim.control)),
                 error = function(e) NULL)
  fits[[i]] <- ff
  utils::setTxtProgressBar(pb, i)
}
save(fits, file = './mpcmpset.RData')

