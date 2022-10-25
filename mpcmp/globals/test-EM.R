rm(list=ls())
source('EM.R')
test <- simData_joint2(n = 250, delta = c(.4,0), 
                       ntms = 10, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0), 2, 2))
data <- test$data

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1
genpois.inits <- T
delta.update.quad <- beta.update.quad <- T
# optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
individual.summax <- T
summax.fn <- NULL
min.summax <- 20
control <- list(verbose = T)

fit <- EM(long.formula, disp.formula, surv.formula, 
          data, control = list(verbose = T), 
          optim.control = list(reltol = .Machine$double.eps^(1/7)),
          genpois.inits = F,
          individual.summax = T,
          summax.fn = function(y) 2 * max(y))

setwd('../../genpois/')
source('EM.R')
fit.gp <- EM(long.formula, disp.formula, surv.formula, data = data,
             genpois.inits = T, control = list(verbose = T))

# Time-varying
rm(list=ls())
source('EM.R')
test <- simData_joint2(n = 250, delta = c(.6, -.1), 
                       ntms = 10, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0), 2, 2))
data <- test$data

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~time
genpois.inits <- T
delta.update.quad <- beta.update.quad <- T
optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
individual.summax <- T
summax.fn <- NULL
min.summax <- 20
control <- list(verbose = T)

fit <- EM(long.formula, disp.formula, surv.formula, 
          data, control = list(verbose = T), 
          optim.control = optim.control, 
          genpois.inits = genpois.inits,
          individual.summax = T,
          summax.fn = function(y) 2 * max(y))

sdat <- test$surv.data
svpart <- coxph(Surv(survtime, status) ~ bin, data = sdat)
coxph(Surv(survtime, status) ~ bin, data= sdat, init = c(fit$coeffs$zeta), 
      control=coxph.control(iter.max = 0))$loglik[2]
fit$logLik
(a <- glmmTMB(long.formula, data, family = genpois()))
adj.ll <- logLik(a) + coxph(Surv(survtime, status) ~ bin, data= test$surv.data, init = c(fit$coeffs$zeta), 
                            control=coxph.control(iter.max = 0))$loglik[2]
-2 * adj.ll + (attr(logLik(a), 'df') + 1) * 2


# Mini simulation study, to be appraised in parse.R
rm(list=ls())
source('EM.R')
data <- replicate(100, simData_joint2(n = 250, delta = c(1., -.1), 
                                      ntms = 10, theta = c(-2, .1), fup = 3,
                                      beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                      D = matrix(c(0.25, 0, 0, 0.05), 2, 2))$data
                  , simplify = F)
long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~time
genpois.inits <- T
delta.update.quad <- beta.update.quad <- T
optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
individual.summax <- T
summax.fn <- NULL
min.summax <- 20
control <- list(verbose = T)

pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  ff <- tryCatch(suppressMessages(EM(long.formula, disp.formula, surv.formula, 
                    d, control = list(verbose = F, beta.update.quad=F), 
                    optim.control = optim.control, 
                    genpois.inits = genpois.inits,
                    individual.summax = T,
                    summax.fn = function(y) 2 * max(y))),
                 error = function(e) NULL)
  fits[[i]] <- ff
  utils::setTxtProgressBar(pb, i)
}
save(fits, file = '/data/c0061461/14oct22/global_intslope2.RData')


