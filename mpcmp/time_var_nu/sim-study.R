rm(list=ls())
source('EM.R')
tests <- replicate(100, simData_joint2(n = 250, delta = c(1, -1.), 
                        ntms = 15, theta = c(-2, .1), fup = 3,
                        beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                        D = matrix(c(0.25, 0, 0, 0), 2, 2)), simplify = F)

save(tests, file = '/data/c0061461/sept-29-sims/data.RData')
fits.inits.quad <- fits.none.none <- fits.none.delta <- vector('list', 100)

# Static stuff...
disp.control <- list(delta.method = 'optim', min.profile.length = 3,
                     truncated = T, max.val = 2.5)
control <- list()
optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
update.deltas <- F

summax.fn <- function(y) max(y) + 10
min.summax <- 20

pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d <- tests[[i]]$data
  # Inits + quadrature on all...
  delta.update.quad <- T
  beta.update.quad <- T
  initialise.delta <- T
  
  fit <- tryCatch(suppressMessages(
    EM(long.formula, surv.formula, d,
       control = control,
       disp.control = disp.control,
       optim.control = optim.control,
       summax.fn = summax.fn, min.summax = min.summax,
       delta.update.quad = delta.update.quad,
       beta.update.quad = beta.update.quad,
       initialise.delta = initialise.delta)
  ), error = function(e) NULL)
  fits.inits.quad[[i]] <- fit
  
  # No inits and no quadrature
  delta.update.quad <- F
  beta.update.quad <- F
  initialise.delta <- F
  
  fit <- tryCatch(suppressMessages(
    EM(long.formula, surv.formula, d,
       control = control,
       disp.control = disp.control,
       optim.control = optim.control,
       summax.fn = summax.fn, min.summax = min.summax,
       delta.update.quad = delta.update.quad,
       beta.update.quad = beta.update.quad,
       initialise.delta = initialise.delta)
  ), error = function(e) NULL)
  fits.none.none[[i]] <- fit
  
  # No inits and delta quadrature
  delta.update.quad <- T
  beta.update.quad <- F
  initialise.delta <- F
  
  fit <- tryCatch(suppressMessages(
    EM(long.formula, surv.formula, d,
       control = control,
       disp.control = disp.control,
       optim.control = optim.control,
       summax.fn = summax.fn, min.summax = min.summax,
       delta.update.quad = delta.update.quad,
       beta.update.quad = beta.update.quad,
       initialise.delta = initialise.delta)
  ), error = function(e) NULL)
  fits.none.delta[[i]] <- fit
  
  utils::setTxtProgressBar(pb, i)
}

save(fits.inits.quad, file = '/data/c0061461/sept-29-sims/fit1.RData')
save(fits.none.none, file = '/data/c0061461/sept-29-sims/fit2.RData')
save(fits.none.delta, file = '/data/c0061461/sept-29-sims/fit3.RData')


# The same with Poisson (from multi-test!) --------------------------------
setwd('../../Multi-test/')
source('EM.R')
fits.poisson <- vector('list', 100)
l.long.formula <- list(long.formula)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d <- tests[[i]]$data
  fit <- tryCatch(suppressMessages(EM(l.long.formula, surv.formula, d, list('poisson'))),
                  error = function(e) NULL)
  fits.poisson[[i]] <- fit
  utils::setTxtProgressBar(pb, i)
}