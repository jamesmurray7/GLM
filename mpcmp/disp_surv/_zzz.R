control <- list(verbose=T, debug = T)
disp.control <- list(delta.method = 'optim', min.profile.length = 3,
                     truncated = T, max.val = 2.25)
optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
summax=3

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~time
update.deltas <- F

summax.fn <- NULL
min.summax <- 20
delta.update.quad <- T
beta.update.quad <- F
initialise.delta <- F
