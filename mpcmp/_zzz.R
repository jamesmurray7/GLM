control <- list(verbose=T)
disp.control <- list()
grid.summax <- 'a'
summax <- 100
delta.init <- NULL

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1