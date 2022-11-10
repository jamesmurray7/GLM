long.formula <- Y ~ time + cont + bin + (1|id)
surv.formula <- Surv(survtime, status) ~ bin
control <- list(verbose = T)