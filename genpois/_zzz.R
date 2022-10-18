long.formula = Y ~ time + cont + bin + (1|id)
surv.formula = Surv(survtime, status) ~ bin
disp.formula = ~ 1
control = list(verbose = T, phi.update.quad = T, beta.update.quad = F, include.all = T)
genpois.inits = T
