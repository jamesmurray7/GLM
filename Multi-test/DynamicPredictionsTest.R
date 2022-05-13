rm(list=ls())
source('EM.R')
source('DynamicPredictions.R')
load('../PBC-case-study/PBC.RData')
pbc$serBilir <- log(pbc$serBilir)

#' Univariate
long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
surv.formula <- Surv(survtime, status) ~ drug
family <- list(gaussian)
my.fit <- EM(long.formulas, surv.formula, pbc, family, control = list(hessian='manual'))

library(joineRML)
jML.fit <- mjoint(
  formLongFixed = list('1' = serBilir ~ drug * time),
  formLongRandom = list('1' = ~ time | id),
  formSurv = Surv(survtime, status) ~ drug,
  data = pbc, timeVar = 'time', control = list(
    type = 'sobol', convCrit = 'rel', tol2 = 1e-2, tol.em = 5e-3
  ), verbose = T
)

test <- pbc[pbc$id == 2 & pbc$time <= 9,]
my.fit$hazard[,1]->ft
u <- ft[ft <= 11 & ft > 9]
jmlsurv <- joineRML::dynSurv(jML.fit, test, u = u, type = 'simulated', M = 200)
mine_normal <- dynSurv2(data = pbc, id = 2, fit = my.fit, u = u, nsim = 200, b.density = 'normal')
mine_t <- dynSurv2(data = pbc, id = 2, fit = my.fit, u = u, nsim = 200, b.density = 't', scale = 2, df = 4)

ROCt <- ROC(my.fit, pbc, 9, 2, control = list(b.density = 't', scale = 2, df = 4, nsim = 25))
ROCn <- ROC(my.fit, pbc, 9, 2, control = list(b.density = 'normal', nsim = 25)) # t seems to perform a little better here!
