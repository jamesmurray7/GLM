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
    type = 'sobol', convCrit = 'rel', tol2 = 1e-1, tol.em = 5e-3
  ), verbose = F
)

test <- pbc[pbc$id == 81 & pbc$time <= 6.9,]
my.fit$hazard[,1]->ft
u <- ft[ft <= 15 & ft > 6.9]
jmlsurv <- joineRML::dynSurv(jML.fit, test, u = u, type = 'simulated', M = 200)
mine_normal <- dynSurv2(data = pbc, id = 81, fit = my.fit, u = u, nsim = 200, b.density = 'normal')
mine_t <- dynSurv2(data = pbc, id = 81, fit = my.fit, u = u, nsim = 200, b.density = 't', scale = 2, df = 4)

jmlsurv;mine_normal$pi;mine_t$pi

ROCt <- ROC(my.fit, pbc, 9, 2, control = list(b.density = 't', scale = 2, df = 4, nsim = 50))
ROCn <- ROC(my.fit, pbc, 9, 2, control = list(b.density = 'normal', nsim = 50)) # t seems to perform a little better here!

plotROC(ROCt, T)
plotROC(ROCn, T)

Ts <- c(3,5,7,9)
ROCs <- setNames(vector('list', 4), Ts)

for(t in seq_along(Ts)){
  ROCs[[t]] <- ROC(my.fit, pbc, Ts[t], delta = 2, 
                   control = list(b.density = 't', scale = 2, df = 4, nsim = 100))
}

par(mfrow=c(2,2))
lapply(ROCs, plotROC, T)
dev.off()


# Bivariate ---------------------------------------------------------------
rm(list=ls())
source('EM.R')
source('DynamicPredictions.R')
load('../PBC-case-study/PBC.RData')
pbc <- pbc[!is.na(pbc$spiders) & !is.na(pbc$serBilir),]
pbc$serBilir <- log(pbc$serBilir)
long.formulas <- list(serBilir ~ drug * time + (1 + time|id),
                      spiders ~ drug * time + (1 + time|id))
surv.formula <- Surv(survtime, status) ~ drug
family <- list(gaussian, binomial)
my.fit <- EM(long.formulas, surv.formula, pbc, family, control = list(hessian='manual',
                                                                      optimiser='optim'))
my.summary(my.fit)



ROCt <- ROC(my.fit, pbc, 9, 2, control = list(b.density = 't', scale = 2, df = 4, nsim = 25))
