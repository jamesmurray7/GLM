rm(list=ls())
source('EM.R')
source('DynamicPredictions.R')
load('../PBC-case-study/PBC.RData')
pbc$serBilir <- log(pbc$serBilir)

#' Univariate
long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
surv.formula <- Surv(survtime, status) ~ drug
family <- list(gaussian)
my.fit <- EM(long.formulas, surv.formula, pbc, family, control = list(hessian='auto'))


#' ROC for interval (3, 5], and (9, 11]
myROC35  <- ROC(my.fit, pbc, Tstart = 3, delta = 2, nsim = 10)
myROC911 <- ROC(my.fit, pbc, Tstart = 9, delta = 2, nsim = 10)

plotROC(myROC35, T)
plotROC(myROC911, T)

#' Univariate but add in some splines
long.formulas2 <- list(serBilir ~ drug * splines::ns(time, df = 3) + (1+splines::ns(time, df = 3)|id))

my.fit2 <- EM(long.formulas2, surv.formula, pbc, family, control = list(hessian='auto'))

#' ROC for interval (3, 5], and (9, 11], any better using splines?
myROC35.spline  <- ROC(my.fit2, pbc, Tstart = 3, delta = 2, nsim = 10)
myROC911.spline <- ROC(my.fit2, pbc, Tstart = 9, delta = 2, nsim = 10)

par(mfrow = c(2,2))
plotROC(myROC35, T); plotROC(myROC35.spline, T)   # No obvious preferential approach.
plotROC(myROC911, T);plotROC(myROC911.spline, T)
dev.off()
  



