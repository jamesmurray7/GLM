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


#' ROC for interval (3, 5], (5, 7] and (9, 11]
myROC35  <- ROC(my.fit, pbc, Tstart = 3, delta = 2, 
                control = list(
                  nsim = 10
                ))
myROC57 <- ROC(my.fit, pbc, Tstart = 5, delta = 2, nsim = 10)
myROC911 <- ROC(my.fit, pbc, Tstart = 9, delta = 2,
                control = list(
                  nsim = 10, b.dens = 't', scale = 2, df = 4
                ))
myROC911.2 <- ROC(my.fit, pbc, Tstart = 9, delta = 2,
                control = list(
                  nsim = 10, b.dens = 'normal', scale = 2, df = 4
                ))

par(mfrow=c(3,1))
plotROC(myROC35, T)
plotROC(myROC57, T)
plotROC(myROC911, T)
dev.off()

#' Univariate with quadratic fixed and random structure.
long.formulas.quad <- list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id))
my.fit.quad <- EM(long.formulas.quad, surv.formula, pbc, family, control = list(hessian='auto'))

#' ROC for interval (3, 5], (5, 7] and (9, 11] for quadratic specification.
myROC35.quad <- ROC(my.fit.quad, pbc, Tstart = 3, delta = 2, 
                    control = list(
                      nsim = 25, b.dens = 't', scale = 2, df = 4
                    ))



myROC57.quad  <- ROC(my.fit.quad, pbc, Tstart = 5, delta = 2, control = list(nsim = 25, b.dens = 't', scale = 2, df = 4))
myROC911.quad <- ROC(my.fit.quad, pbc, Tstart = 9, delta = 2, control = list(nsim = 25, b.dens = 't', scale = 2, df = 4))
myROC911.quad2 <- ROC(my.fit.quad, pbc, Tstart = 9, delta = 2, control = list(nsim = 25, b.dens = 'normal'))

par(mfrow = c(3, 2))
plotROC(myROC35, T); plotROC(myROC35.quad, T)
plotROC(myROC57, T); plotROC(myROC57.quad, T)
plotROC(myROC911, T); plotROC(myROC911.quad, T)
dev.off()

#' Univariate with polynomial (natural cubic splines) time specification in fixed and random effects.
#' long.formulas2 <- list(serBilir ~ drug * splines::ns(time, df = 3) + (1+splines::ns(time, df = 3)|id))
#' 
#' my.fit.ns <- EM(long.formulas2, surv.formula, pbc, family, control = list(hessian='auto'))
#' 
#' #' ROC for interval (3, 5], (5, 7] and (9, 11] with this cubic spline specification.
#' myROC35.ns  <- ROC(my.fit.ns, pbc, Tstart = 3, delta = 2, nsim = 10)
#' myROC57.ns  <- ROC(my.fit.ns, pbc, Tstart = 5, delta = 2, nsim = 10)
#' myROC911.ns <- ROC(my.fit.ns, pbc, Tstart = 9, delta = 2, nsim = 10)
#' 
#' par(mfrow = c(3,3))
#' plotROC(myROC35, T);  plotROC(myROC35.quad, T);   plotROC(myROC35.ns, T)
#' plotROC(myROC57, T);  plotROC(myROC57.quad, T);   plotROC(myROC57.ns, T)
#' plotROC(myROC911, T); plotROC(myROC911.quad, T);  plotROC(myROC911.ns, T)
#' dev.off()
  

# Bivariate any good? ----------------------------------------------------
pbc$SGOT <- log(pbc$SGOT)
long.formulas <- list(
  serBilir ~ drug*time + (1 + time|id),
  albumin ~ drug * time + (1 + time|id),
  SGOT ~ drug * time + (1 + time|id)
)
families <- list(gaussian, gaussian, gaussian)
fit <- EM(long.formulas, surv.formula, pbc, families, control = list(hessian = 'auto'))  

myROC <- ROC(fit, pbc, Tstart = 8, delta = 4, control = list(nsim = 1, b.dens = 'normal'))
myROC2 <- ROC(fit, pbc, Tstart = 8, delta = 4, control = list(nsim = 1, b.dens = 't', scale = 2, df = 4))

normals <- ts <- vector('list', 4)
for(i in 1:4){
  normals[[i]] <- ROC(fit, pbc, Tstart = 8, delta = i, control = list(nsim = 10, b.density = 'normal'))
  ts[[i]] <- ROC(fit, pbc, Tstart = 8, delta = i, control = list(nsim = 10, b.density = 't'))
}

par(mfrow=c(2,2))
lapply(normals, plotROC, T)
lapply(ts, plotROC, T)
