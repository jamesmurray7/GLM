#' ########
#' Dynpreds for the model fit fit
#' with Gaussian model for platelet count in place of negbin.
#' ########
rm(list=ls())
setwd('~/Documents/PhD/GLM/mix/_02/')
source("EM.R")

load('../_01/PBC/pbc-hepatomegaly.RData')
data <- pbcdata$pbc
survdata <- pbcdata$survdata
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)

nlmefit <- nlme::lme(Y.3 ~ time + cont + bin, random = ~time|id , data = data)
plot(nlmefit)

fitG <- EM(data, ph, survdata, gh = 3, prepopulateD = T)

source('dynSurv/draw.R')
source('dynSurv/dynSurv.R')
source('dynSurv/prepdata.R')
sourceCpp('dynSurv/helper.cpp')

ROCs <- list(); p <- 1
for(r in c(3,5,7,9)){
  ROCS[[p]] <- ROC(fitG, data, Tstart = r, Delta = 2, what = 'last', nsim = 10)
  p <- p + 1
}
