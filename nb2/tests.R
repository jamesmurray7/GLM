rm(list=ls())
source('EM.R')
d <- simData()
data <- d$data
survdata <- data[!duplicated(data[,'id']),]
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)
fit <- EM(data, ph, survdata, dispformula = ~time)
fit$coeffs; fit$SE

rm(list=ls())
source('EM.R')
d <- simData(dispFormula = ~1, alpha = 1)
data <- d$data
survdata <- data[!duplicated(data[,'id']),]
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)
fit <- EM(data, ph, survdata, dispformula = ~1)
fit$coeffs; fit$SE
