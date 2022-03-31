rm(list=ls())
source('EM.R')
d <- simData_joint()
data <- d$data
survdata <- d$surv.data
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)

fitNQ <- EM(data, ph, survdata, gh = 3, quad = F, nb = F, verbose = T)
fitQ <- EM(data, ph, survdata, gh = 9, quad = T, nb = F, verbose = T)
  