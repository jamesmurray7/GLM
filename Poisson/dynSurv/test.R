setwd('~/Documents/GLMM/Poisson/')
rm(list=ls())
source('EM.R')

dd <- simData_joint()
data <- dd$data
survdata <- dd$surv.data

ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)
fit <- EM(data, ph, survdata, gh = 3)

data %>% filter(status == 1) %>% distinct(id) %>% pull -> fail.id
data %>% filter(id %in% fail.id)

source('dynSurv/dynSurv.R')
source('dynSurv/draw.R')
source('dynSurv/prepdata.R')
sourceCpp('dynSurv/helper.cpp')

ds <- dynSurv(fit, data, 2, u = c(2,3,4,5,6,7,8,9))
ds <- do.call(rbind, ds)
ds
max(pmax(ds[,2],ds[,3]))
min(pmin(ds[,2],ds[,3]))

plot(ds[,1]~c(2:9), type = 'o', ylim = c(min(pmin(ds[,2],ds[,3])),
                                         max(pmax(ds[,2],ds[,3]))),
     col = 'red', ylab = 'S', xlab = 'u')
lines(ds[,3]~c(2:9), lty = 3)
lines(ds[,2]~c(2:9), lty = 3)
