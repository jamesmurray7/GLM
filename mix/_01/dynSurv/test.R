rm(list=ls())
setwd('~/Documents/GLMM/mix/_01')
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

ds <- dynSurv(fit, data, 1, u = c(4, 5, 6, 7, 8, 9))
ds <- do.call(rbind, ds)
max(pmax(ds[,2],ds[,3]))
min(pmin(ds[,2],ds[,3]))

plot(ds[,1]~c(4:9), type = 'o', ylim = c(min(pmin(ds[,2],ds[,3])),
                                         max(pmax(ds[,2],ds[,3]))),
     col = 'red', ylab = 'S', xlab = 'u')
lines(ds[,3]~c(4:9), lty = 3)
lines(ds[,2]~c(4:9), lty = 3)
abline(v = unique(data[data$id == 1, 'survtime']), lty = 5, col = 'red')
