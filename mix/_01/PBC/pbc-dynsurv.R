setwd('~/Documents/GLMM/mix/_01')
load('pbc.RData')
pbcdata

source('EM.R')

data <- pbcdata$pbc
survdata <- pbcdata$survdata
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)

fit <- EM(data, ph, survdata, gh = 3)
fit

source('dynSurv/draw.R')
source('dynSurv/prepdata.R')
source('dynSurv/dynSurv.R')
sourceCpp('dynSurv/helper.cpp')

data[data$id == 1, ]
data %>% filter(status == 1) %>% distinct(survtime) %>% max
args(dynSurv)
ds <- dynSurv(fit, data, id = 1, u = c(2,3,4,5,6,7,8,9))
ds <- do.call(rbind, ds)
ds

plot(ds[,1] ~ c(2:9),
     ylim = c(min(pmin(ds[, 2], ds[, 3])), max(pmax(ds[,2], ds[,3]))),
     type = 'o', col = 'red',
     xlab = 'u', ylab = 'S')
lines(ds[,2] ~ c(2:9), lty = 5)
lines(ds[,3] ~ c(2:9), lty = 5)

# bit of a longer profile

data[data$status == 1, ]
data %>% filter(status == 1) %>% distinct(survtime) %>% max
args(dynSurv)
ds <- dynSurv(fit, data, id = 8, u = c(7,8,9,10,11))
ds <- do.call(rbind, ds)

plot(ds[,1] ~ c(7:11),
     ylim = c(min(pmin(ds[, 2], ds[, 3])), max(pmax(ds[,2], ds[,3]))),
     type = 'o', col = 'red',
     xlab = 'u', ylab = 'S')
lines(ds[,2] ~ c(7:11), lty = 5)
lines(ds[,3] ~ c(7:11), lty = 5)

# Someone who survives

data[data$status == 0, ]
ds <- dynSurv(fit, data, id = 2, u = c(9, 10, 11, 12, 13))
ds <- do.call(rbind, ds)
plot(ds[,1] ~ c(9:13),             # So model would expect this person to fail?
     ylim = c(min(pmin(ds[, 2], ds[, 3])), max(pmax(ds[,2], ds[,3]))),
     type = 'o', col = 'red',
     xlab = 'u', ylab = 'S')
lines(ds[,2] ~ c(9:13), lty = 5)
lines(ds[,3] ~ c(9:13), lty = 5)