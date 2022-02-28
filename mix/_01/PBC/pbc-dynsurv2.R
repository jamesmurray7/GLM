rm(list=ls())

load('./PBC/myfit-hepatomegaly.RData')
load('./PBC/pbc-hepatomegaly.RData')

source('dynSurv/dynSurv.R')
source('dynSurv/draw.R')
source('dynSurv/prepdata.R')
sourceCpp('dynSurv/helper.cpp')

data <- pbcdata$pbc

dynSurv(fit, data, 1, 3)


prepdata
dynSurv

load('~/Downloads/jmfit-hepa.RData')
jmROC <- JMbayes2:::tvROC.jm(jmfits[[1]], data, Tstart = 8, Dt = 2)
myROC <- ROC(fit, data, Tstart = 8, Delta = 2, what = 'last', nsim = 1)
myROC_simlower <- ROC(fit, data, Tstart = 8, Delta = 2, what = 'lowest', nsim = 10)

myROC
myROC_simlower

AUC(myROC, reduced = F)

par(mfrow=c(1,3))
plotROC(myROC)
plotROC(myROC_simlower)
plot(jmROC, main = 'JMbayes2')
legend('bottomright', paste('AUC', round(tvAUC(jmROC)$auc, 2)), bty = 'n')
dev.off()
