rm(list=ls())

load('./PBC/myfit-hepatomegaly.RData')
load('./PBC/pbc-hepatomegaly.RData')

source('dynSurv/dynSurv.R')
source('dynSurv/draw.R')
source('dynSurv/prepdata.R')
sourceCpp('dynSurv/helper.cpp')

data <- pbcdata$pbc

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

# Looks to be in order...

# Can we do rolling 2-Dt windows?

JMB <- mine <- list()
p <- 1
for(t in c(3, 5, 7, 9)){
  # Stick with one draw?
  myROCt <- ROC(fit, data, Tstart = t, Delta = 2, nsim = 1)
  jmROCt <- tvROC(jmfits[[1]], data, Tstart = t, Dt = 2)
  JMB[[p]] <- jmROCt; mine[[p]] <- myROCt
  p <- p + 1
}

save(JMB, file = '~/Downloads/JMB.RData')
save(mine, file = '~/Downloads/mine.RData')

lapply(mine, AUC)
lapply(JMB, function(x) tvAUC(x)$auc)

pdf('~/Downloads/aa.pdf', paper ='a4r', width = 144, height = 60)
par(mfrow = c(2,4))
lapply(mine, plotROC, legend = T)
lapply(JMB, function(x){
  plot(x)
  legend('bottomright', paste('AUC', round(tvAUC(x)$auc, 3)))
  })
dev.off()


mine2 <- list(); p <- 1
for(t in c(3, 5, 7, 9)){
  # Stick with one draw?
  myROCt <- ROC(fit, data, Tstart = t, Delta = 2, nsim = 50, what = 'lowest')
  mine2[[p]] <- myROCt
  p <- p + 1
}

pdf('~/Downloads/bb.pdf', paper ='a4r', width = 144, height = 60)
par(mfrow = c(2,4))
lapply(mine,  plotROC, legend = T)
lapply(mine2, plotROC, legend = T) # Sig. better at increased computational cost.
dev.off()
