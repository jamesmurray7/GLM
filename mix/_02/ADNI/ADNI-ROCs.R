#' #########
#' ADNI-ROCs
#' Comparing ROC curve EM/JMbayes2. WD: mix/_01
#' #########

rm(list=ls())
source('dynSurv/dynSurv.R')
source('dynSurv/draw.R')
source('dynSurv/prepdata.R')
source('EM.R')
sourceCpp('dynSurv/helper.cpp')

# >= 80% ------------------------------------------------------------------
load('ADNI/fit-ADNI80.RData')
load('ADNI/jm-ADNI80.RData')
load('ADNI/ADNI80.RData')
data <- ADNI$data
ROC80 <- ROC(fit80, data, Tstart = 6, Delta = 3, nsim = 5)
jmROC80 <- JMbayes2::tvROC(jm80, data, Tstart = 6, Dt = 3)
save(ROC80, file = 'ADNI/ROC-ADNI80.RData')
save(jmROC80, file = 'ADNI/jm-ROC-ADNI80.RData')

# >= 100% -----------------------------------------------------------------
load('ADNI/fit-ADNI100.RData')
load('ADNI/jm-ADNI100.RData')
load('ADNI/ADNI100.RData')
data <- ADNI$data
ROC100.2 <- ROC(fit100, data, Tstart = 6, Delta = 3, nsim = 5, smoothed = T)
jmROC100 <- JMbayes2::tvROC(jm100, data, Tstart = 6, Dt = 3)
save(ROC100, file = 'ADNI/ROC-ADNI100.RData')
save(jmROC100, file = 'ADNI/jm-ROC-ADNI100.RData')

par(mfrow = c(2, 1))
plotROC(ROC100, legend = T)
plotROC(ROC100.2, legend = T)
plot(jmROC100); JMbayes2::tvAUC(jmROC100)


# Rolling windows for Delta = 2 -------------------------------------------
# >= 80%
data <- data[!data$id %in% c(286, 333), ]
mine <- JMB <- list(); p <- 1
for(t in c(2, 4, 5, 7)){
  # Stick with one draw?
  myROCt <- ROC(fit80, data, Tstart = t, Delta = 2, nsim = 1, what = 'last')
  jmROCt <- tvROC(jm80, data, Tstart = t, Dt = 2)
  JMB[[p]] <- jmROCt; mine[[p]] <- myROCt
  p <- p + 1
}

pdf('~/Downloads/cc.pdf', paper ='a4r', width = 144, height = 60)
par(mfrow = c(2, 4))
lapply(mine, plotROC, legend = T)
lapply(JMB, function(x){
  plot(x)
  legend('bottomright', paste('AUC', round(tvAUC(x)$auc, 3)))
})
dev.off()





