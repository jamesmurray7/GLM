#' #############
#' Tests of EM fits.
#  args(simData_joint) = 
#  beta = c(1, 0.1, 0.33, -0.5), theta = 1.5
#  D = matrix(c(0.5, 0, 0, 0.1), 2, 2), gamma = 0.5, surv.eta = c(0.05, -0.3)
#' #########

rm(list=ls())
source('EM.R')

# One fit...
dd <- simData_joint(thetaDisp = 1)
ph <- coxph(Surv(survtime, status) ~ cont + bin, dd$surv.data)
fit2 <- EM(dd$data, ph, dd$surv.data, verbose = T)

# Many fits...
rm(list=ls())
source('EM.R')
# Simulate some different datasets
data1 <- replicate(100, simData_joint(ntms = 15, thetaDisp = 0.25), simplify = F)
data2 <- replicate(100, simData_joint(ntms = 15, thetaDisp = 0.50), simplify = F)
data3 <- replicate(100, simData_joint(ntms = 15, thetaDisp = 0.75), simplify = F)
data4 <- replicate(100, simData_joint(ntms = 15, thetaDisp = 1.00), simplify = F)
data5 <- replicate(100, simData_joint(ntms = 15, thetaDisp = 1.50), simplify = F)
data6 <- replicate(100, simData_joint(ntms = 15, thetaDisp = 2.00), simplify = F)

fits1 <- fits2 <- fits3 <- fits4 <- fits5 <- fits6 <- list()
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  ph1 <- coxph(Surv(survtime, status) ~ cont + bin, data1[[i]]$surv.data)
  ph2 <- coxph(Surv(survtime, status) ~ cont + bin, data2[[i]]$surv.data)
  ph3 <- coxph(Surv(survtime, status) ~ cont + bin, data3[[i]]$surv.data)
  ph4 <- coxph(Surv(survtime, status) ~ cont + bin, data4[[i]]$surv.data)
  ph5 <- coxph(Surv(survtime, status) ~ cont + bin, data5[[i]]$surv.data)
  ph6 <- coxph(Surv(survtime, status) ~ cont + bin, data6[[i]]$surv.data)
  fits1[[i]] <- tryCatch(suppressMessages(EM(data1[[i]]$data, ph1, data1[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  fits2[[i]] <- tryCatch(suppressMessages(EM(data2[[i]]$data, ph2, data2[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  fits3[[i]] <- tryCatch(suppressMessages(EM(data3[[i]]$data, ph3, data3[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  fits4[[i]] <- tryCatch(suppressMessages(EM(data4[[i]]$data, ph4, data4[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  fits5[[i]] <- tryCatch(suppressMessages(EM(data5[[i]]$data, ph5, data5[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  fits6[[i]] <- tryCatch(suppressMessages(EM(data6[[i]]$data, ph6, data6[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}
fits <- list(fits1, fits2, fits3, fits4, fits5, fits6)
save(fits, file = '~/Downloads/nbfits.RData')