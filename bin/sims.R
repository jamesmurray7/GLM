#' #############
#' Tests of EM fits.
# args(simData_joint) = 
#   beta = c(1, 0.1, 0.33, -0.5), 
#   D = matrix(c(0.5, 0, 0, 0.1), 2, 2), gamma = 0.5, surv.eta = c(0.05, -0.3)
#' #########

rm(list=ls())
source('EM.R')

# One fit...
dd <- simData_joint()
ph <- coxph(Surv(survtime, status) ~ cont + bin, dd$surv.data)
fit <- EM(dd$data, ph, dd$surv.data, verbose = T)
fit$coeffs

# Many fits...
rm(list=ls())
source('EM.R')
# Simulate some different datasets
data1 <- replicate(100, simData_joint(), simplify = F)
data2 <- replicate(100, simData_joint(ntms = 15), simplify = F)
data3 <- replicate(100, simData_joint(gamma = -1, theta = c(-3, 0.2)), simplify = F)

fits1 <- fits2 <- fits3 <- list()
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  ph1 <- coxph(Surv(survtime, status) ~ cont + bin, data1[[i]]$surv.data)
  ph2 <- coxph(Surv(survtime, status) ~ cont + bin, data2[[i]]$surv.data)
  ph3 <- coxph(Surv(survtime, status) ~ cont + bin, data3[[i]]$surv.data)
  fits1[[i]] <- tryCatch(suppressMessages(EM(data1[[i]]$data, ph1, data1[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  fits2[[i]] <- tryCatch(suppressMessages(EM(data2[[i]]$data, ph2, data2[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  fits3[[i]] <- tryCatch(suppressMessages(EM(data3[[i]]$data, ph3, data3[[i]]$surv.data, verbose = F)),
                         error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}
fits <- list(fits1, fits2, fits3)
save(fits, file = '~/Downloads/fits011221.RData')