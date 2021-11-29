rm(list=ls())
setwd('~/Documents/GLMM/zip/zhu')
source('EM.R')
data6 <- replicate(100, 
                  simData_zip_joint(n = 250), simplify = F)
data10 <- replicate(100, 
                  simData_zip_joint(n = 250, ntms = 10), simplify = F)
data14 <- replicate(100, 
                  simData_zip_joint(n = 250, ntms = 14), simplify = F)

pb <- utils::txtProgressBar(max=100, style = 3)
fits6 <- fits10 <- fits14 <- list()
for(i in 1:100){
  ph6 <- coxph(Surv(survtime, status) ~ bin, data = data6[[i]]$surv.data)
  ph10 <- coxph(Surv(survtime, status) ~ bin, data = data10[[i]]$surv.data)
  ph14 <- coxph(Surv(survtime, status) ~ bin, data = data14[[i]]$surv.data)
  fits6[[i]] <- tryCatch(suppressMessages(EM(data6[[i]]$data, ph6, data6[[i]]$surv.data)),
                    error = function(e) NULL)
  fits10[[i]] <- tryCatch(suppressMessages(EM(data10[[i]]$data, ph10, data10[[i]]$surv.data)),
                    error = function(e) NULL)
  fits14[[i]] <- tryCatch(suppressMessages(EM(data14[[i]]$data, ph14, data14[[i]]$surv.data)),
                    error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}

fits <- list(fits6=fits6,fits10=fits10,fits14=fits14)
save(fits, file = '~/Downloads/fits291121.RData')
