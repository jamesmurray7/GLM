#' #############
#' Tests of EM fits.
#' Y_1: Gaussian; Y_2: binary; Y_3: Poisson 
#'  args(simData_joint) = 
#'  beta =      [,1] [,2] [,3] [,4]
#'        [1,]  0.2 -0.1  0.1 -0.2
#'        [2,]  1.0 -1.0  1.0 -1.0
#'        [3,]  3.0 -0.1  0.1 -0.2
#'  diag(D) = 0.16 0.09 0.25 0.25 0.16
#'  gamma = 0.5, 0.3, -0.2
#'  var.e = 0.16, 
#' #########

rm(list=ls())
source('EM.R')

rm(list=ls())
N <- 100
source('EM.R')
# Simulate some different datasets
# data1 <- replicate(N, simData_joint(n = 250, theta = c(-4.4, 0.1)), simplify = F) #20%
# data2 <- replicate(N, simData_joint(n = 500, theta = c(-4.4, 0.1)), simplify = F)
data1 <- replicate(N, simData_joint(n = 500, theta = c(-2.85, 0.2)), simplify = F) # 40%
data2 <- replicate(N, simData_joint(n = 500, theta = c(-3.75, 0.2)), simplify = F) # 20%

# check average failure rate...
quantile(do.call(c, lapply(data1, function(x) sum(x$surv$status == 1)/nrow(x$surv)))) #20%
quantile(do.call(c, lapply(data2, function(x) sum(x$surv$status == 1)/nrow(x$surv)))) #40%

fit1 <- fit2 <- vector('list', N)
pb <- utils::txtProgressBar(max = N, style = 3)
for(i in 1:N){
  ph1 <- coxph(Surv(survtime, status) ~ 1, data1[[i]]$surv.data)
  ph2 <- coxph(Surv(survtime, status) ~ 1, data2[[i]]$surv.data)

  # fit
  fit1[[i]] <- tryCatch(suppressMessages(EM(data1[[i]]$data, ph1, data1[[i]]$surv.data, verbose = F, gh = 3, quad = T)),
                         error = function(e) NULL)
  fit2[[i]] <- tryCatch(suppressMessages(EM(data2[[i]]$data, ph2, data2[[i]]$surv.data, verbose = F, gh = 3, quad = T)),
                         error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}

fits <- list(`40%` = fit1, `20%` = fit2)
save(fits, file = '~/Downloads/fits-rustand-failure2040.RData')
save(fitsQ, file = '~/Downloads/Quad-fits-rustand-ntms.RData')
save(fitsNQ, file = '~/Downloads/Noquad-fits-rustand-ntms.RData')
