#' #############
#' Tests of EM fits.
#'  args(simData_joint) = 
#'  beta =      [,1] [,2] [,3]  [,4]
#'         [1,]    0 1.00 0.50 -0.25   // Gauss
#'         [2,]    1 0.10 0.33 -0.50   // bin
#'         [3,]    0 0.05 0.01  0.30   // count
#'  diag(D) = 0.25 0.06 0.50 0.04 0.25 0.05
#'  gamma = 0.5, -0.25, 0.40,
#'  eta = 0.05, -0.3,
#'  var.e = 0.25, 
#' #########

rm(list=ls())
source('EM.R')

rm(list=ls())
source('EM.R')
# Simulate some different datasets
data1 <- replicate(100, simData_joint(n = 250, theta = c(-4.4, 0.1)), simplify = F)
data2 <- replicate(100, simData_joint(n = 250, theta = c(-3.8, 0.2)), simplify = F)

# check average failure rate...
quantile(do.call(c, lapply(data1, function(x) sum(x$surv$status == 1)/nrow(x$surv)))) #20%
quantile(do.call(c, lapply(data2, function(x) sum(x$surv$status == 1)/nrow(x$surv)))) #40%

fits1Q <- fits2Q <- fits1NQ <- fits2NQ <- list()
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  ph1 <- coxph(Surv(survtime, status) ~ 1, data1[[i]]$surv.data)
  ph2 <- coxph(Surv(survtime, status) ~ 1, data2[[i]]$surv.data)

  fits1Q[[i]] <- tryCatch(suppressMessages(EM(data1[[i]]$data, ph1, data1[[i]]$surv.data, verbose = F, gh = 3, quad = T)),
                         error = function(e) NULL)
  fits2Q[[i]] <- tryCatch(suppressMessages(EM(data2[[i]]$data, ph2, data2[[i]]$surv.data, verbose = F, gh = 3, quad = T)),
                         error = function(e) NULL)
  fits1NQ[[i]] <- tryCatch(suppressMessages(EM(data1[[i]]$data, ph1, data1[[i]]$surv.data, verbose = F, gh = 3, quad = F)),
                         error = function(e) NULL)
  fits2NQ[[i]] <- tryCatch(suppressMessages(EM(data2[[i]]$data, ph2, data2[[i]]$surv.data, verbose = F, gh = 3, quad = F)),
                         error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}
fitsQ <- list(fits1Q, fits2Q)#, fits3)#, fits4, fits5, fits6)
fitsNQ <- list(fits1NQ, fits2NQ)

save(fitsQ, file = '~/Downloads/Quad-fits-rustand.RData')
save(fitsNQ, file = '~/Downloads/Noquad-fits-rustand.RData')
