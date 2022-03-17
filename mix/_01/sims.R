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

#' ####
#' 23/2/22
#' --------
#' Simulating 50 sets for comparison across univariate GLMM fits and multivariate ones
#' ####

simdata <- replicate(50, simData_joint(), simplify = F)
save(simdata, file = 'simdata.RData')

#' #####
#' 23/2/22 
#' -----
#' END
#' #####

# One fit...
dd <- simData_joint(thetaDisp = 1.5)
ph <- coxph(Surv(survtime, status) ~ cont + bin, dd$surv.data)
fitP <- EM(dd$data, ph, dd$surv.data, verbose = T, gh =3, nb = F)
fitNB <- EM(dd$data, ph, dd$surv.data, verbose = T, gh =3, nb = T)

# Many fits...
rm(list=ls())
source('EM.R')
diag(true.D) <- c(.25, .06, .50, .04, .25, .05)
true.D <- as.matrix(Matrix::nearPD(true.D)$mat)
# Simulate some different datasets
data1 <- replicate(100, simData_joint(n = 250), simplify = F)
data2 <- replicate(100, simData_joint(n = 500), simplify = F)
data3 <- replicate(100, simData_joint(n = 1000), simplify = F)
# data3 <- replicate(100, simData_joint(thetaDisp = 2), simplify = F)
# data3 <- replicate(100, simData_joint(ntms = 15, theta = c(-6, 0.25)), simplify = F)

fits1 <- fits2 <- fits3 <- list() # fits4 <- fits5 <- fits6 <- list()
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  ph1 <- coxph(Surv(survtime, status) ~ cont + bin, data1[[i]]$surv.data)
  ph2 <- coxph(Surv(survtime, status) ~ cont + bin, data2[[i]]$surv.data)
  ph3 <- coxph(Surv(survtime, status) ~ cont + bin, data3[[i]]$surv.data)
  # ph4 <- coxph(Surv(survtime, status) ~ cont + bin, data4[[i]]$surv.data)
  # ph5 <- coxph(Surv(survtime, status) ~ cont + bin, data5[[i]]$surv.data)
  # ph6 <- coxph(Surv(survtime, status) ~ cont + bin, data6[[i]]$surv.data)
  fits1[[i]] <- tryCatch(suppressMessages(EM(data1[[i]]$data, ph1, data1[[i]]$surv.data, verbose = F, gh = 3)),
                         error = function(e) NULL)
  fits2[[i]] <- tryCatch(suppressMessages(EM(data2[[i]]$data, ph2, data2[[i]]$surv.data, verbose = F, gh = 3)),
                         error = function(e) NULL)
  fits3[[i]] <- tryCatch(suppressMessages(EM(data3[[i]]$data, ph3, data3[[i]]$surv.data, verbose = F, gh = 3)),
                         error = function(e) NULL)
  # fits4[[i]] <- tryCatch(suppressMessages(EM(data4[[i]]$data, ph4, data4[[i]]$surv.data, verbose = F)),
  #                        error = function(e) NULL)
  # fits5[[i]] <- tryCatch(suppressMessages(EM(data5[[i]]$data, ph5, data5[[i]]$surv.data, verbose = F)),
  #                        error = function(e) NULL)
  # fits6[[i]] <- tryCatch(suppressMessages(EM(data6[[i]]$data, ph6, data6[[i]]$surv.data, verbose = F)),
  #                        error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}
fits <- list(fits1, fits2, fits3)#, fits4, fits5, fits6)
save(fits, file = '~/Downloads/mixfits-2.RData')
