#' ####
#' Gamma vs Gaussian
#' ####

rm(list=ls())
source('EM.R')

# Larger distns -----------------------------------------------------------
test <- replicate(100,
                  simData_joint(shape = 5,
                                beta = c(5, -1, 0.1, -0.2),
                                D = matrix(c(0.5, 0, 0, .1), 2, 2))$data,
                  simplify = F)

pb <- utils::txtProgressBar(max = 100, style = 3)
fits.Gamma <- vector('list', 100)
for(i in 1:100){
  d <- test[[i]]
  fits.Gamma[[i]] <- suppressMessages(
    EM(long.formula = Y ~ time + cont + bin + (1 + time|id), 
       surv.formula = Surv(survtime, status) ~ bin, 
       d)
  )
  utils::setTxtProgressBar(pb, i)
}
save(fits.Gamma, file = '/data/c0061461/gammagauss/Gamma_bigger.RData')
save(test, file = '/data/c0061461/gammagauss/bigger_data.RData')

test <- replicate(100,
                  simData_joint(shape = .5,
                                beta = c(2, -0.1, 0.1, -0.2),
                                D =  matrix(c(0.1, 0, 0, 0), 2, 2))$data,
                  simplify = F)

pb <- utils::txtProgressBar(max = 100, style = 3)
fits.Gamma <- vector('list', 100)
for(i in 1:100){
  d <- test[[i]]
  fits.Gamma[[i]] <- suppressMessages(
    EM(long.formula = Y ~ time + cont + bin + (1|id), 
       surv.formula = Surv(survtime, status) ~ bin, 
       d)
  )
  utils::setTxtProgressBar(pb, i)
}
save(fits.Gamma, file = '/data/c0061461/gammagauss/Gamma_smaller.RData')
save(test, file = '/data/c0061461/gammagauss/smaller_data.RData')



# Multi-test/Gaussian -----------------------------------------------------
rm(list=ls())
setwd('../Multi-test/')
source('EM.R')

# Smaller
load('/data/c0061461/gammagauss/smaller_data.RData')
fits.Gauss <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d <- test[[i]]
  ff <- tryCatch(suppressMessages(
    EM(long.formula = list(Y ~ time + cont + bin + (1|id)), 
       surv.formula = Surv(survtime, status) ~ bin, 
       d, family= list('gaussian'))
  ), error = function(e) NULL)
  fits.Gauss[[i]] <- ff
  utils::setTxtProgressBar(pb, i)
}
save(fits.Gauss, file = '/data/c0061461/Gaussian_smaller.RData')
# Larger
load('/data/c0061461/gammagauss/bigger_data.RData')
fits.Gauss <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d <- test[[i]]
  ff <- tryCatch(suppressMessages(
    EM(long.formula = list(Y ~ time + cont + bin + (1+time|id)), 
       surv.formula = Surv(survtime, status) ~ bin, 
       d, family= list('gaussian'))
  ), error = function(e) NULL)
  fits.Gauss[[i]] <- ff
  utils::setTxtProgressBar(pb, i)
}
save(fits.Gauss, file = '/data/c0061461/Gaussian_bigger.RData')

