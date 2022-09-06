setwd('~/Documents/GLMM/mpcmp/hybrid2')
source('EM.R')
N <- 50
tests <- replicate(N, simData_joint(n = 250, delta = c(.4,0), 
                       ntms = 10, theta = c(-2, .1), fup = 3,
                       beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                       D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data, simplify = F)

save(tests, file = '/data/c0061461/a_bunch_of_tests_4.RData')

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

fits1 <- vector('list', N)
pb <- utils::txtProgressBar(max = N, style = 3)
for(i in 1:N){
  d <- tests[[i]]
  f1 <- tryCatch(suppressMessages(EM(long.formula, disp.formula, surv.formula, d, summax = 100,
                   control = list(verbose = F), grid.summax = 'aaa')),
                 error = function(e) NULL
  )
  # f2 <- tryCatch(suppressMessages(EM(long.formula, disp.formula, surv.formula, d, summax = 100,
  #                   control = list(verbose = F), grid.summax = 'aaa')),
  #                error = function(e) NULL
  # )
  fits1[[i]] <- f1
  utils::setTxtProgressBar(pb, i)
}

save(fits1, file = '/data/c0061461/hybrid2_fits4.RData')


# no-grids ----------------------------------------------------------------
rm(list=ls())
N <- 50
setwd('../no-grids')
source('EM.R')
load('/data/c0061461/a_bunch_of_tests_4.RData')
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

fits2 <- vector('list', N)
pb <- utils::txtProgressBar(max = N, style = 3)
for(i in 1:N){
  d <- tests[[i]]
  f2 <- tryCatch(suppressMessages(EM(long.formula, disp.formula, surv.formula, d)),
                 error = function(e) NULL
  )
  fits2[[i]] <- f2
  utils::setTxtProgressBar(pb, i)
}

save(fits2, file = '/data/c0061461/no-grids_fits4.RData')
