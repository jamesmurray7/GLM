rm(list=ls())
source('EM.R')
long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin
disp.formula <- ~time

data <- simData_joint(n = 100, delta = c(-0.5, 0.1), ntms = 15, theta = c(-3, .25))
data <- data$data

test <- EM(long.formula, disp.formula, surv.formula, data, 
           control = list(verbose = T, summax.override = T, tol = 5e-2, gh.nodes = 3),
           summax=100)
plot.stepmat(test)

data <- replicate(50, simData_joint(n=200, delta = c(-0.5, 0.15), ntms = 10, theta = c(-3, 0.25)), simplify = F)
data <- lapply(data, function(x) x$data)

pb <- utils::txtProgressBar(max=50, style = 3)
fits <- vector('list', 50)
for(i in 1:50){
  d <- data[[i]]
  test <- tryCatch(suppressMessages(EM(long.formula, disp.formula, surv.formula, d, 
                      control = list(verbose = F, summax.override = T, tol = 0.05, gh.nodes = 3),
                      summax=100)),
                   error = function(e) NULL)
  fits[[i]] <- test
  utils::setTxtProgressBar(pb, i)
}
save(fits, file = '/data/c0061461/cmpfits.RData')
