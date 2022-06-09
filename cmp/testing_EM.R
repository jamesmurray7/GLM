rm(list=ls())
source('EM.R')
long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin
disp.formula <- ~time

data <- simData_joint(n = 150)
data <- data$data

test <- EM(long.formula, disp.formula, surv.formula, data, control = list(verbose = T),summax=50)
plot(test$stepmat)



s <- replicate(100, simData_joint(n=150), simplify = F)
fits <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- s[[i]]$data
  test <- tryCatch(EM(long.formula, disp.formula, surv.formula, d, control = list(verbose = T), summax = 100),
                   error = function(e) NULL)
  fits[[i]] <- test
  utils::setTxtProgressBar(pb, i)
}
  
fits[[1]]$coeffs
test$coeffs

