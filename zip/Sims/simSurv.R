setwd('~/Documents/GLMM/zip/')
source('EM.R')

beta <- c(.5, -0.2, 0.1, 0.2)
alpha <- c(-1, -0.1)
D <- diag(c(0.5, .1))
n <- 250
ntms <- 10
gamma <- -0.5
data1 <- replicate(100, 
                   simData_zip_joint(n = 250, ntms = 5, beta, alpha, D, gamma = gamma, surv.eta = c(0.05, -0.3), theta = c(-3, .2)), simplify = F)
data2 <- replicate(100, 
                   simData_zip_joint(n = 250, ntms = 10, beta, alpha, D, gamma = gamma, surv.eta = c(0.05, -0.3), theta = c(-2.5, .2)), simplify = F)
data3 <- replicate(100, 
                   simData_zip_joint(n = 250, ntms = 15, beta, alpha, D, gamma = gamma, surv.eta = c(0.05, -0.3), theta = c(-2.5, .2)), simplify = F)
data4 <- replicate(100, 
                   simData_zip_joint(n=250, ntms=15, beta, alpha, D, gamma = gamma, surv.eta = c(0.05, -0.3), theta = c(-3, .2)), simplify = F)

pb <- utils::txtProgressBar(max = 100, style = 3)
fits1 <- fits2 <- fits3 <- fits4 <- list() # <- fits5 <- fits6 <- fits7 <- list();
for(m in 1:100){
  ph1 <- coxph(Surv(survtime, status) ~ cont  + bin, data = data1[[m]]$surv.data)
  ph2 <- coxph(Surv(survtime, status) ~ cont  + bin, data = data2[[m]]$surv.data)
  ph3 <- coxph(Surv(survtime, status) ~ cont  + bin, data = data3[[m]]$surv.data)
  ph4 <- coxph(Surv(survtime, status) ~ cont  + bin, data = data4[[m]]$surv.data)
  fits1[[m]] <- tryCatch(suppressMessages(EM(data1[[m]]$data, ph1, data1[[m]]$surv.data, gh = 9)), 
                           error = function(e) NULL)
  fits2[[m]] <- tryCatch(suppressMessages(EM(data2[[m]]$data, ph2, data2[[m]]$surv.data, gh = 9)), 
                           error = function(e) NULL)
  fits3[[m]] <- tryCatch(suppressMessages(EM(data3[[m]]$data, ph3, data3[[m]]$surv.data, gh = 9)), 
                         error = function(e) NULL)
  fits4[[m]] <- tryCatch(suppressMessages(EM(data4[[m]]$data, ph4, data4[[m]]$surv.data, gh = 9)), 
                         error = function(e) NULL)
  utils::setTxtProgressBar(pb, m)
}
fits <- list(fits1, fits2, fits3, fits4)#, fits5, fits6, fits7)
save(fits,file='~/Downloads/test.RData')
