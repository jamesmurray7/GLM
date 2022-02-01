#' ####
#' Simulations under a quadratic setup
#' ####


# True params
beta <- rbind(c(1, -0.2, 0.01, 0.33, -0.50),
              c(0, -0.5, 0.05, -0.33, 0.50),
              c(3, 0.1, -0.05, 0.5, 0.1))

D <- as.matrix(Matrix::bdiag(replicate(3, diag(c(0.5^2, .2^2, .05^2)), simplify = F)))

# Many fits...
rm(list=ls())
source('EM.R')
# Simulate some different datasets
data1 <- replicate(100, quad.simData(n = 250, ntms = 10, beta = beta, D = D, eta = c(0.05, -0.3), gamma = c(0.50, -0.25, 0.40),
                                      theta0 = -4.5), simplify = F)
data2 <- replicate(100, quad.simData(n = 250, ntms = 6, beta = beta, D = D, eta = c(0.05, -0.3), gamma = c(0.50, -0.25, 0.40),
                                      theta0 = -3), simplify = F)
data3 <- replicate(100, quad.simData(n = 250, ntms = 14, beta = beta, D = D, eta = c(0.05, -0.3), gamma = c(0.50, -0.25, 0.40),
                                     theta0 = -5), simplify = F)

fits1 <- fits2 <- fits3 <- list()
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  ph1 <- coxph(Surv(survtime, status) ~ cont + bin, data1[[i]]$survdat)
  ph2 <- coxph(Surv(survtime, status) ~ cont + bin, data2[[i]]$survdat)
  ph3 <- coxph(Surv(survtime, status) ~ cont + bin, data3[[i]]$survdat)
  fits1[[i]] <- tryCatch(suppressMessages(EM(data1[[i]]$dat, ph1, data1[[i]]$survdat, verbose = F, gh = 3, MVLME = T)),
                         error = function(e) NULL)
  fits2[[i]] <- tryCatch(suppressMessages(EM(data2[[i]]$dat, ph2, data2[[i]]$survdat, verbose = F, gh = 3, MVLME = T)),
                         error = function(e) NULL)
  fits3[[i]] <- tryCatch(suppressMessages(EM(data3[[i]]$dat, ph3, data3[[i]]$survdat, verbose = F, gh = 3, MVLME = T)),
                         error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}
fits <- list(fits1, fits2, fits3)
save(fits, file = '~/Downloads/quadfits.RData')
