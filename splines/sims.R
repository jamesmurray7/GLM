#' ####
#' Simulations under a quadratic setup
#' ####

rm(list=ls())
# True params -> underlying Quadratic set-up 
beta <- rbind(c(1, -0.2, 0.01, 0.33, -0.50),
              c(0, -0.5, 0.05, -0.33, 0.50),
              c(3, 0.1, -0.05, 0.5, 0.1))

D <- as.matrix(Matrix::bdiag(replicate(3, diag(c(0.5^2, .2^2, .05^2)), simplify = F)))

# Many fits...
source('EM.R')
# Simulate some different datasets
data1 <- replicate(100, poly.simData(n = 200, ntms = 10, beta = beta, D = D, theta0 = -4.5), simplify = F)
# Shorter profile availability
data2 <- replicate(100, poly.simData(n = 200, ntms = 6, beta = beta, D = D, theta0 = -2.75), simplify = F)

fits1 <- fits2 <- fits3 <- fits4 <- list()
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  dat1 <- data1[[i]]$dat; survdata1 <- data1[[i]]$survdat
  ph1 <- coxph(Surv(survtime, status) ~ cont + bin, survdata1)
  dat2 <- data2[[i]]$dat; survdata2 <- data2[[i]]$survdat
  ph2 <- coxph(Surv(survtime, status) ~ cont + bin, survdata2)
  fits1[[i]] <- tryCatch(suppressMessages(EM(dat1, ph1, survdata1, verbose = F, gh = 3, MVLME = T, degree = 3)),
                         error = function(e) NULL)
  fits2[[i]] <- tryCatch(suppressMessages(EM(dat1, ph1, survdata1, verbose = F, gh = 3, MVLME = T, degree = 6)),
                         error = function(e) NULL)
  fits3[[i]] <- tryCatch(suppressMessages(EM(dat2, ph2, survdata2, verbose = F, gh = 3, MVLME = T, degree = 3)),
                         error = function(e) NULL)
  fits4[[i]] <- tryCatch(suppressMessages(EM(dat2, ph2, survdata2, verbose = F, gh = 3, MVLME = T, degree = 6)),
                         error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}
fits <- list(fits1, fits2, fits3, fits4)
save(fits, file = '~/Downloads/splinefits.RData')
