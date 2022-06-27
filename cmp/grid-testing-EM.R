# Test one ----------------------------------------------------------------
rm(list=ls())
source('grid-simData.R')
test <- simData_joint(n = 250, delta = c(.8,0), 
                      ntms = 10, theta = c(-3, .25), fup = 3,
                      beta = c(0.0, -0.1, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                      D = matrix(c(0.16, 0, 0, 0.00), 2, 2))

# Testing dispersion of outputted 'Y'.
summary(with(test$data, tapply(Y, id, var))/with(test$data, tapply(Y, id, mean)))

range(test$data$Y)
data <- test$data
rm(to.remove)
N <- 1e4; pete.flag <- T
source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

fit <- EM(long.formula, disp.formula, surv.formula, data, N = N,
          control = list(verbose = T, summax.override = T, tol = 1e-2, gh.nodes = 9, debug = T),
          summax = 100)

# Produce many (intercept-only) datasets ----------------------------------
rm(list=ls())
source('grid-simData.R')
data <- replicate(100,
                  simData_joint(n = 250, delta = c(0.4,-0.2), ntms = 10, theta = c(-3, .25), fup = 3,
                                beta = c(0.0, -0.1, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                                D = matrix(c(0.33, 0, 0, 0.00), 2, 2)),
                  simplify = F)
data <- lapply(data, '[[', 1)
save(data, file = '/data/c0061461/cmp-data-intonly-23-06-22-new.RData')
plot(do.call(rbind, lapply(lapply(data, function(x) x$Y), range))[,2]); abline(h=9.99,lty=5)
rm(to.remove)
# And fit them!
N <- 1e4; pete.flag <- T
source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~time


pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  fit <- tryCatch(suppressMessages(
                EM(long.formula, disp.formula, surv.formula, d, N = N,
                   control = list(verbose = F, summax.override = T, tol = 1e-2, gh.nodes = 3, debug = F),
                   summax=100)
                ), 
                error = function(e) NULL)
  fits[[i]] <- fit
  utils::setTxtProgressBar(pb, i)
}

save(fits, file = '/data/c0061461/cmp-fits-intonly-23-06-22-new.RData')


# Intercept-and-slope -----------------------------------------------------

# Test one //
rm(list=ls())
source('grid-simData.R')
test <- simData_joint(n = 250, delta = c(0.5,-0.1), ntms = 10, theta = c(-3, .25), fup = 3,
                      beta = c(0.0, -0.1, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                      D = matrix(c(0.20, 0, 0, 0.05), 2, 2))

range(test$data$Y)
data <- test$data
rm(to.remove)
N <- 1e4; pete.flag <- T
source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1 + time|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~time

fit <- EM(long.formula, disp.formula, surv.formula, data, N = N, 
          control = list(verbose = T, summax.override = T, tol = 1e-2, gh.nodes = 9, debug = T),
          summax = 100)

# Produce many (intercept-and-slope) datasets -----------------------------
rm(list=ls())
source('grid-simData.R')
data <- replicate(100,
                  simData_joint(n = 250, delta = c(0.5,-0.1), ntms = 10, theta = c(-3, .25), fup = 3,
                                beta = c(0.0, -0.1, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                                D = matrix(c(0.20, 0, 0, 0.05), 2, 2)),
                  simplify = F)
data <- lapply(data, '[[', 1)
save(data, file = '/data/c0061461/cmp-data-intslope-23-06-22-new.RData')
plot(do.call(rbind, lapply(lapply(data, function(x) x$Y), range))[,2]); abline(h=9.99,lty=5)
rm(to.remove)
# And fit them!
N <- 1e4; pete.flag <- T
source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1+time|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~time


pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  fit <- tryCatch(suppressMessages(
    EM(long.formula, disp.formula, surv.formula, d, N = N,
       control = list(verbose = F, summax.override = T, tol = 1e-2, gh.nodes = 3, debug = F),
       summax=100)
  ), 
  error = function(e) NULL)
  fits[[i]] <- fit
  utils::setTxtProgressBar(pb, i)
}

save(fits, file = '/data/c0061461/cmp-fits-intslope-23-06-22-new.RData')

