# Test one ----------------------------------------------------------------
rm(list=ls())
source('grid-simData.R')
test <- simData_joint(n = 250, delta = c(.8,0), 
                      ntms = 10, theta = c(-2, .1), fup = 3,
                      beta = c(0.0, -0.1, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                      D = matrix(c(0.25, 0, 0, 0.00), 2, 2))

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
          control = list(verbose = T, summax.override = T, tol = 1e-2, gh.nodes = 3, debug = T),
          summax = 100)

# Produce many (intercept-only) datasets ----------------------------------
rm(list=ls())
source('grid-simData.R')
data <- replicate(100,
                  simData_joint(n = 250, delta = c(0.8,0), ntms = 10, theta = c(-3, .25), fup = 3,
                                beta = c(0.0, -0.1, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                                D = matrix(c(0.16, 0, 0, 0.00), 2, 2)),
                  simplify = F)
data <- lapply(data, '[[', 1)
save(data, file = '/data/c0061461/cmp-data-intonly-28-06-22-new.RData')
plot(do.call(rbind, lapply(lapply(data, function(x) x$Y), range))[,2]); abline(h=9.99,lty=5)
# plot dispersion.
plot(
  do.call(c, lapply(data, function(x) tapply(x$Y, x$id, var)))/
  do.call(c, lapply(data, function(x) tapply(x$Y, x$id, mean)))
)
rm(to.remove)
# And fit them!
N <- 1e4; pete.flag <- T
source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1


pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  fit <- tryCatch(suppressMessages(
                EM(long.formula, disp.formula, surv.formula, d, N = N,
                   control = list(verbose = F, summax.override = T, tol = 1e-2, gh.nodes = 9, debug = F),
                   summax=100)
                ), 
                error = function(e) NULL)
  fits[[i]] <- fit
  utils::setTxtProgressBar(pb, i)
}

save(fits, file = '/data/c0061461/cmp-fits-intonly-28-06-22-with-QUAD.RData')

# Intercept-only time-varying nu ------------------------------------------
# Test one //
rm(list=ls())
source('grid-simData.R')
test <- simData_joint(n = 250, delta = c(0.8,-0.2), ntms = 10, theta = c(-3, .25), fup = 3,
                      beta = c(-0.5, -0.1, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                      D = matrix(c(0.16, 0, 0, 0.00), 2, 2))
range(test$data$Y)
data <- test$data
rm(to.remove)
N <- 1e4; pete.flag <- T
summary(with(test$data, tapply(Y, id, var))/with(test$data, tapply(Y, id, mean))) # check var(Y) < mean(Y)

source('grid-EM.R')
long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~time
Longit.inits(long.formula, disp.formula, data) # Make sure it doesn't error; keeps happening!
fit <- glmmTMB(long.formula, family = poisson, ziformula = ~1, dispformula = ~time,
               data = data, control = glmmTMBControl(optimizer='optim', optArgs =  list(method="BFGS") ))
glmmTMB::VarCorr(fit)$c$id; dimD <- dim(D) # Quite badly wrong...

fit <- EM(long.formula, disp.formula, surv.formula, data, N = N, 
          control = list(verbose = T, summax.override = T, tol = 1e-2, gh.nodes = 3, debug = T),
          summax = 100)

# Produce many (time-varying-nu) datasets ---------------------------------
rm(list=ls())
source('grid-simData.R')
data <- replicate(100,
                  simData_joint(n = 250, delta = c(0.8,-0.2), ntms = 10, theta = c(-3, .25), fup = 3,
                                beta = c(-0.5, -0.1, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                                D = matrix(c(0.16, 0, 0, 0.00), 2, 2)),
                  simplify = F)
data <- lapply(data, '[[', 1)
save(data, file = '/data/c0061461/cmp-data-intonly-timevarnu.RData')
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

save(fits, file = '/data/c0061461/cmp-fits-intonly-timevarnu.RData')

