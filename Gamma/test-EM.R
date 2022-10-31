# Gamma
rm(list=ls())
source('EM.R')
test <- simData_joint()
data <- test$data
source('_zzz.R')

fit <- EM(long.formula, surv.formula, 
          data, control = control)

# Intslope ----------------------------------------------------------------
rm(list=ls())
source('EM.R')
test <- simData_joint(D = matrix(c(.2, 0,0,.05), 2, 2),
                      beta = c(0.5, 0.05, 0.1, -0.1))
data <- test$data
source('_zzz.R')
long.formula <- Y~time+cont+bin+(1+time|id)

fit <- EM(long.formula, surv.formula, 
          data, control = control)

my.summary(fit)


# Playing around with shape/scale -----------------------------------------
fn <- function(shape = 2., 
               betas = c(2, -0.1, 0.1, -0.2), 
               D = matrix(c(.25, 0, 0, .05), 2, 2)){
  a <- suppressMessages(simData_joint(shape = shape, D = D, beta = betas))
  d <- a$data; b <- a$tr
  shape = shape
  mu <- lapply(1:250, function(i){
    X <- model.matrix(~time+cont+bin, d[d$id == i,])
    Z <- model.matrix(~time, d[d$id == i,])
    exp(X %*% betas + Z %*% b[i,])
  })
  mu <- do.call(c, mu)
  plot(density(d$Y))
  lines(density(rgamma(1000, shape = shape, scale = mu/shape)), col = 'red')
  return(invisible(1+1))
}

# Larger _range_ of values, but still predominantly located in one central range.
fn(shape = 5,
   betas = c(5, -1, 0.1, -0.2),
   D =  matrix(c(0.5, 0, 0, .1), 2, 2))

test <- simData_joint(shape = 5,
                      beta = c(5, -1, 0.1, -0.2),
                      D =  matrix(c(0.5, 0, 0, .1), 2, 2))$data
fit <- EM(long.formula, surv.formula, 
          test, control = control)

# Much smaller
fn(shape = .5,
   betas = c(2, -0.1, 0.1, -0.2),
   D =  matrix(c(0.1, 0, 0, 0), 2, 2))

test <- simData_joint(shape = .5,
                      beta = c(2, -0.1, 0.1, -0.2),
                      D =  matrix(c(0.1, 0, 0, 0), 2, 2))$data
fit <- EM(Y ~ time + cont + bin + (1|id), surv.formula, 
          test, control = control)
