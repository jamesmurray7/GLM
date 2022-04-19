#' ###
#' Individual univariate testing
#' ###


rm(list=ls())
source('EM WIP.R')

ff <- function(){
  rm(list=setdiff(ls(), 'ff'))
  source('EM WIP.R')
}

# Gaussian ----------------------------------------------------------------
d <- simData()
data <- d$data
long.formula <- Y ~ time + cont + bin + (1 + time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin

EM(long.formula, surv.formula, data, gaussian, control = list(verbose = T))


# Poisson -----------------------------------------------------------------
ff()
d <- simData(family = poisson)
data <- d$data
long.formula <- Y ~ time + cont + bin + (1 + time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin

EM(long.formula, surv.formula, data, poisson, control = list(verbose = T))

# Binomial ----------------------------------------------------------------
ff()
d <- simData(family = binomial)
data <- d$data
long.formula <- Y ~ time + cont + bin + (1 + time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin

EM(long.formula, surv.formula, data, binomial, control = list(verbose = T))


# Negative binomial -------------------------------------------------------
ff()
d <- simData(family = "negative.binomial", disp = 1.5)
data <- d$data
long.formula <- Y ~ time + cont + bin + (1 + time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin

EM(long.formula, surv.formula, data, "negative.binomial", control = list(verbose = T))

# Messing around with different forms
ff()
d <- simData(family = "negative.binomial", disp = 1.5)
data <- d$data
long.formula <- Y ~ time + (1|id)
surv.formula <- Surv(survtime, status) ~ bin

EM(long.formula, surv.formula, data, "negative.binomial", control = list(verbose = T))
