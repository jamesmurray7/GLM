rm(list=ls())
source('EM.R')
dataDir <- paste0(getwd(), '/Simulations/data')

# Define fitting functions ------------------------------------------------
emfit <- function(data, k, family){
  long.formulas <- vector('list', k)
  for(kk in 1:k) long.formulas[[kk]] <- as.formula(paste0('Y.', kk, '~ time + cont + bin + (1 + time|id)'))
  surv.formula <- Surv(survtime, status) ~ bin
  families <- as.list(rep(family, k))
  fit <- tryCatch(suppressMessages(
    EM(long.formulas, surv.formula, data = data, family = families,
       control = list(hessian = 'auto'))),
    error = function(e) NULL
  )
  fit
}

joineRMLfit <- function(data, k){
  long <- rand <- vector('list', k)
  for(kk in 1:k){
    long[[kk]] <- as.formula(paste0('Y.', kk, '~time + cont + bin'))
    rand[[kk]] <- ~ 1 + time | id
  }
  fit <- suppressMessages(joineRML::mjoint(
    formLongFixed = long,
    formLongRandom = rand,
    formSurv = Surv(survtime, status) ~ bin,
    timeVar = 'time', data = data, 
    control = list(
      convCrit = 'rel',
      tol.em = 5e-3,
      type = 'sobol', tol2 = 1e-2
    )
  ))
  summary(fit)
}

JMbayes2fit <- function(data, k, family){
  M <- vector('list', k)
  if(family == 'gaussian'){
    for(kk in 1:k){
      long <- as.formula(paste0('Y.', kk, '~time + cont + bin'))
      random <- ~time|id
      M[[kk]] <- lme(fixed = long, random = random, 
                     data = data)
    }
  }else{
    for(kk in 1:k){
      long <- as.formula(paste0('Y.', kk, '~time + cont + bin'))
      random <- ~time|id
      M[[kk]] <- mixed_model(fixed = long, random = random,
                             data = data, family = family)
    }
  }
  
  survdata <- data[!duplicated(data[, 'id']), ]
  ph <- coxph(Surv(survtime, status) ~ bin, survdata)
  
  fit <- tryCatch(
    jm(ph, M, time_var = 'time', id_var = 'id', data_Surv = survdata,
            control = list(
              n_chains = 1,
              n_iter = 1500,
              n_burnin = 500, cores = 1
            )),error=function(e) NULL)
  if(!is.null(fit)) rtn <- summary(fit) else rtn <- NULL
  rtn
}

.loader <- function(file){
  assign('data', get(load(file)))
  lapply(data, function(x) x$data)
}

# K=1 ---------------------------------------------------------------------
data1 <- .loader('Simulations/data/gaussianK-1.RData')
fit1 <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data1[[i]]
  fit1[[i]] <- emfit(d, 1, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit1, file =  'Simulations/fits/gaussianK-1.RData')

# K=2 ---------------------------------------------------------------------
data2 <- .loader('Simulations/data/gaussianK-2.RData')
fit2 <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data2[[i]]
  fit2[[i]] <- emfit(d, 2, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit2, file = 'Simulations/fits/gaussianK-2.RData')

# K=3 ---------------------------------------------------------------------
data3 <- .loader('Simulations/data/gaussianK-3.RData')
fit3 <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data3[[i]]
  fit3[[i]] <- emfit(d, 3, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit3, file = 'Simulations/fits/gaussianK-3.RData')


# joineRML fits -----------------------------------------------------------
# K=1 joineRML ------------------------------------------------------------
data1 <- .loader('Simulations/data/gaussianK-1.RData')
fit1.jML <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data1[[i]]
  fit1.jML[[i]] <- joineRMLfit(d, 1)
  utils::setTxtProgressBar(pb, i)
}
save(fit1.jML, file =  'Simulations/fits/gaussianK-1_joineRML.RData')

# K=2 joineRML ------------------------------------------------------------
data2 <- .loader('Simulations/data/gaussianK-2.RData')
fit2.jML <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data2[[i]]
  fit2.jML[[i]] <- joineRMLfit(d, 2)
  utils::setTxtProgressBar(pb, i)
}
save(fit2.jML, file =  'Simulations/fits/gaussianK-2_joineRML.RData')

# K=3 joineRML ------------------------------------------------------------
data3 <- .loader('Simulations/data/gaussianK-3.RData')
fit3.jML <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data3[[i]]
  fit3.jML[[i]] <- joineRMLfit(d, 3)
  utils::setTxtProgressBar(pb, i)
}
save(fit3.jML, file =  'Simulations/fits/gaussianK-3_joineRML.RData')


# JMbayes2 ----------------------------------------------------------------
# K=1 JMbayes 2------------------------------------------------------------
data1 <- .loader('Simulations/data/gaussianK-1.RData')
fit1.JMb <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data1[[i]]
  fit1.JMb[[i]] <- JMbayes2fit(d, 1, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit1.JMb, file = 'Simulations/fits/gaussianK-1_JMbayes2.RData')

# K=2 JMbayes 2------------------------------------------------------------
data2 <- .loader('Simulations/data/gaussianK-2.RData')
fit2.JMb <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data2[[i]]
  fit2.JMb[[i]] <- JMbayes2fit(d, 2, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit2.JMb, file = 'Simulations/fits/gaussianK-2_JMbayes2.RData')

# K=3 JMbayes 2------------------------------------------------------------
data3 <- .loader('Simulations/data/gaussianK-3.RData')
fit3.JMb <- vector('list', 100)
pb <- utils::txtProgressBar(max=100,style=3)
for(i in 1:100){
  d <- data3[[i]]
  fit3.JMb[[i]] <- JMbayes2fit(d, 3, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}

save(fit3.JMb, file = 'Simulations/fits/gaussianK-3_JMbayes2.RData')

# END GAUSSIAN --------
# Poisson -----------------------------------------------------------------
rm(data1, data2, data3)
data1 <- .loader('Simulations/data/poissonK-1.RData')
data2 <- .loader('Simulations/data/poissonK-2.RData')
data3 <- .loader('Simulations/data/poissonK-3.RData')
fit1 <- fit2 <- fit3 <- fit1.JMb <- fit2.JMb <- fit3.JMb <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d1 <- data1[[i]]; d2 <- data2[[i]]; d3 <- data3[[i]]
  fit1[[i]] <- emfit(d1, 1, 'poisson')
  fit2[[i]] <- emfit(d2, 2, 'poisson')
  fit3[[i]] <- emfit(d3, 3, 'poisson')
  utils::setTxtProgressBar(pb, i)
}
save(fit1, file = 'Simulations/fits/poissonK-1.RData')
save(fit2, file = 'Simulations/fits/poissonK-2.RData')
save(fit3, file = 'Simulations/fits/poissonK-3.RData')

# JMbayes2
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d1 <- data1[[i]]; d2 <- data2[[i]]; d3 <- data3[[i]]
  fit1.JMb[[i]] <- JMbayes2fit(d1, 1, 'poisson')
  fit2.JMb[[i]] <- JMbayes2fit(d2, 2, 'poisson')
  fit3.JMb[[i]] <- JMbayes2fit(d3, 3, 'poisson')
  utils::setTxtProgressBar(pb, i)
}
save(fit1.JMb, file = 'Simulations/fits/poissonK-1_JMbayes2.RData')
save(fit2.JMb, file = 'Simulations/fits/poissonK-2_JMbayes2.RData')
save(fit3.JMb, file = 'Simulations/fits/poissonK-3_JMbayes2.RData')

# END COUNTS --------
# Binary ------------------------------------------------------------------
