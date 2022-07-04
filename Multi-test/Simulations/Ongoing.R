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
       control = list(hessian = 'manual'))),
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
                     data = data, method = 'ML',
                     control = nlme::lmeControl(opt = "optim"))
    }
  }else{
    for(kk in 1:k){
      long <- as.formula(paste0('Y.', kk, '~time + cont + bin'))
      random <- ~time|id
      TMBfit <- glmmTMB(as.formula(paste0('Y.',kk,'~time+cont+bin+(1+time|id)')),
                        data = data, family = family)
      inits <- list(
        betas = glmmTMB::fixef(TMBfit)$cond,
        D = matrix(glmmTMB::VarCorr(TMBfit)$cond$id, 2, 2)
      )
      M[[kk]] <- mixed_model(fixed = long, random = random,
                             data = data, family = family,
                             initial_values = inits)
    }
  }
  
  survdata <- data[!duplicated(data[, 'id']), ]
  ph <- coxph(Surv(survtime, status) ~ bin, survdata)
  
  fit <- tryCatch(
    jm(ph, M, time_var = 'time', id_var = 'id', data_Surv = survdata,
       control = list(
         n_chains = 3,
         n_iter = 1500,
         n_burnin = 500, cores = 3
       )),error=function(e) NULL)
  if(!is.null(fit)){
    s <- summary(fit)
    rtn <- list()
    rtn$Outcome1 <- s$Outcome1
    rtn$Outcome2 <- s$Outcome2
    rtn$Outcome3 <- s$Outcome3
    rtn$survival <- s$Survival
    rtn$comp.time <- s$time
    return(rtn)
  }else{ 
    return(NA)
  }
}

INLAfit <- function(data, k, family, cl = T){
  long.formulas <- vector('list', k)
  for(kk in 1:k) long.formulas[[kk]] <- as.formula(paste0('Y.', kk, '~ time + cont + bin + (1 + time|id)'))
  assocs <- as.list(rep('SRE', k))
  survdata <- data[!duplicated(data[, 'id']), c('survtime', 'status', 'id')]
  data <- data[,c('id', 'time', 'cont', 'bin', paste0('Y.', 1:k))]
  fail <<- inla.surv(time = c(survdata$survtime), event = c(survdata$status))
  JMINLA <- suppressMessages(joint(
    formLong = long.formulas,
    formSurv = list(fail ~ bin),
    dataLong = data,
    id = 'id', timeVar = 'time', corLong = cl,
    family = rep(family, k), basRisk = 'rw2',
    assoc = assocs,
    control = list(int.strategy='eb',
                   priorRandom = list(r = 2*k, R = 1),
                   priorAssoc = list(mean = 0, prec = 0.16),
                   priorFixed = list(mean = 0, prec = 0.16, 
                                     mean.intercept = 0,
                                     prec.intercept = 0.16))
  ))
  s <- summary(JMINLA)
  longi <- do.call(rbind, s$FixedEff)
  survi <- s$SurvEff[[1]]
  gammas <- s$AssocLS
  time <- s$cpu.used   # Just running time?
  
  list(
    fixed = longi,
    survival = survi,
    gamma = gammas,
    comp.time = time
  )
}

.loader <- function(file){
  assign('data', get(load(file)))
  lapply(data, function(x) x$data)
}

library(JMbayes2)
library(INLA)
library(INLAjoint)
inla.setOption(inla.mode="experimental")
# Gaussian {2,3} JMbayes2 --------------------------------------------------
data2 <- .loader('Simulations/data/gaussianK-2.RData')
data3 <- .loader('Simulations/data/gaussianK-3.RData')
fit2.JMb <- fit3.JMb <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d2 <- data2[[i]]; d3 <- data3[[i]]
  fit2.JMb[[i]] <- JMbayes2fit(d2, 2, 'gaussian')
  fit3.JMb[[i]] <- JMbayes2fit(d3, 3, 'gaussian')
  utils::setTxtProgressBar(pb, i)
}
save(fit2.JMb, file = 'Simulations/fits/gaussianK-2_JMbayes2.RData')
save(fit3.JMb, file = 'Simulations/fits/gaussianK-3_JMbayes2.RData')

rm(data2, data3, fit2.JMb, fit3.JMb)


# Binomial {1, 2} JMbayes2 ------------------------------------------------
data1 <- .loader('Simulations/data/binomialK-1.RData')
data2 <- .loader('Simulations/data/binomialK-2.RData')
fit1.JMb <- fit2.JMb <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)

for(i in 1:100){
  d1 <- data1[[i]]; d2 <- data2[[i]]
  fit1.JMb[[i]] <- JMbayes2fit(d1, 1, 'binomial')
  fit2.JMb[[i]] <- JMbayes2fit(d2, 2, 'binomial')
  utils::setTxtProgressBar(pb, i)
}
save(fit1.JMb, file = 'Simulations/fits/binomialK-1_JMbayes2.RData')
save(fit2.JMb, file = 'Simulations/fits/binomialK-2_JMbayes2.RData')

rm(data1, data2, fit1.JMb, fit2.JMb)

# Poisson {3} INLA --------------------------------------------------------
data3 <- .loader('Simulations/data/poissonK-3.RData')
fit3.INLA <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 91:100){
  d3 <- data3[[i]]
  fit3.INLA[[i]] <- tryCatch(suppressMessages(INLAfit(d3, 3, 'poisson')), error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}
save(fit3.INLA, file = 'Simulations/fits/poissonK-3_INLA.RData')

rm(data3, fit3.INLA)

# Binomial {3} INLA -------------------------------------------------------
data3 <- .loader('Simulations/data/binomialK-3.RData')
fit3.INLA <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d3 <- data3[[i]]
  fit3.INLA[[i]] <- tryCatch(suppressMessages(INLAfit(d3, 3, 'binomial')), error = function(e) NULL)
  utils::setTxtProgressBar(pb, i)
}
save(fit3.INLA, file = 'Simulations/fits/binomialK-3_INLA.RData')

# Binomial {3} JMbayes2 ---------------------------------------------------
data3 <- .loader('Simulations/data/binomialK-3.RData')
fit3.JMb <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d3 <- data3[[i]]
  fit3.JMb[[i]] <- JMbayes2fit(d3, 3, 'binomial')
  utils::setTxtProgressBar(pb, i)
}
save(fit3.JMb, file = 'Simulations/fits/binomialK-3_JMbayes2.RData')
