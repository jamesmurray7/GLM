#' #########
#' ADNI-Fitter
#' Fitting both EM and JMbayes2 to ADNI.
#' #########

rm(list=ls())
source('EM.R')

# Binary >= 80% -----------------------------------------------------------
load('ADNI/ADNI80.RData')
data <- ADNI$data
survdata <- ADNI$survdata
ph <- ADNI$ph
fit80 <- EM(data, ph, survdata, gh = 3)
save(fit80, file = 'ADNI/fit-ADNI80.RData')


# Binary >= 100% ----------------------------------------------------------
load('ADNI/ADNI100.RData')
data <- ADNI$data
survdata <- ADNI$survdata
ph <- ADNI$ph
fit100 <- EM(data, ph, survdata, gh = 3)
save(fit100, file = 'ADNI/fit-ADNI100.RData') 


# JMbayes2 ----------------------------------------------------------------
rm(list=ls())
library(JMbayes2)

# >= 80%
load('ADNI/ADNI80.RData')
data <- ADNI$data
survdata <- ADNI$survdata
ph <- ADNI$ph

m1 <- lme(Y.1 ~ cont + time + bin, random = ~ time | id, data = data)
m2 <- mixed_model(Y.2 ~ cont + time + bin, random = ~ time | id, data = data, family = binomial())
m3 <- mixed_model(Y.3 ~ cont + time + bin, random = ~ time | id, data = data, family = poisson())

M <- list(m1, m2, m3)

jm80 <- jm(ph, M, time_var = 'time',
           n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)
save(jm80, file = 'ADNI/jm-ADNI80.RData')

# >= 100%
load('ADNI/ADNI100.RData')
data <- ADNI$data
survdata <- ADNI$survdata
ph <- ADNI$ph

m1 <- lme(Y.1 ~ cont + time + bin, random = ~ time | id, data = data)
m2 <- mixed_model(Y.2 ~ cont + time + bin, random = ~ time | id, data = data, family = binomial())
m3 <- mixed_model(Y.3 ~ cont + time + bin, random = ~ time | id, data = data, family = poisson())

M <- list(m1, m2, m3)

jm100 <- jm(ph, M, time_var = 'time',
            n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)
save(jm100, file = 'ADNI/jm-ADNI100.RData')


