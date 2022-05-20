# Comparing intercept and slope fits via 
rm(list=ls())
library(tidyverse)
source('EM.R')
source('Simulations/funs.R')
dataDir <- paste0(getwd(), '/Simulations/fits')

# Tabulating estimates ----------------------------------------------------
#' Gaussian ----
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^gaussian')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}

# -- K = 1 --
G1 <- left_join(tabulate_wrapper(fit1, 'aEM', 'gaussian', 1),
                tabulate_wrapper(fit1.jML, 'jML', 'gaussian', 1), c('param', 'target')) %>% 
      left_join(., tabulate_wrapper(fit1.JMb, 'JMb', 'gaussian', 1), c('param', 'target')) %>% 
      left_join(., tabulate_wrapper(fit1.INLA, 'INLA', 'gaussian', 1), c('param', 'target'))

# -- K = 2 --
G2 <- left_join(tabulate_wrapper(fit2, 'aEM', 'gaussian', 2),
                tabulate_wrapper(fit2.jML, 'jML', 'gaussian', 2), c('param', 'target')) %>% 
  left_join(., tabulate_wrapper(fit2.JMb, 'JMb', 'gaussian', 2), c('param', 'target')) %>% 
  left_join(., tabulate_wrapper(fit2.INLA, 'INLA', 'gaussian', 2), c('param', 'target'))

# -- K = 3 --
G3 <- left_join(tabulate_wrapper(fit3, 'aEM', 'gaussian', 3),
                tabulate_wrapper(fit3.jML, 'jML', 'gaussian', 3), c('param', 'target')) %>% 
  left_join(., tabulate_wrapper(fit3.JMb, 'JMb', 'gaussian', 3), c('param', 'target')) %>% 
  left_join(., tabulate_wrapper(fit3.INLA, 'INLA', 'gaussian', 3), c('param', 'target'))

#' Poisson ----
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^poisson')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}

P1 <- left_join(tabulate_wrapper(fit1, 'aEM', 'poisson', 1),
                tabulate_wrapper(fit1.JMb, 'JMb', 'poisson', 1), c('param', 'target')) %>% 
  left_join(., tabulate_wrapper(fit1.INLA, 'INLA', 'poisson', 1), c('param', 'target'))

# -- K = 2 --
P2 <- left_join(tabulate_wrapper(fit2, 'aEM', 'poisson', 2),
                tabulate_wrapper(fit2.JMb, 'JMb', 'poisson', 2), c('param', 'target')) %>% 
  left_join(., tabulate_wrapper(fit2.INLA, 'INLA', 'poisson', 2), c('param', 'target'))

# -- K = 3 --
P3 <- left_join(tabulate_wrapper(fit3, 'aEM', 'poisson', 3),
                   tabulate_wrapper(fit3.JMb, 'JMb', 'poisson', 3), c('param', 'target')) %>% 
  left_join(., tabulate_wrapper(fit3.INLA, 'INLA', 'poisson', 3), c('param', 'target'))

