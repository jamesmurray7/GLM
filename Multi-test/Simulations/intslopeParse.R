# Comparing intercept and slope fits across fitting methods.
rm(list=ls())
library(tidyverse)
source('EM.R')
dataDir <- paste0(getwd(), '/Simulations/fits')

# Gaussians ---------------------------------------------------------------
for(file in dir(dataDir, '^gaussian')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}

# Functions for elapsed time ----------------------------------------------
elapsed.aEM <- function(fit, what = 'EM'){
  if(what == 'EM') return(fit$EMtime + fit$postprocess.time)
  else if(what == 'total') return(fit$totaltime)
}

elapsed.jML <- function(fit, what = 'EM'){
  if(what == 'EM'){
    x <- fit$comp.time[2]
    if(attr(x, 'units') == 'mins') x <- x * 60
    rtn <- as.numeric(x)
  }else if(what == 'total'){
    x <- fit$comp.time[1]
    if(attr(x, 'units') == 'mins') x <- x * 60
    rtn <- as.numeric(x)
  }
  rtn
}

elapsed.JMb <- function(fit){
  if(!is.null(fit)) return(fit$time[3]) else return(NA)
}

elapsed.inj <- function(fit){
  
}

# Functions for object parsing --------------------------------------------



# Plotting elapsed times --------------------------------------------------
#' Gaussian
jML.times <- data.frame(method = 'joineRML',
                   `K1` = do.call(c, lapply(fit1.jML, elapsed.jML)),
                   `K2` = do.call(c, lapply(fit2.jML, elapsed.jML)),
                   `K3` = do.call(c, lapply(fit3.jML, elapsed.jML)))
aEM.times <- data.frame(method = 'Approximate EM',
                        `K1` = do.call(c, lapply(fit1, elapsed.aEM)),
                        `K2` = do.call(c, lapply(fit2, elapsed.aEM)),
                        `K3` = do.call(c, lapply(fit3, elapsed.aEM)))
JMb.times <- data.frame(method = 'JMbayes2',
                        `K1` = do.call(c, lapply(fit1.JMb, elapsed.JMb)),
                        `K2` = do.call(c, lapply(fit2.JMb, elapsed.JMb)),
                        `K3` = do.call(c, lapply(fit3.JMb, elapsed.JMb)))
inj.times <- data.frame(method = 'INLAjoint')

elapsed.times <- rbind(aEM.times, jML.times, JMb.times)

elapsed.times %>% 
  pivot_longer(K1:K3, names_to = 'K', values_to = 'elapsed') %>% 
  mutate_at('K', parse_number) %>% 
  mutate(
    s = ifelse(K == 1, '', 's'),
    K_label = paste0('K = ', K, ' Gaussian response', s)
  ) %>% 
  ggplot(aes(x = method, y = elapsed)) + 
  geom_boxplot() + 
  facet_wrap(~ K_label) +
  scale_y_log10() + 
  theme_light() + 
  labs(
    x = NULL, y = 'Elapsed time (seconds, log10)'
  )

#' Poisson
for(file in dir(dataDir, '^poisson')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}

aEM.times <- data.frame(method = 'Approximate EM',
                        `K1` = do.call(c, lapply(fit1, elapsed.aEM)),
                        `K2` = do.call(c, lapply(fit2, elapsed.aEM)),
                        `K3` = do.call(c, lapply(fit3, elapsed.aEM)))
JMb.times <- data.frame(method = 'JMbayes2',
                        `K1` = do.call(c, lapply(fit1.JMb, elapsed.JMb)),
                        `K2` = do.call(c, lapply(fit2.JMb, elapsed.JMb)),
                        `K3` = do.call(c, lapply(fit3.JMb, elapsed.JMb)))
inj.times <- data.frame(method = 'INLAjoint')

elapsed.times <- rbind(aEM.times, JMb.times)

elapsed.times %>% 
  pivot_longer(K1:K3, names_to = 'K', values_to = 'elapsed') %>% 
  mutate_at('K', parse_number) %>% 
  mutate(
    s = ifelse(K == 1, '', 's'),
    K_label = paste0('K = ', K, ' Poisson response', s)
  ) %>% 
  ggplot(aes(x = method, y = elapsed)) + 
  geom_boxplot() + 
  facet_wrap(~ K_label) +
  scale_y_log10() + 
  theme_light() + 
  labs(
    x = NULL, y = 'Elapsed time (seconds, log10)'
  )

# Tabulating estimates ----------------------------------------------------


