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
    if(attr(x, 'units') == 'mins') x <- x * 60 else x <- x
    rtn <- as.numeric(x)
  }else if(what == 'total'){
    x <- fit$comp.time[1]
    if(attr(x, 'units') == 'mins') x <- x * 60
    rtn <- as.numeric(x)
  }
  rtn
}

elapsed.JMb <- function(fit){
  if(!is.null(fit)) return(fit$comp.time[3]) else return(NA)
}

elapsed.INLA <- function(fit){
  if(!is.null(fit)) return(unname(fit$comp.time[2])) else return(NA)
}

# Functions for object parsing --------------------------------------------

parse.aEM <- function(fit){
  co <- fit$coeffs
  beta <- co$beta
  if(any(unlist(co$sigma) != 0)) sig <- unlist(co$sigma) else sig <- NULL
  gam <- co$gamma
  ze <- co$zeta
  SE <- fit$SE[!grepl('^D', names(fit$SE))]
  Omega <- setNames(c(beta, sig, gam, ze), names(SE))
  
  lb <- Omega - qnorm(.975) * SE; ub <- Omega + qnorm(.975) * SE
  return(cbind(Omega, SE, lb, ub))
}

parse.JMb <- function(fit){
  o1 <- fit$Outcome1
  if(!is.null(fit$Outcome2)) o2 <- fit$Outcome2 else o2 <- NULL 
  if(!is.null(fit$Outcome3)) o3 <- fit$Outcome3 else o3 <- NULL
  o <- rbind(o1, o2, o3)
  s <- rbind(fit$survival[-1, ], fit$survival[1,]) # swap
  as.matrix(rbind(o, s)[,c('Mean', 'StDev', '2.5%', '97.5%')])
}

parse.INLA <- function(fit){
  f <- fit$fixed
  z <- fit$survival[-c(1,2), ]
  g <- fit$gamma
  method <- 'INLA'
  cbind(rbind(f, g, z), method)[,c('mean', 'sd', '0.025quant', '0.975quant', 'method')]
}

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
inj.times <- data.frame(method = 'INLAjoint',
                        `K1` = do.call(c, lapply(fit1.INLA, elapsed.INLA)),
                        `K2` = do.call(c, lapply(fit2.INLA, elapsed.INLA)),
                        `K3` = do.call(c, lapply(fit3.INLA, elapsed.INLA)))

elapsed.times <- rbind(aEM.times, jML.times, JMb.times, inj.times)

elapsed.times %>% 
  filter(method!='JMbayes2') %>% 
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
  ) ->gfits

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
inj.times <- data.frame(method = 'INLAjoint',
                        `K1` = do.call(c, lapply(fit1.INLA, elapsed.INLA)),
                        `K2` = do.call(c, lapply(fit2.INLA, elapsed.INLA)),
                        `K3` = do.call(c, lapply(fit3.INLA, elapsed.INLA)))

elapsed.times <- rbind(aEM.times, JMb.times, inj.times)

elapsed.times %>% 
  pivot_longer(K1:K3, names_to = 'K', values_to = 'elapsed') %>% 
  mutate_at('K', parse_number) %>% 
  mutate(
    s = ifelse(K == 1, '', 's'),
    K_label = paste0('K = ', K, ' count response', s)
  ) %>% 
  ggplot(aes(x = method, y = elapsed)) + 
  geom_boxplot() + 
  facet_wrap(~ K_label) +
  scale_y_log10() + 
  theme_light() + 
  labs(
    x = NULL, y = 'Elapsed time (seconds, log10)'
  ) -> counts

# Plotting Estimates ------------------------------------------------------

#' K = 1 -----
nm <- c(paste0('beta[', 0:3, ']'), 'gamma[1]', 'zeta') # setting names
mine <- as.data.frame(do.call(rbind,lapply(fit1, parse.aEM)))
JMb <- as.data.frame(do.call(rbind, lapply(fit1.JMb, parse.JMb)))
inla <- do.call(rbind, lapply(fit1.INLA, parse.INLA))

#' Mine
mine$method <- 'Approximate EM'
mine$param <- rep(nm, 100)  
mine <- mine %>% select(method, param, est = Omega)

#' JMbayes2
JMb$method <- 'JMbayes2'
JMb$param <- rep(nm, 100)
JMb <- JMb %>% select(method, param, est = Mean)

#' INLA
inla$param <- rep(nm, 100) 
inla <- inla %>% select(method, param, est = mean)

all <- rbind(mine, JMb, inla)

ggplot(all, aes(x = est, colour = method)) + 
  geom_density() + 
  facet_wrap(~param, labeller = label_parsed, scales = 'free')


#' K = 2 ----
nm <- c(paste0(rep(c('beta[1', 'beta[2'), each = 4), 0:3, ']'), 
        'gamma[1]', 'gamma[2]', 'zeta')
mine <- as.data.frame(do.call(rbind,lapply(fit2, parse.aEM)))
JMb <- as.data.frame(do.call(rbind, lapply(fit2.JMb, parse.JMb)))
inla <- do.call(rbind, lapply(fit2.INLA, parse.INLA))

#' Mine
mine$method <- 'Approximate EM'
mine$param <- rep(nm, 100)  
mine <- mine %>% select(method, param, est = Omega)

#' JMbayes2
JMb$method <- 'JMbayes2'
JMb$param <- rep(nm, 100)
JMb <- JMb %>% select(method, param, est = Mean)

#' INLA
inla$param <- rep(nm, 100) 
inla <- inla %>% select(method, param, est = mean)

all <- rbind(mine, JMb, inla)

ggplot(all, aes(x = est, colour = method)) + 
  geom_density() + 
  facet_wrap(~param, labeller = label_parsed, scales = 'free')

#' K = 3 ----
nm <- c(paste0(rep(c('beta[1', 'beta[2', 'beta[3'), each = 4), 0:3, ']'), 
        'gamma[1]', 'gamma[2]', 'gamma[3]', 'zeta')
mine <- as.data.frame(do.call(rbind,lapply(fit3, parse.aEM)))
JMb <- as.data.frame(do.call(rbind, lapply(fit3.JMb, parse.JMb)))
inla <- do.call(rbind, lapply(fit3.INLA, parse.INLA))

#' Mine
mine$method <- 'Approximate EM'
mine$param <- rep(nm, 100)  
mine <- mine %>% select(method, param, est = Omega)

#' JMbayes2
JMb$method <- 'JMbayes2'
JMb$param <- rep(nm, 100)
JMb <- JMb %>% select(method, param, est = Mean)

#' INLA
inla$param <- rep(nm, 100) 
inla <- inla %>% select(method, param, est = mean)

all <- rbind(mine, JMb, inla)

ggplot(all, aes(x = est, colour = method)) + 
  geom_density() + 
  facet_wrap(~param, labeller = label_parsed, scales = 'free')