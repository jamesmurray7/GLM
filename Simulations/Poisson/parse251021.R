#' ####
#' Parsing simulations
#' ####
library(tidyverse)

setwd('~/Documents/GLMM')

files.to.load <- as.list(dir('./Poisson', pattern = 'sim\\d\\.RData'))

in.dir <- paste0(getwd(), '/Poisson/')

extract.estimates <- function(fit){
  # D
  D <- diag(fit$coeffs$D)
  D.31 <- fit$coeffs$D[3, 1]
  D <- c(D, D.31); names(D) <- c('D1,1', 'D2,2', 'D3,3', 'D4,4', 'D3,1')
  # \beta
  beta <- fit$coeffs$beta; names(beta) <- paste0(rep(c('beta1_', 'beta2_'), each  = 4), gsub('\\(|\\)', '', names(beta)))
  # \eta
  eta <- fit$coeffs$eta; names(eta) <- paste0('surv_', names(eta))
  # \gamma
  gamma <- fit$coeffs$gamma; names(gamma) <- c('gamma_1', 'gamma_2')
  # Timings
  EM <- fit$EMtime; names(EM) <- 'EM time'
  tot <- fit$comp.time; names(tot) <- 'Elapsed time'
  iter <- max(fit$history$iter); names(iter) <- 'Iterations'
  c(D, beta, gamma, eta, EM, tot, iter)
}

#' sim1: 30% failure rate on 15ntms; 
#' sim2: 40% failure rate on 10ntms.
#' sim3: 50% failure rate on 6ntms

loadfn <- function(f){
  message(f)
  load(paste0(in.dir, f))
  if(grepl('fits3', f)){
    gh <- 3
    ests <- as.data.frame(do.call(rbind, lapply(fits3, extract.estimates)))
  }else{
    gh <- 9
    ests <- as.data.frame(do.call(rbind, lapply(fits9, extract.estimates)))
  }
  
  if(grepl('sim1', f)){
    prof <- 'long'
  }else if(grepl('sim2', f)){
    prof <- 'medium'
  }else if(grepl('sim3', f)){
    prof <- 'short'
  }else{
    stop('something wrong...', f)
  }
  
  ests$gh <- gh
  ests$prof <- prof
  return(ests)
}

ests <- as.data.frame(do.call(rbind, lapply(files.to.load, loadfn)))

targets <- data.frame(
  parameter = names(ests[-c(18:22)]),
  targets = c(.5^2, .2^2, .5^2, .2^2, -0.125,
              0.75, -0.05, -0.2, 0.2,
              1, 0.05,  0.5, -0.80,
              -0.5, 1, 0.05, -0.30)
)


ests %>% 
  pivot_longer(cols = `D1,1`:`surv_bin`,
               names_to = 'parameter', values_to = 'estimate') %>% 
  left_join(., targets, 'parameter') %>% 
  ggplot(aes(x = estimate, colour = prof, linetype = as.factor(gh))) +
  geom_vline(aes(xintercept = targets), lty = 5, colour = 'grey50') + 
  geom_density() +
  scale_colour_manual(values = c('black', 'red', 'blue')) + 
  scale_linetype_manual(values = c(1, 3)) + 
  facet_wrap(~parameter, scales = 'free') + 
  labs(y = NULL, x = 'Estimate', lty = '# GH nodes', colour = 'Profile length') + 
  theme_csda()

# Elapsed time and iterations
ests %>% 
  select(gh, prof, `EM time`, `Elapsed time`, `Iterations`) %>% 
  pivot_longer(cols = `EM time`:`Iterations`) %>% 
  ggplot(aes(x = name, y = value, fill = as.factor(gh)))+
  geom_boxplot(outlier.alpha = .5) + 
  facet_wrap(~prof, scales = 'free') + 
  theme_csda() + 
  labs(x=NULL, y = NULL)
ggsave('~/Downloads/K2-Poisson-Elapsed.png')
