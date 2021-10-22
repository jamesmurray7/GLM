setwd('~/Documents/GLMM')

files.to.load <- as.list(dir('./Simulations/Poisson/', pattern = 'RData'))

in.dir <- paste0(getwd(), '/Simulations/Poisson/')

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
  c(D, beta, gamma, eta, EM, tot)
}

loadfn <- function(f){
  message(f)
  load(paste0(in.dir, f))
  ests <- as.data.frame(do.call(rbind, lapply(fits, extract.estimates)))
  aa <- stringr::str_extract(f, 'sim\\d')
  if(aa == 'sim1') rate <- "lower" else rate <- "higher"
  bb <- stringr::str_extract(f, '\\d')
  if(bb == '3') gh <- 3 else gh <- 9
  
  ests$rate <- rate
  ests$gh <- gh
  return(ests)
}

ests <- as.data.frame(do.call(rbind, lapply(files.to.load, loadfn)))

targets <- data.frame(
  parameter = names(ests[-c(18:21)]),
  targets = c(.5^2, .2^2, .5^2, .2^2, -0.125,
              0.75, -0.05, -0.2, 0.2,
              1, 0.05,  0.5, -0.80,
              -0.25, 0.6, -0.1, 0.25)
)
  

ests %>% 
  pivot_longer(cols = `D1,1`:`Elapsed time`,
               names_to = 'parameter', values_to = 'estimate') %>% 
  left_join(., targets, 'parameter') %>% 
  ggplot(aes(x = estimate, lty = as.factor(gh), col = rate)) + 
  geom_vline(aes(xintercept = targets), lty = 5, colour = 'grey50') + 
  geom_density() +
  scale_colour_manual(values = c('black', 'red')) + 
  scale_linetype_manual(values = c(1, 3)) + 
  facet_wrap(~parameter, scales = 'free')
