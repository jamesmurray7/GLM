setwd("~/Documents/GLMM/mix/_01")
load('./pbc-jmb2.RData')
load('myfit2.RData')
load('~/Documents/GLMM/mix/_01/pbc.RData')
source('EM.R')

jmb.fit$running_time[3]/60 # > 1 hour
my.fit3$totaltime # < 20s minutes


# Extraction/compilation functions ----------------------------------------

# Extraction from JMbayes2 fit
jmb.extract <- function(fit){
  sj <- summary(jmb.fit)
  G <- sj$Outcome1
  B <- sj$Outcome2
  P <- sj$Outcome3
  surv <- sj$Surv
  
  # Setting names to match my fit for ease of comparison.
  beta_rhs <- c('_(Intercept)', '_time', '_cont', '_bin')
  # 1. Gaussian
  G_beta <- G[1:4,1]; names(G_beta) <- paste0('G', beta_rhs)
  var.e <- G[5,1]^2; names(var.e) <- 'var.e'
  # CI
  G_beta95 <- G[1:4, 3:4]; row.names(G_beta95) <- names(G_beta)
  var.e95 <- G[5, 3:4]^2; row.names(var.e95) <- 'var.e'
  # 2. Binomial
  B_beta <- B[1:4, 1]; names(B_beta) <- paste0('B', beta_rhs)
  # CI
  B_beta95 <- B[1:4, 3:4];  row.names(B_beta95) <- names(B_beta)
  # 3. Poisson
  P_beta <- P[1:4, 1]; names(P_beta) <- paste0('P', beta_rhs)
  # CI
  P_beta95 <- P[1:4, 3:4];  row.names(P_beta95) <- names(P_beta)
  # 4. Survival
  eta <- surv[1:2, 1]; names(eta) <- c('surv_cont', 'surv_bin')
  gamma <- surv[3:5, 1]; names(gamma) <- paste0('gamma_', 1:3)
  # CI
  eta95 <- surv[1:2, 3:4]; row.names(eta95) <- names(eta)
  gamma95 <- surv[3:5, 3:4]; row.names(gamma95) <- names(gamma)
  # 5. Covariance D
  D <- vech(sj$D)
  names(D) <- paste0('[', apply(which(lower.tri(sj$D, T), arr.ind = T), 1, paste0, collapse = ', '), ']')
  
  # Pack and return
  estimates <- list(
    D,
    c(G_beta, B_beta, P_beta),var.e,
    gamma, eta
  )
  CIs <- list(
    rbind(G_beta95, B_beta95, P_beta95), var.e95,
    gamma95, eta95
  )
  return(list(estimates, CIs))
}

# Extraction from my AEM fit
my.extract <- function(x){
  xc <- x$coeffs
  vD <- vech(xc$D)
  names(vD) <- paste0('[' , apply(which(lower.tri(xc$D, T), arr.ind = T), 1, paste0, collapse = ', '), ']')
  beta <- xc$beta
  names(beta) <- paste0(rep(c('G_', 'B_', 'P_'), each = 4), c('(Intercept)', 'time', 'cont', 'bin'))
  out <- c(vD, beta, 'var.e' = xc$var.e, xc$gamma, xc$eta); names(out)[38:39] <- c('surv_cont', 'surv_bin')
  out
}



# Point-wise comparison ---------------------------------------------------
library(tidyverse)
library(xtable)

pointwise <- left_join(
  my.extract(my.fit3) %>% t %>% as_tibble() %>% pivot_longer(everything(), values_to = 'Approximate EM'),
  do.call(c, el(jmb.extract(jmb.fit))) %>% t %>% as_tibble %>% pivot_longer(everything(), values_to = 'JMbayes2'),
  'name'
)

# Comparison via 95% CIs from JMbayes2 ------------------------------------

cis <- do.call(rbind, el(jmb.extract(jmb.fit), 2)) %>% rownames_to_column('name') %>% 
  left_join(pointwise, ., 'name') %>% 
  mutate(flag = ifelse(`Approximate EM` < `97.5%` & `Approximate EM` > `2.5%`, 1 ,0))

# Plot?
cis %>% 
  filter(!grepl('^\\[', name)) %>% 
  ggplot(aes(x = fct_inorder(name), y = JMbayes2)) + 
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`), size = .25) + 
  geom_point(aes(y = `Approximate EM`), colour = 'red') + 
  coord_flip() + 
  theme_light() + 
  labs(x = '', y = '')

# ggsave('~/Downloads/jm-vs-me2.png')

# Tabulation
is_in <- function(x) paste0('textcolour{darkgreen}{', as.character(x), '}')
is_not_in <- function(x) paste0('textcolour{red}{', as.character(x), '}')

table <- cis %>% 
  filter(!grepl('^\\[', name)) %>% 
  mutate_at(vars(`Approximate EM`: `97.5%`), ~ format(round(.x, 3), nsmall = 3)) %>% 
  mutate(out = ifelse(flag == 1, is_in(`Approximate EM`), is_not_in(`Approximate EM`)),
         ci = paste0('[', `2.5%`, ', ', `97.5%`, ']')) %>% 
  select(name, out, JMbayes2, ci)
xtable(table)
