# Comparing intercept and slope fits across fitting methods visually.
rm(list=ls())
library(tidyverse)
source('EM.R')
source('Simulations/funs.R')
theme_set(theme_bw())
dataDir <- paste0(getwd(), '/Simulations/fits')

# Plotting elapsed times --------------------------------------------------
#' Gaussian
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^gaussian')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}

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
  ) -> Gplot

#' Poisson
rm(list = ls()[grepl('fit',ls())])
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
  ) -> Pplot

#' Binomial
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^binomial')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}

# Plotting Estimates ------------------------------------------------------

#' Gaussian
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^gaussian')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}

# K=1
plot.ests(fit1, fit1.jML, fit1.JMb, fit1.INLA, f = 'gaussian', K = 1)
plot.ests(fit2, fit2.jML, fit2.JMb, fit2.INLA, f = 'gaussian', K = 2)
plot.ests(fit3, fit3.jML, fit3.JMb, fit3.INLA, f = 'gaussian', K = 3)

#' Poisson
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^poisson')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}
plot.ests(fit1, NULL, fit1.JMb, fit1.INLA, f = 'poisson', K = 1)
plot.ests(fit2, NULL, fit2.JMb, fit2.INLA, f = 'poisson', K = 2)
plot.ests(fit3, NULL, fit3.JMb, fit3.INLA, f = 'poisson', K = 3)

#' Binomial
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^binom')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  load(paste0(dataDir, '/', file))
}
plot.ests(fit1, NULL, NULL, fit1.INLA, f = 'binomial', K = 1)
plot.ests(fit2, NULL, NULL, fit2.INLA, f = 'binomial', K = 2)
plot.ests(fit3, NULL, NULL, NULL, f = 'binomial', K = 3)

