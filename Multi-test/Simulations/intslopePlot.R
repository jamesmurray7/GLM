# Comparing intercept and slope fits across fitting methods visually.
rm(list=ls())
library(tidyverse)
source('EM.R')
source('Simulations/funs.R')
theme_set(theme_bw())
dataDir <- paste0(getwd(), '/Simulations/fits')

# Function to create med[iqr] to compare
make.et.summary <- function(et){
  et %>%
    pivot_longer(K1:K3) %>% 
    group_by(method, name) %>%
    summarise(avg_time = median(value, na.rm = T),
              p25 = quantile(value, probs = .25, na.rm = T),
              p75 = quantile(value, probs = .75, na.rm = T),
              .groups = 'keep') %>% ungroup %>%
    mutate(across(avg_time:p75, ~ format(round(.x, 2), nsmall = 2))) %>%
    mutate(out = paste0(avg_time, ' [', p25, ', ', p75, ']')) %>%
    select(method, name, out) %>%
    pivot_wider(id_cols = method, names_from = name, values_from = out) %>%
    as.data.frame() %>%
    xtable %>%
    print(., include.rownames = T)
}


# Plotting elapsed times --------------------------------------------------
#' Gaussian
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^gaussian')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  # Ensure we take second run at JMbayes2
  if(grepl('JMbayes2', file)){
    if(grepl('2\\-2', file)){
      load(paste0(dataDir, '/', file))
    }else{
      cat('--> Not loading first pass of JMbayes\n')
    }
  }else{
    load(paste0(dataDir, '/', file))
  }
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

make.et.summary(elapsed.times)

elapsed.times %>% 
  pivot_longer(K1:K3, names_to = 'K', values_to = 'elapsed') %>% 
  mutate_at('elapsed', log10) %>% # Different to scale_y_log10?
  mutate_at('K', parse_number) %>% 
  mutate(
    s = ifelse(K == 1, '', 's'),
    K_label = paste0('K = ', K)#, ' Gaussian response', s)
  ) %>% 
  filter(!(method == 'joinerRML' & elapsed > 250)) %>%
  ggplot(aes(x = K_label, y = elapsed)) + # log10(elapsed))) + 
  geom_boxplot(outlier.alpha = .50) + 
  stat_summary(fun=median, geom="line", aes(group=1), lty = 5, alpha = .5) +
  stat_summary(fun=median, geom="point", size = 2) + 
  facet_wrap(~ method, scales = 'free_y') +
  theme_csda() + 
  theme(
    strip.text = element_text(vjust = 1.25) 
  ) +
  labs(
    x = NULL, y = expression(underline(bold(Gaussian)))
  ) -> Gplot

#' Poisson
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^poisson')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  # Ensure we take second run at JMbayes2
  if(grepl('JMbayes2', file)){
    if(grepl('2\\-2', file)){
      load(paste0(dataDir, '/', file))
    }else{
      cat('--> Not loading first pass of JMbayes\n')
    }
  }else{
    load(paste0(dataDir, '/', file))
  }
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
make.et.summary(elapsed.times)

elapsed.times %>% 
  pivot_longer(K1:K3, names_to = 'K', values_to = 'elapsed') %>% 
  mutate_at('K', parse_number) %>% 
  mutate(
    s = ifelse(K == 1, '', 's'),
    K_label = paste0('K = ', K)#, ' count response', s)
  ) %>% 
  ggplot(aes(x = K_label, y = log10(elapsed))) + 
  geom_boxplot() + 
  stat_summary(fun=median, geom="line", aes(group=1), lty = 5, alpha = .5) +
  stat_summary(fun=median, geom="point", size = 2) + 
  facet_wrap(~ method, scales = 'free_y') +
  # scale_y_log10() +
  labs(
    x = NULL, y = expression(underline(bold(Poisson)))
  ) + 
  theme_csda() + 
  theme(
    strip.text = element_text(vjust = 1.25) 
  ) -> Pplot

#' Binomial
rm(list = ls()[grepl('fit',ls())])
for(file in dir(dataDir, '^binomial')){ # load in all data -> Gaussian
  cat(paste0(file, '\n'))
  # Ensure we take second run at JMbayes2
  if(grepl('JMbayes2', file)){
    if(grepl('2\\-2', file)){
      load(paste0(dataDir, '/', file))
    }else{
      cat('--> Not loading first pass of JMbayes\n')
    }
  }else{
    load(paste0(dataDir, '/', file))
  }
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
make.et.summary(elapsed.times)
elapsed.times %>% 
  pivot_longer(K1:K3, names_to = 'K', values_to = 'elapsed') %>% 
  mutate_at('K', parse_number) %>% 
  mutate(
    s = ifelse(K == 1, '', 's'),
    K_label = paste0('K = ', K)#, ' count response', s)
  ) %>% 
  ggplot(aes(x = K_label, y = log10(elapsed))) + 
  geom_boxplot() + 
  stat_summary(fun=median, geom="line", aes(group=1), lty = 5, alpha = .5) +
  stat_summary(fun=median, geom="point", size = 2) + 
  facet_wrap(~ method, scales = 'free_y') +
  labs(
    x = NULL, y = expression(underline(bold(Binomial)))
  ) + 
  theme_csda() + 
  theme(
    strip.text = element_text(vjust = 1.25) 
  )-> Bplot

#' Combine?
all.plots <- ggpubr::ggarrange(Gplot,Pplot,Bplot,ncol=1)
ggpubr::annotate_figure(all.plots, left = 'Elapsed time (seconds, log10)')
ggsave('~/Downloads/all_elapsed.png',
       height = 240, width = 190, units = 'mm')
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
p <- plot.ests(fit3, fit3.jML, fit3.JMb, fit3.INLA, f = 'gaussian', K = 3)

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
plot.ests(fit1, NULL, fit1.JMb, fit1.INLA, f = 'binomial', K = 1)
plot.ests(fit2, NULL, fit2.JMb, fit2.INLA, f = 'binomial', K = 2)
plot.ests(fit3, NULL, fit3.JMb, fit3.INLA, f = 'binomial', K = 3)

