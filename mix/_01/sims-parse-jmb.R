library(tidyverse)
vech <- function(x) x[lower.tri(x, T)]
load('~/Downloads/jmbfits.RData')
load('~/Downloads/mixests-1.RData')
out <- out %>% as.data.frame %>%
  rename(surv_bin=bin, surv_cont=cont) %>% 
  pivot_longer(everything()) %>% 
  mutate(a = 'Approximate EM')

fits <- do.call(rbind, lapply(lapply(jmb.fits,'[[', 1), function(s) do.call(c, s))) %>% 
  as.data.frame %>% 
  pivot_longer(everything()) %>% 
  mutate(a = 'JMbayes2') %>% 
  rbind(out, .)

targets <- data.frame(
  name = unique(fits$name),
  target = c(vech(true.D), 
             do.call(c, lapply(1:3, function(i) true.beta[i,])),
             0.25, true.gamma, true.eta)
)


fits %>% 
  left_join(., targets, 'name') %>% 
  filter(grepl('^\\[', name)) %>% 
  ggplot(aes(x = value, colour = a)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~fct_inorder(name), scales = 'free', nc = 4) + 
  theme_bw()
ggsave('~/Downloads/jm-vs-me-ntms10simsD.png')
