#' ####
#' Parsing list of fits (from 18/11/21, gh = 3, random intercept in both count and ZI process)
#' ####
load('~/Downloads/fits.RData')
library(tidyverse)
# True values
beta <- c(1.5, 0.05, 0.33, 0.50)
alpha <- c(-0.5, 0.25)
D <- diag(c(.5^2, .15^2))

# functions
vech <- function(x) x[lower.tri(x, diag = T)]
extract.ests <- function(x){ # x a list
  X <- x$coeffs; Y <- x$inits
  beta <- c(X$beta); names(beta) <- names(Y$beta)
  alpha <- c(X$alpha); names(alpha) <- names(Y$alpha)
  vD <- diag(X$D); names(vD) <- c('D11', 'D22')
  yD <- diag(Y$D); names(yD) <- c('D11', 'D22')
  my <- c(vD, beta, alpha)
  inits <- c(yD, Y$beta, Y$alpha)
  list(my = my, inits = inits)
}

timings <- function(x) list(EMtime = x$EMtime, totaltime = x$totaltime, iterations = x$iter)


# Focus on just my ests first ---------------------------------------------
ests <- lapply(fits, extract.ests)
my.ests <- as.data.frame(do.call(rbind, lapply(ests, '[[', 1)))

# Targets
targets <- data.frame(
  name = names(my.ests),
  target = c(.5^2, .15^2, 1.5, 0.05, 0.33, 0.50, -.5, .25)
)

my.ests %>% 
  pivot_longer(D11:alpha_time) %>% 
  left_join(., targets, 'name') %>% 
  mutate(bias = value - target) %>% 
  ggplot(aes(x = value)) + 
  geom_vline(aes(xintercept = target), lty = 3) + 
  geom_density() + 
  facet_wrap(~name, scales='free') + 
  theme_light()
  
# Overlaying my estimates with inits - improvement? -----------------------
inits <- as.data.frame(do.call(rbind, lapply(ests, '[[', 2)))
inits <- inits %>% 
  pivot_longer(D11:alpha_time) %>% 
  mutate(a = 'Initial estimate')

my.ests %>% 
  pivot_longer(D11:alpha_time) %>% 
  mutate(a = 'approx. EM') %>% 
  rbind(., inits) %>% 
  left_join(., targets, 'name') %>% 
  ggplot(aes(x = value, colour = as.factor(a))) + 
  geom_vline(aes(xintercept = target), lty = 3) + 
  geom_density() + 
  facet_wrap(~name, scales='free') + 
  theme_light() + 
  labs(colour = NULL)

  
ggsave('~/Downloads/EMvsglmmTMB.png')
