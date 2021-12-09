# Extract out individual fits and store separately ------------------------
save.location <- '~/Downloads/'
vech <- function(X) X[lower.tri(X, diag = T)]
# Functions
extract.coeffs <- function(x){ # x a 'sub-list'
  xc <- x$coeffs
  apply(which(lower.tri(xc$D, T), arr.ind = T), 1, paste0, collapse = ', ')
  vD <- vech(xc$D)
  names(vD) <- paste0('[' , apply(which(lower.tri(xc$D, T), arr.ind = T), 1, paste0, collapse = ', '), ']')
  beta <- xc$beta
  names(beta) <- paste0(rep(c('G_', 'B_', 'P_'), each = 4), c('(Intercept)', 'time', 'cont', 'bin'))
  out <- c(vD, beta, xc$gamma, xc$eta)
  out
}

# Currently not keeping inits due to large object size...
# extract.inits <- function(x){ # x a 'sub-list'
#   xi <- x$inits
#   D <- vech(xi$D); names(D) <-  c('D_11', 'D_21', 'D_22')
#   c(D, xi$beta, xi$theta)
# }

ests <- as.data.frame(do.call(rbind, lapply(fit, extract.coeffs)))


# Now we can parse --------------------------------------------------------
source('./simData.R')
library(tidyverse)

df <- as_tibble(ests) %>% pivot_longer(everything())

targets <- data.frame(
  name = unique(df$name),
  target = c(vech(true.D), do.call(c, lapply(1:3, function(i) true.beta[i,])), true.gamma, true.eta)
)

df %>% 
  filter(!str_detect(name, '^\\[')) %>% 
  left_join(., targets, 'name') %>% 
  ggplot(aes(x = value)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~name, scales = 'free') + 
  theme_bw()

ggsave('~/Downloads/mix-ests-noD.png')

df %>% 
  filter(!str_detect(name, '^\\[', negate = T)) %>% 
  left_join(., targets, 'name') %>% 
  ggplot(aes(x = value)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~fct_inorder(name), scales = 'free') + 
  theme_bw()

ggsave('~/Downloads/mix-ests-D.png')

# Timings -----------------------------------------------------------------

aa <- function(x) list(EM = x$EMtime, total = x$totaltime)
ets <- lapply(fit, aa)

EM <- tibble(elapsed = do.call(c, lapply(ets, '[[', 1)), name = 'EM only')
tot <- tibble(elapsed = do.call(c, lapply(ets, '[[', 2)), name = 'Total time')

ets <- rbind(EM, tot)

ets %>% 
  ggplot(aes(x = name, y= elapsed)) + 
  geom_boxplot()
ggsave('~/Downloads/mix-ets.png')
