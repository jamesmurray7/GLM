# Extract out individual fits and store separately ------------------------
save.location <- '~/Downloads/'
vech <- function(X) X[lower.tri(X, diag = T)]
# Functions
extract.coeffs <- function(x){ # x a 'sub-list'
  xc <- x$coeffs
  out <- c(vech(xc$D), xc$beta,  xc$theta, xc$gamma, xc$eta)
  names(out)[1:3] <- c('D_11', 'D_21', 'D_22')
  names(out)[4:7] <- c('beta_(Intercept)', 'beta_time', 'beta_cont', 'beta_bin')
  out
}

extract.inits <- function(x){ # x a 'sub-list'
  xi <- x$inits
  D <- vech(xi$D); names(D) <-  c('D_11', 'D_21', 'D_22')
  c(D, xi$beta, xi$theta)
}

# Unpacking function
.unpack <- function(x){ # x a list of lists.
  for(i in 1:6){
    this <- x[[i]]
    ests <- do.call(rbind, lapply(this, function(x){
      if(!is.null(x)) extract.coeffs(x)
    }))
    inits <- do.call(rbind, lapply(this, function(x){
      if(!is.null(x)) extract.inits(x)
    }))
    out <- list(ests=ests, inits=inits)
    message('\nfits[[',i,']] had ', nrow(ests), ' successful fits out of 100.')
    message('\nSaving in ', save.location, 'nbests3-', i, '.RData')
    save(out, file = paste0(save.location, 'nbests3-', i, '.RData'))
  }
}

.unpack(fits)


# Now we can parse --------------------------------------------------------
source('./_Functions.R')
library(tidyverse)
#' ## Table of what set is what ##
#' ## beta = c(1, 0.1, 0.33, -0.5), 
#' ## D = matrix(c(0.5, 0, 0, 0.1), 2, 2), 
#' ##gamma = 0.5, surv.eta = c(0.05, -0.3)
#' ## ------------------------- ##
#' ## Fits1: from ntms = 10; 2: ntms = 15; 3: ntms = 10, gamma = -1

# Function to load one by one, rbind and store
loader <- function(i){
  load(paste0('~/Downloads/nbests3-',i,'.RData'))
  df <- as_tibble(out$ests) %>% pivot_longer(everything()) %>% mutate(a = as.character(i))
  message(i)
  df
}

df <- list()
for(i in 1:6){df[[i]] <- loader(i)}
df <- do.call(rbind, df)

df <- df %>% 
  mutate(
    description = case_when(
      a == '1' ~ 'theta = 0.25',
      a == '2' ~ 'theta = 0.50',
      a == '3' ~ 'theta = 0.75',
      a == '4' ~ 'theta = 1.00',
      a == '5' ~ 'theta = 1.50',
      a == '6' ~ 'theta = 2.00',
      T ~ 'AA'
    )
  )

targets <- data.frame(
  name = unique(df$name),
  target = c(0.5, 0, 0.1, 1, 0.1, 0.33, -0.50, NA, 0.5, 0.05, -0.3)
)

df %>% 
  left_join(., targets, 'name') %>% 
  ggplot(aes(x = value, colour = description)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~name, scales = 'free') + 
  theme_bw()

ggsave('~/Downloads/nb-ests.png')

# Thetas only
glimpse(df)
df %>% 
  filter(name == 'theta') %>% 
  separate(description, c('dummy', 'target'), sep = '\\=') %>% 
  mutate(dummy = paste0('var', dummy, '[', a, ']')) %>% 
  ggplot(aes(x = value)) + 
  geom_vline(aes(xintercept = as.numeric(target)), lty = 3, colour = 'steelblue') + 
  geom_density() + 
  facet_wrap(~dummy, scales = 'free', labeller = label_parsed) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 15)
  )
ggsave('~/Downloads/nb-ests-theta.png')

# Comparing starting values to  final -------------------------------------

loader2 <- function(i){
  load(paste0('~/Downloads/nbests3-',i,'.RData'))
  df <- as_tibble(out$inits) %>% pivot_longer(everything()) %>% mutate(a = as.character(i), b = 'Initial estimate')
  message(i)
  df
}

df.inits <- list()
for(i in 1:6) df.inits[[i]] <- loader2(i)
df.inits <- do.call(rbind, df.inits) %>% 
  as_tibble %>% 
  mutate(
    description = case_when(
      a == '1' ~ 'theta = 0.25',
      a == '2' ~ 'theta = 0.50',
      a == '3' ~ 'theta = 0.75',
      a == '4' ~ 'theta = 1.00',
      a == '5' ~ 'theta = 1.50',
      a == '6' ~ 'theta = 2.00',
      T ~ 'AA'
    )
  ) 


df$b <- 'aEM'
df.all <- rbind(df, df.inits)

df.all %>% 
  left_join(., targets, 'name') %>% 
  filter(grepl('beta', name) | grepl('alpha', name) | name %in% c('D_11', 'D_22')) %>% 
  ggplot(aes(x = value, colour = description, lty = b)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~name, scales = 'free') + 
  labs(lty = NULL, colour = 'Profile', x = 'Estimate') +
  theme_bw()
ggsave('~/Downloads/nb3-vs-inits.png')

# Thetas
df.all %>% 
  filter(grepl('theta', name) | name == 'V8') %>% 
  separate(description, c('dummy', 'target'), '\\s\\=\\s') %>% 
  mutate(
    target = as.numeric(target),
    dummy = paste0('var', dummy, '[', a, ']'),
    name = ifelse(name == 'V8', 'theta', name)
  ) %>% 
  ggplot(aes(x = value, colour = as.factor(b))) + 
  geom_vline(aes(xintercept = target), colour = 'steelblue', lty = 5) + 
  geom_density() + 
  facet_wrap(~dummy, scales = 'free', labeller = label_parsed) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 15)
  ) + 
  labs(colour = NULL)
ggsave('~/Downloads/nb3-vs-inits-theta.png')


# Timings -----------------------------------------------------------------
# reload the fits
load('~/Downloads/nbfits2.RData')

EM <- as.data.frame(sapply(1:6, function(i) do.call(c, lapply(fits[[i]], '[[', 4))))
tot <- as.data.frame(sapply(1:6, function(i) do.call(c, lapply(fits[[i]], '[[', 6))))
iter <- as.data.frame(sapply(1:6, function(i) do.call(c, lapply(fits[[i]], '[[', 5))))
all(EM < tot)

EM$a <- 'EM'; tot$a <- 'Total'; iter$a <- 'Iterations'

ets <- rbind(EM, tot, iter)
names(ets)[1:6] <- paste0('vartheta[', 1:6, ']')

ets %>% 
  as_tibble %>% 
  pivot_longer(`vartheta[1]`:`vartheta[6]`) %>% 
  arrange(name) %>% 
  ggplot(aes(x = fct_inorder(name), y = value)) + 
  geom_boxplot(outlier.alpha = .250) + 
  facet_wrap(~fct_inorder(a), scales = 'free') + 
  scale_x_discrete(labels = ggplot2:::parse_safe) + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) + 
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  ) + 
  labs(x = NULL, y = NULL)
ggsave('~/Downloads/ETs.png')

