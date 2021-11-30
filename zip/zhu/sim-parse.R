# Extract out individual fits and store separately ------------------------
save.location <- '~/Downloads/'
vech <- function(X) X[lower.tri(X, diag = T)]
# Functions
extract.coeffs <- function(x){ # x a 'sub-list'
  xc <- x$coeffs
  out <- c(vech(xc$D), xc$beta, xc$alpha, xc$gamma, xc$eta)
  names(out)[1:3] <- c('D11', 'D21', 'D22')
  out
}

extract.inits <- function(x){ # x a 'sub-list'
  xi <- x$inits
  D <- vech(xi$D); names(D) <-  c('D11', 'D21', 'D22')
  c(D, xi$beta, xi$alpha)
}
# Unpacking function
.unpack <- function(x){ # x a list of lists.
  for(i in 1:3){
    this <- x[[i]]
    ests <- do.call(rbind, lapply(this, function(x){
      if(!is.null(x)) extract.coeffs(x)
    }))
    inits <- do.call(rbind, lapply(this, function(x){
      if(!is.null(x)) extract.inits(x)
    }))
    out <- list(ests=ests, inits=inits)
    message('\nfits[[',i,']] had ', nrow(ests), ' successful fits out of 100.')
    message('\nSaving in ', save.location, 'ests291121-new-', i, '.RData')
    save(out, file = paste0(save.location, 'ests291121-new-', i, '.RData'))
  }
}

.unpack(fits)


# Now we can parse --------------------------------------------------------
source('./_Functions.R')
library(tidyverse)
#' ## Table of what set is what ##
#' ## D = zhuD(); beta = zhubeta(); alpha = zhualpha();
#' ## gamma = (-0.6, 0.4), eta = -1.50
#' ## 1: ntms = 6; 2: ntms = 10, 3: ntms = 14.
#' ## ------------------------- ##

# Function to load one by one, rbind and store
loader <- function(i){
  load(paste0('~/Downloads/ests291121-new-',i,'.RData'))
  df <- as_tibble(out$ests) %>% pivot_longer(everything()) %>% mutate(a = as.character(i))
  message(i)
  df
}

df <- list()
for(i in 1:3){df[[i]] <- loader(i)}
df <- do.call(rbind, df)

df <- df %>% 
  mutate(
    description = case_when(
      a == '1' ~ 'ntms = 6',
      a == '2' ~ 'ntms = 10',
      a == '3' ~ 'ntms = 14',
      T ~ 'AA'
    )
  )

targets <- data.frame(
  name = unique(df$name),
  target = c(vech(zhuD()), zhubeta(), zhualpha(), -0.6, 0.4, -1.5)
)

df %>% 
  left_join(., targets, 'name') %>% 
  ggplot(aes(x = value, colour = description)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~name, scales = 'free') + 
  theme_bw()

ggsave('~/Downloads/zhu-like-new.png')

# Comparing starting values to  final -------------------------------------

loader2 <- function(i){
  load(paste0('~/Downloads/ests291121-new-',i,'.RData'))
  df <- as_tibble(out$inits) %>% pivot_longer(everything()) %>% mutate(a = as.character(i), b = 'Initial estimate')
  message(i)
  df
}

df.inits <- list()
for(i in 1:3) df.inits[[i]] <- loader2(i)
df.inits <- do.call(rbind, df.inits) %>% 
  as_tibble %>% 
  mutate(
    description = case_when(
      a == '1' ~ 'ntms = 6',
      a == '2' ~ 'ntms = 10',
      a == '3' ~ 'ntms = 14',
      T ~ 'AA'
    )
  ) 

df$b <- 'aEM'
df.all <- rbind(df, df.inits)

df.all %>% 
  left_join(., targets, 'name') %>% 
  filter(grepl('beta', name) | grepl('alpha', name) | name %in% c('D11', 'D22')) %>% 
  ggplot(aes(x = value, colour = description, lty = b)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~name, scales = 'free') + 
  labs(lty = NULL, colour = 'Profile', x = 'Estimate') +
  theme_bw()
