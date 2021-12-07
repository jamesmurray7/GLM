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
    message('\nSaving in ', save.location, 'nbests2-', i, '.RData')
    save(out, file = paste0(save.location, 'nbests2-', i, '.RData'))
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
  load(paste0('~/Downloads/nbests2-',i,'.RData'))
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

# with neg
targets2 <- data.frame(
  name = unique(df$name),
  target = c(0.5, 0, 0.1, 1, 0.1, 0.33, -0.50, 1.5, -1, 0.05, -0.3)
)
df %>% 
  filter(a %in% c('3')) %>% 
  left_join(., targets2, 'name') %>% 
  ggplot(aes(x = value, colour = description)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~name, scales = 'free') + 
  theme_bw()


# Comparing starting values to  final -------------------------------------

loader2 <- function(i){
  load(paste0('~/Downloads/nbests-',i,'.RData'))
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
      a == '1' ~ 'ntms = 10',
      a == '2' ~ 'ntms = 15',
      a == '3' ~ 'ntms = 10, gamma = -1',
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
ggsave('~/Downloads/nb-vs-inits.png')


# Timings -----------------------------------------------------------------

aa <- function(x) list(EM = x$EMtime, total = x$totaltime)
ntms10 <- lapply(fits[[1]], aa)
ntms15 <- lapply(fits[[2]], aa)

EM <- data.frame(ntms10 = do.call(c, lapply(ntms10, '[[', 1)),
                 ntms15 = do.call(c, lapply(ntms15, '[[', 1)))
tot <- data.frame(ntms10 = do.call(c, lapply(ntms10, '[[', 2)),
                  ntms15 = do.call(c, lapply(ntms15, '[[', 2)))
par(mfrow = c(1,2))
boxplot(EM)
boxplot(tot)
dev.off()
