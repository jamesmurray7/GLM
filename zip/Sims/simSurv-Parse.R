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
  for(i in 1:4){
    this <- x[[i]]
    ests <- do.call(rbind, lapply(this, function(x){
      if(!is.null(x)) extract.coeffs(x)
    }))
    inits <- do.call(rbind, lapply(this, function(x){
      if(!is.null(x)) extract.inits(x)
    }))
    out <- list(ests=ests, inits=inits)
    message('\nfits[[',i,']] had ', nrow(ests), ' successful fits out of 100.')
    message('\nSaving in ', save.location, 'ests2-', i, '.RData')
    save(out, file = paste0(save.location,'ests2-',i,'.RData'))
  }
}

.unpack(fits)


# Now we can parse --------------------------------------------------------
library(tidyverse)
#' ## Table of what set is what ##
#' ## 1: n = 500, ntms = 5,  gamma = -0.50, surv.eta = c(0.05, -0.3), theta = c(-4, .2))
#' ## 2: n = 500, ntms = 10, gamma = -0.50, surv.eta = c(0.05, -0.3), theta = c(-4, .2))
#' ## 3: n = 500, ntms = 15, gamma = -0.50, surv.eta = c(0.05, -0.3), theta = c(-4, .2))
#' ## 4: n = 500, ntms = 15, gamma = -0.50, surv.eta = c(0.05, -0.3), theta = c(-6, .2))
#' ## 5: n = 500, ntms = 15, gamma = -0.50, surv.eta = c(0.05, -0.3), theta = c(-4, .2)) (same as 3, but D lower variance terms)
#' ## 6: n = 500, ntms = 15, gamma = -0.50, surv.eta = c(0.05, -0.3), theta = c(-6, .2)) (same as 4, but D lower variance terms)
#' ## 7: n = 500, ntms = 15, gamma =  0.50, surv.eta = c(0.05, -0.3), theta = c(-6, .2)) (same as 6, but positive gamma)
#' ## ------------------------- ##

# Function to load one by one, rbind and store
loader <- function(i){
  load(paste0('~/Downloads/ests2-',i,'.RData'))
  df <- as_tibble(out$ests) %>% pivot_longer(everything()) %>% mutate(a = as.character(i))
  message(i)
  df
}

df <- list()
for(i in 1:4){df[[i]] <- loader(i)}
df <- do.call(rbind, df)

df <- df %>% 
  mutate(
    description = case_when(
      a == '1' ~ 'ntms = 5',
      a == '2' ~ 'ntms = 10',
      a == '3' ~ 'ntms = 15',
      a == '4' ~ 'ntms = 15, lower failure',
      a == '5' ~ 'ntms = 15, lower D',
      a == '6' ~ 'ntms = 15, lower failure, lower D',
      a == '7' ~ 'ntms = 15, lower failure, lower D, +ve gamma',
      T ~ 'AA'
    )
  )

targets1 <- data.frame(
  name = unique(df$name),
  target = c(.5, 0.0, 0.1, 0.5, -0.2, 0.1, 0.2, -1.0, -0.1, -0.50, 0.05, -0.30)
)

targets2 <- data.frame(
  name = unique(df$name),
  target = c(.5^2, 0.0, 0.05^2, 0.5, -0.2, 0.1, 0.2, -1.0, -0.1, -0.50, 0.05, -0.30)
)

df %>% 
  filter(a %in% c('1','2','3','4')) %>% #, name != 'gamma') %>%
  
  left_join(., targets1, 'name') %>% 
  ggplot(aes(x = value, colour = description)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~name, scales = 'free') + 
  theme_bw()

df %>% 
  filter(a %in% c('5', '6')) %>% 
  left_join(., targets2, 'name') %>% 
  ggplot(aes(x = value, colour = description)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~name, scales = 'free') + 
  theme_bw()
