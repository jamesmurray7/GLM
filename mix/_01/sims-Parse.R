# Extract out individual fits and store separately ------------------------
save.location <- '~/Downloads/'
vech <- function(X) X[lower.tri(X, diag = T)]
# Functions
extract.coeffs <- function(x){ # x a 'sub-list'
  xc <- x$coeffs
  vD <- vech(xc$D)
  names(vD) <- paste0('[' , apply(which(lower.tri(xc$D, T), arr.ind = T), 1, paste0, collapse = ', '), ']')
  beta <- xc$beta
  names(beta) <- paste0(rep(c('G_', 'B_', 'P_'), each = 4), c('(Intercept)', 'time', 'cont', 'bin'))
  out <- c(vD, beta, 'var.e' = xc$var.e, xc$gamma, xc$eta)
  out
}

# Unpacking function
.unpack <- function(x){ # x a list of lists.
  for(i in 1:3){
    this <- x[[i]]
    ests <- do.call(rbind, lapply(this, function(x){
      if(!is.null(x)) extract.coeffs(x)
    }))
    out <- ests
    message('\nfits[[',i,']] had ', nrow(ests), ' successful fits out of 100.')
    message('\nSaving in ', save.location, 'mixests-', i, '.RData')
    save(out, file = paste0(save.location, 'mixests-', i, '.RData'))
  }
}

.unpack(fits)

# Now we can parse --------------------------------------------------------
source('./simData.R')
library(tidyverse)
#' ## Table of what set is what ##
#' ## beta = c(1, 0.1, 0.33, -0.5), 
#' ## D = matrix(c(0.5, 0, 0, 0.1), 2, 2), 
#' ##gamma = 0.5, surv.eta = c(0.05, -0.3)
#' ## ------------------------- ##
#' ## Fits1: from ntms = 10; 2: ntms = 15; 3: ntms = 10, gamma = -1

# Function to load one by one, rbind and store
loader <- function(i){
  load(paste0('~/Downloads/mixests-',i,'.RData'))
  df <- as_tibble(out) %>% pivot_longer(everything()) %>% mutate(a = as.character(i))
  message(i)
  df
}

df <- list()
for(i in 1:3){df[[i]] <- loader(i)}
df <- do.call(rbind, df)

df <- df %>% 
  mutate(
    description = case_when(
      a == '1' ~ 'ntms = 10',
      a == '2' ~ 'ntms = 15',
      a == '3' ~ 'ntms = 15, lower failure',
      T ~ 'AA'
    )
  )

targets <- data.frame(
  name = unique(df$name),
  target = c(vech(true.D), 
             do.call(c, lapply(1:3, function(i) true.beta[i,])),
             0.25, true.gamma, true.eta)
)

df %>% 
  left_join(., targets, 'name') %>% 
  filter(!grepl('^\\[', name)) %>% 
  ggplot(aes(x = value, colour = description)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~fct_inorder(name), scales = 'free', nc = 4) + 
  theme_bw()

ggsave('~/Downloads/mix-ests.png')

# Covariance only with variance

df %>% 
  left_join(., targets, 'name') %>% 
  filter(grepl('^\\[|^var', name)) %>% 
  ggplot(aes(x = value, colour = description)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~fct_inorder(name), scales = 'free', nc = 4) + 
  theme_bw()

ggsave('~/Downloads/mix-ests-D.png')

# Timings -----------------------------------------------------------------

aa <- function(x) list(EM = x$EMtime, total = x$totaltime)
ntms10 <- lapply(fits[[1]], aa)
ntms15 <- lapply(fits[[2]], aa)
ntms152 <- lapply(fits[[3]], aa)

EM <- data.frame(`ntms10` = do.call(c, lapply(ntms10, '[[', 1)),
                 `ntms15` = do.call(c, lapply(ntms15, '[[', 1)),
                 `ntms = 15, lower failure` = do.call(c, lapply(ntms152, '[[', 1)))
total <- data.frame(`ntms10` = do.call(c, lapply(ntms10, '[[', 2)),
                 `ntms15` = do.call(c, lapply(ntms15, '[[', 2)),
                 `ntms = 15, lower failure` = do.call(c, lapply(ntms152, '[[', 2)))

EM$a <- "Approximate EM"
total$a <- "Total Computation Time"

both <- rbind(EM, total) %>% as_tibble 
both %>% 
  rename(`ntms = 15, lower failure` = "ntms...15..lower.failure",
         `ntms = 10` = ntms10, `ntms = 15` = ntms15) %>% 
  pivot_longer( `ntms = 10`:`ntms = 15, lower failure`) %>% 
  ggplot(aes(x = name, y = value, fill = a)) + 
  geom_boxplot() + 
  labs(fill = NULL, x = NULL, y = "Elapsed time (s)") + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) + 
  scale_fill_brewer(palette  = 'Spectral') + 
  theme_bw() + 
  theme(legend.position = 'bottom') 
ggsave('ET.png')
