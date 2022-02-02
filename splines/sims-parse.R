# Extract out individual fits and store separately ------------------------
save.location <- '~/Downloads/'
vech <- function(X) X[lower.tri(X, diag = T)]
# Functions
extract.coeffs <- function(x){ # x a 'sub-list'
  xc <- x$coeffs
  vD <- vech(xc$D)
  names(vD) <- paste0('[' , apply(which(lower.tri(xc$D, T), arr.ind = T), 1, paste0, collapse = ', '), ']')
  beta <- xc$beta
  names(beta) <- paste0(rep(paste0('beta', 1:3, '_'), each = 5),
                        c('(Intercept)', 'time', 'time^2', 'cont', 'bin'))
  var.e <- xc$var.e; names(var.e) <- paste0('var.e_', 1:3)
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
    message('\nSaving in ', save.location, 'quadests-', i, '.RData')
    save(out, file = paste0(save.location, 'quadests-', i, '.RData'))
  }
}

load('~/Downloads/quadfits.RData')
.unpack(fits)

# Now we can parse --------------------------------------------------------
source('./simData.R')
library(tidyverse)
#' Fits1: from ntms = 10; 2: ntms = 6; 3: ntms = 14

true.beta <- c(1, -0.2, 0.01, 0.33, -0.50,
               0, -0.5, 0.05, -0.33, 0.50,
               3, 0.1, -0.05, 0.5, 0.1)
true.var <- rep(.25, 3)
true.D <- as.matrix(Matrix::bdiag(replicate(3, diag(c(0.5^2, .2^2, .05^2)), simplify = F)))
true.gamma <- c(0.50, -0.25, 0.40)
true.eta  <-  c(0.05, -0.3)
# Function to load one by one, rbind and store
loader <- function(i){
  load(paste0('~/Downloads/quadests-',i,'.RData'))
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
      a == '2' ~ 'ntms = 6',
      a == '3' ~ 'ntms = 14',
      T ~ 'AA'
    )
  )

targets <- data.frame(
  name = unique(df$name),
  target = c(vech(true.D), true.beta, true.var, true.gamma, true.eta)
)

df %>% 
  left_join(., targets, 'name') %>% 
  filter(!grepl('^\\[', name), a != 2) %>% 
  ggplot(aes(x = value, colour = description)) +
  geom_vline(aes(xintercept = target)) +
  geom_density(lwd = 0.8, alpha = .75) + 
  facet_wrap(~fct_inorder(name), scales = 'free', nc = 4) + 
  theme_bw()

ggsave('~/Downloads/quad-ests.png')

# Timings -----------------------------------------------------------------

aa <- function(x) list(EM = x$EMtime, total = x$totaltime)
ntms10 <- lapply(fits[[1]], aa)
ntms6 <- lapply(fits[[2]], aa)
ntms14 <- lapply(fits[[3]], aa)

EM <- data.frame(`ntms10` = do.call(c, lapply(ntms10, '[[', 1)),
                 #`ntms6` = do.call(c, lapply(ntms6, '[[', 1)),
                 `ntms = 14` = do.call(c, lapply(ntms14, '[[', 1)))
total <- data.frame(`ntms10` = do.call(c, lapply(ntms10, '[[', 2)),
                  #  `ntms6` = do.call(c, lapply(ntms6, '[[', 2)),
                    `ntms = 14` = do.call(c, lapply(ntms14, '[[', 2)))

EM$a <- "Approximate EM"
total$a <- "Total Computation Time"

both <- rbind(EM, total) %>% as_tibble 
both %>% 
  rename(`ntms = 14` = "ntms...14",
         `ntms = 10` = ntms10) %>% 
  pivot_longer(`ntms = 10`:`ntms = 14`) %>% 
  ggplot(aes(x = name, y = value, fill = a)) + 
  geom_boxplot() + 
  labs(fill = NULL, x = NULL, y = "Elapsed time (s)") + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) + 
  scale_fill_brewer(palette  = 'Spectral') + 
  theme_bw() + 
  theme(legend.position = 'bottom') 
ggsave('ET.png')
