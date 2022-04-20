rm(list=ls(
  
))
source('EM.R')
load('~/Downloads/allfits.RData')

# Functions ---------------------------------------------------------------
get.estimates <- function(y){ # x a sub-list
    if(!is.null(y)){
      z <- c(vech(y$coeff$D), c(y$co$beta), y$co$gamma, y$co$zeta)
      if(y$co$sigma != 0) z <- c(z, y$co$sigma)
      z <- setNames(z, names(y$SE))
    }else{
      z <- NA
    }
    z
}

get.SE <- function(x) if(!is.null(x)) x$SE else NA

library(tidyverse)
for(f in c('G', 'B', 'P', 'N')){
  switch(f,
         'G' = family <- 'Gaussian', 
         'B' = family <- 'Binomial',
         'P' = family <- 'Poisson',
         'N' = family <- 'Negative binomial')
  
  # Estimates
  ests <- do.call(rbind, lapply(fits[[f]], get.estimates))
  rS <- rowSums(ests)
  ests <- ests[!is.na(rS), ]
  
  # Standard Errors
  SE <- do.call(rbind, lapply(fits[[f]], get.SE))
  SE <- SE[!is.na(rS), ]
  
  # Lower and upper bounds
  lb <- ests - qnorm(.975) * SE; ub <- ests + qnorm(.975) * SE
  # Targets
  targets <- c(0.5, 0.0, 0.1, 1.0, 0.1, 0.33, -0.5, 0.5, 0.05, -0.30)
  if(f == 'G') targets <- c(targets, 0.16)
  if(f == 'N') targets <- c(targets, 2.00)
  names(targets) <- colnames(SE)
  targets.mat <- apply(t(targets), 2, rep, nrow(SE))
  # Coverage
  CP <- colSums(lb <= targets.mat & ub >= targets.mat)/nrow(SE)
  
  ests.long <- pivot_longer(as.data.frame(ests), `D[1,1]`:colnames(ests)[ncol(ests)], values_to = 'estimate', names_to = 'parameter')
  SE.long <- pivot_longer(as.data.frame(SE), `D[1,1]`:colnames(SE)[ncol(SE)], values_to = 'SE', names_to = 'parameter')
  
  ests <- cbind(ests.long, SE.long[, -1])
  ests <- left_join(ests, as.data.frame(CP) %>% rownames_to_column('parameter'), 'parameter')
  ests <- left_join(ests, as.data.frame(targets) %>% rownames_to_column('parameter'), 'parameter')
  ests <- ests %>% 
    mutate(bias = estimate - targets) %>% 
    group_by(parameter) %>% 
    mutate(mbias = mean(bias),
           sbias = sd(bias)) %>% 
    ungroup()
  
 ests$param <- gsub(',', '', ests$parameter)
 ests$param <- gsub('\\_', '[', ests$param)
 ests$param <- gsub('\\(Intercept\\)', '"(Intercept)"]', ests$param)
 ests$param <- gsub('time', '"time"]', ests$param)
 ests$param <- gsub('cont', '"cont"]', ests$param)
 ests$param <- gsub('bin', '"bin"]', ests$param)
 
 ests$CP.annotate <- paste0('CP: ', format(round(ests$CP, 2), nsmall = 2))

  ests %>% 
    ggplot(aes(x = estimate)) + 
    geom_vline(aes(xintercept = targets), lty = 5) +
    geom_density() + 
    facet_wrap(~param, scales = 'free', labeller = label_parsed) +
    geom_text(data = distinct(ests, parameter, param, CP.annotate),
              aes(x = Inf, y = Inf, label = CP.annotate),
              colour = 'navy',
              hjust = 1.0,
              vjust = 1.0) + 
    labs(title = family, y = '', x = 'Estimate') + 
    theme_bw() + 
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size=11),
      plot.title = element_text(hjust=.4)
    )
  ggsave(file = paste0('~/Downloads/', family, '-ests.png'))
  
}
