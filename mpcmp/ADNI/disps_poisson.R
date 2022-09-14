#' #########
#' Getting handles dispersions of ADNI variables.
#' #########
rm(list=ls())
if(unname(Sys.info()['user']=='c0061461')) setwd('~/Documents/GLMM/Multi-test/') else setwd('~/Documents/PhD/GLM/Multi-test/')
load('../PBC-case-study/ADNI.RData') # Same across Linux/Mac(!)
source('EM.R')
library(dplyr)
# Define a function that 'resets' the data --------------------------------
newadni <- function(Y){
  y <- adni[,Y]
  inds.to.remove <- is.na(y) | y < 0
  if(length(inds.to.remove) > 0){
    rtn <- adni[!inds.to.remove, ]
    rtn <- rtn %>% 
      group_by(id) %>% 
      mutate(new_id = cur_group_id()) %>% 
      ungroup %>% 
      select(-id) %>% 
      rename(id = new_id) %>% 
      as.data.frame()
  }else{
    rtn <- adni
  }
  rtn
}

TMBs <- function(Y){
  d <- newadni(Y)
  one <- glmmTMB(
    as.formula(paste0(Y, ' ~ time + age_scaled + bin + (1 + time|id)')),
    data = d, family = 'poisson'
  )
  two <- glmmTMB(
    as.formula(paste0(Y, ' ~ time + age_scaled + bin + (1|id)')),
    data = d, family = 'poisson'
  )
  c(intslope = AIC(one), int = AIC(two))
}

adni$MMSE.reverse <- 30 - adni$MMSE

responses <- c('ADAS11', 'ADAS13', 'MMSE', 'MMSE.reverse', 'RAVLT.immediate', 'RAVLT.learning', 'RAVLT.forgetting', 'FAQ')

tmbs <- sapply(responses, TMBs)
apply(tmbs, 2, which.min)
# ALL int-slope besides MMSE (which i don't think we can fit a model to!) and RAVLT.forgetting...

# EM fits (Poisson) -------------------------------------------------------

surv.formula <- Surv(survtime, status) ~ bin
family <- list('poisson')

f <- function(Y){
  d <- newadni(Y)
  if(Y == 'RAVLT.forgetting') RE <- '(1|id)' else RE <- '(1+time|id)'
  
  
  file.name <- paste0('../mpcmp/ADNI/logs/poisson/', gsub('\\.', '_', Y), '_poisson.log')
  
  long.formula <- list(as.formula(paste0(Y, ' ~ time + age_scaled + bin +', RE)))
  
  message('Starting fit for ', Y)
  fit <- suppressMessages(EM(long.formula, surv.formula, d, family))
  
  sink(file.name)
  my.summary(fit, T)
  sink()
  
  message('Saved in ', file.name, '\n')
  invisible(1+1)
}

# ADAS11
f('ADAS11')

# ADAS13
f('ADAS13')

# MMSE (reversed)
f('MMSE.reverse')

# RAVLT.immediate
f('RAVLT.immediate')

# RAVLT.learning
f('RAVLT.learning')

# RAVLT.forgetting
f('RAVLT.forgetting')

# FAQ
f('FAQ')
