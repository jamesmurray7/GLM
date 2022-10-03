#' #########
#' Approximate EM fits using half_and_half2 MPCMP implementation.
#' #########

rm(list=ls())
if(unname(Sys.info()['user']=='c0061461')) setwd('~/Documents/GLMM/mpcmp/half_and_half2/') else setwd('~/Documents/PhD/GLM/mpcmp/half_and_half2/')
load('../../PBC-case-study/ADNI.RData') # Same across Linux/Mac(!)
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

# EM fits (MPCMP) ---------------------------------------------------------

surv.formula <- Surv(survtime, status) ~ bin

control <- list(verbose=F, debug = T)
disp.control <- list(delta.method = 'optim', min.profile.length = 2,
                     truncated = T, max.val = 2.5)
optimiser.arguments <- optim.control <- list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3)
surv.formula <- Surv(survtime, status) ~ bin

summax.fn <- function(y) max(y) + 10
min.summax <- 20

f <- function(Y){
  d <- newadni(Y)
  if(Y == 'RAVLT.forgetting') RE <- '(1|id)' else RE <- '(1+time|id)'
  
  file.name <- paste0('./ADNI/logs/', gsub('\\.', '_', Y), '_mpcmp2.log')
  
  long.formula <- as.formula(paste0(Y, ' ~ time + age_scaled + bin +', RE))
  
  message('Starting fit for ', Y)
  fit <- EM(long.formula, surv.formula, d,
             control = control,
             disp.control = disp.control,
             optim.control = optim.control,
             summax.fn = summax.fn, min.summax = min.summax,
             delta.update.quad = T,
             beta.update.quad = F,
             initialise.delta = F)
  
  sink(file.name)
  my.summary(fit, T, T)
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
