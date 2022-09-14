#' #########
#' Getting handles dispersions of ADNI variables.
#' #########
rm(list=ls())
if(unname(Sys.info()['user']=='c0061461')) setwd('~/Documents/GLMM/mpcmp/no-grids/') else setwd('~/Documents/PhD/GLM/mpcmp/no-grids/')
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


# EM fits -----------------------------------------------------------------
surv.formula <- Surv(survtime, status) ~ APOE4
disp.formula <- ~ 1
# ADAS11
d <- newadni('ADAS11')
long.formula <- ADAS11 ~ time + age_scaled + APOE4 + (1 + time|id)
ADAS11.ng <- EM(long.formula, disp.formula, surv.formula,
                d, control = list(verbose = T), disp.control = list(delta.method = 'bobyqa'),
                optim.control = list(optimiser = 'optim', Hessian = 'obj'))
sink('../ADNI/logs/ADAS11_no-grids.log')
my.summary(ADAS11.ng, T)
sink()

# ADAS13
d <- newadni('ADAS13')
long.formula <- ADAS13 ~ time + age_scaled + bin + (1+time|id)
ADAS13.ng <- EM(long.formula, disp.formula, surv.formula, summax = 100,
                data = d, control = list(verbose = T, auto.summax = T), optimiser = 'bobyqa', Hessian = 'obj')
sink('../ADNI/logs/ADAS13_no-grids.log')
my.summary(ADAS13.ng, T)
sink()

# RAVLT.forgetting
d <- newadni('RAVLT.forgetting')
TMBs('RAVLT.forgetting') # intonly margianlly preferable.
long.formula <- RAVLT.forgetting ~ time + age_scaled + APOE4 + (1|id)
RAVLT.forgetting.ng <- EM(long.formula, disp.formula, surv.formula,
                d, control = list(verbose = T), disp.control = list(delta.method = 'bobyqa', max.val = Inf),
                optim.control = list(optimiser = 'bobyqa', Hessian = 'obj'))
sink('../ADNI/logs/RAVLT_Forgetting_no-grids.log')
my.summary(RAVLT.forgetting.ng, T)
sink()

# MMSE
d <- newadni('MMSE') # Erroring... random effects are absolutely tiny.
d$MMSE.reverse <- 30 - d$MMSE
TMBs('MMSE.reverse')
long.formula <- MMSE ~ time + age_scaled + APOE4 + (1|id)
io <- glmmTMB(long.formula, d, 'poisson')
io2 <- glmmTMB(update(long.formula, MMSE.reverse ~ .), d, 'poisson')
# io2b <- glmmTMB(MMSE ~ time + age_scaled + APOE4 + (1+time|id), d, 'poisson') # intonly still...
MMSE.ng <- EM(update(long.formula, MMSE.reverse ~ .), disp.formula, surv.formula,
              d, control = list(verbose = T), disp.control = list(delta.method = 'bobyqa', max.val = Inf),
              optim.control = list(optimiser = 'bobyqa', Hessian = 'obj'))
sink('../ADNI/logs/MMSE_Reversed_no-grids.log')
my.summary(MMSE.ng, T)
sink()
