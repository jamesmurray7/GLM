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



# EM fits -----------------------------------------------------------------
surv.formula <- Surv(survtime, status) ~ APOE4
disp.formula <- ~ 1
# ADAS11
d <- newadni('ADAS11')
long.formula <- ADAS11 ~ time + age_scaled + APOE4 + (1 + time|id)
ADAS11.ng <- EM(long.formula, disp.formula, surv.formula,
                     d, control = list(verbose = T))
sink('../ADNI/ADAS11_no-grids.log')
my.summary(ADAS11.ng, T)
sink()

# ADAS13
d <- newadni('ADAS13')
long.formula <- ADAS13 ~ time + age_scaled + bin + (1+time|id)
ADAS13.ng <- EM(long.formula, disp.formula, surv.formula, summax = 100,
                data = d, control = list(verbose = T, auto.summax = T), optimiser = 'bobyqa', Hessian = 'obj')
# Issue with betas...
sink('../ADNI/logs/ADAS13_no-grids.log')
my.summary(ADAS13.ng, T)
sink()

# MMSE
d <- newadni('MMSE') # Erroring... random effects are absolutely tiny.
long.formula <- MMSE ~ time + age_scaled + APOE4 + (1|id)
MMSE.hybrid2 <- EM(long.formula, disp.formula, surv.formula, d, 
                   control = list(verbose = T), 
                   disp.control = list(interval = c(-2, 2)), 
                   grid.summax = 'same')

# RAVLT.immediate
d <- newadni('RAVLT.immediate')
long.formula <- RAVLT.immediate ~ time + age_scaled + APOE4 + (1|id)
RAVLT.immediate.hybrid2 <- EM(
  long.formula, disp.formula, surv.formula, d,
  control = list(verbose = T), grid.summax = 'a', optimiser = 'bobyqa'
)



