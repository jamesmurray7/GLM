#' ###
#' Program which fits all univariate GLMM joint models to PBC data.
#' ###
rm(list=ls())
library(dplyr)
source('EM.R')
source('univariate-compare.R')
load('ADNI.RData')


# Define a function that 'resets' the data --------------------------------
newadni <- function(Y){
  y <- adni[,Y]
  inds.to.remove <- is.na(y)
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

# Establishing a 'global' proportional hazards model
coxph(Surv(survtime, status) ~ APOE4 + age_scaled + PTGENDER, adni)
surv.formula <- Surv(survtime, status) ~ APOE4 + age_scaled + gender

#' ####
#' GAUSSIANS
#' i.e. (log-)normal or Box-Cox transformed Normals
#' ####

# Mid temp ----------------------------------------------------------------
midtemp <- newadni('MidTemp')
midtemp$MidTemp <- c(scale(midtemp$MidTemp))
# Transform PTGENDER to be numeric
midtemp$gender <- as.numeric(midtemp$PTGENDER) - 1 # 1: Female
long.formula <- MidTemp ~ APOE4 + time + age_scaled + PTGENDER + (1 + time|id)
summary(glmmTMB(long.formula, midtemp))

my.fit <- EM(long.formula, surv.formula, midtemp, family = gaussian)
my.summary(my.fit)

# ADAS13 ------------------------------------------------------------------
adas13 <- newadni('ADAS13')
plot.long(adas13, ADAS13, PTGENDER)
long.formula <- ADAS13 ~ APOE4 * time + age_scaled + PTGENDER + (1 + time|id)
my.fit <- EM(long.formula, Surv(survtime, status) ~ APOE4, adas13, poisson,
             control = list(hessian = 'auto'))

adas13$drug <- adas13$APOE4 # lazy
jm.fit <- JMbayes2.fit(adas13, poisson, ADAS13 ~ drug * time, ~time|id, Surv(survtime, status) ~ drug)

my.summary(my.fit)
summary(jm.fit)

# ADAS11 ------------------------------------------------------------------

adas11 <- newadni('ADAS11')
plot.long(adas13, ADAS11, APOE4)
long.formula <- ADAS11 ~ APOE4 * time + age_scaled + PTGENDER + (1 + time|id)
my.fit <- EM(long.formula, Surv(survtime, status) ~ APOE4, adas11, poisson,
             control = list(hessian = 'auto'))

adas11$drug <- adas11$APOE4 # lazy
jm.fit <- JMbayes2.fit(adas11, poisson, ADAS11 ~ drug * time + age_scaled + PTGENDER, ~time|id, Surv(survtime, status) ~ drug)

summary(jm.fit)
my.summary(my.fit)

# MMSE --------------------------------------------------------------------
mmse <- newadni("MMSE")
plot.long(mmse, MMSE, as.factor(status))
long.formula <- MMSE ~ APOE4 * time + age_scaled + PTGENDER + (1|id)

summary(glmmTMB(long.formula, mmse, family = poisson)) # Seems quite unstable

my.fit <- EM(long.formula, Surv(survtime, status) ~ APOE4, mmse, poisson,
             control = list(hessian = 'auto'))


# RAVLT-Immediate ---------------------------------------------------------
rav <- newadni("RAVLT.immediate")
plot.long(rav, RAVLT.immediate, as.factor(APOE4))

long.formula <- RAVLT.immediate ~ APOE4 * time + age_scaled + PTGENDER + (1 + time|id)

one <- glmmTMB(long.formula, rav, family = poisson)
two <- update(one, . ~ . -APOE4)
anova(one, two)

my.fit <- EM(RAVLT.immediate ~ time + APOE4:time + age_scaled + PTGENDER + (1 + time|id),
             Surv(survtime, status) ~ APOE4, rav, family = poisson, control = list(hessian = 'auto'))

rav$drug <- rav$APOE4
jm.fit <- JMbayes2.fit(rav, poisson,
                       RAVLT.immediate ~ time + drug:time + age_scaled + PTGENDER, ~time|id,
                       Surv(survtime, status) ~ drug)
summary(jm.fit)
my.summary(my.fit) # My ests more in-keeping with initial GLMMTMB ests?



