#' #####
#' ADNI-Export
#' ----
#' Creating subset of ADNI data to use with mixture GLMM-JM fit.
#' #####

rm(list=ls())
library(dplyr)
load("../mix/_01/ADNI/ADNI_JM.RData")

ADNI_JM$month <- as.numeric(as.character(ADNI_JM$Month))
ADNI_JM$Year <- ADNI_JM$month/12

set <- ADNI_JM %>% 
  mutate_at('RAVLT.perc.forgetting', ~ abs(.x)) %>% 
  filter(EventTime >= Year) 

set$gender = as.numeric(set$PTGENDER) - 1
set$bin <- ifelse(set$APOE4 >= 1, 1, 0)
set$RAVLT.perc.forgetting[which(set$RAVLT.perc.forgetting>100)] <- NA
set$RAVLT.80pc <- ifelse(set$RAVLT.perc.forgetting >= 80, 1, 0)
set$RAVLT.100pc <- ifelse(set$RAVLT.perc.forgetting == 100, 1, 0)
set$ADAS11 <- floor(set$ADAS11)
set$ADAS13 <- floor(set$ADAS13)
set$APOE4 <- ifelse(set$APOE4 > 1, 1, 0)

set2 <- set %>% 
  rename(time = Year,
         status = cens,
         survtime = EventTime)

survdata <- set2 %>% distinct(RID, AGE, survtime, status)
survdata$age_scaled <- c(scale(survdata$AGE))
survdata$id <- 1:nrow(survdata)

set3 <- left_join(set2, survdata %>% select(RID, id, age_scaled), 'RID')

data <- set3

# Make sure there's a period of time between final t and T.
data <- data %>% group_by(id) %>% mutate(tempflag = any(time == survtime)) %>% 
  ungroup %>% mutate(survtime = ifelse(tempflag, survtime + 1e-2, survtime)) %>% 
  select(-tempflag)
# Ensuring no times are duplicated: Staggering by some small time 1e-1.
data <- data %>% group_by(id) %>% mutate(ltime = lag(time)) %>% 
  mutate(time2 = ifelse(time == ltime, time + 1e-1, time),
         time = ifelse(!is.na(time2), time2, time)) %>% ungroup %>% select(-time2, -ltime)

adni <- as.data.frame(data)
save(adni, file = 'ADNI.RData')
