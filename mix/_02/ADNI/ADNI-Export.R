#' #####
#' ADNI-Export
#' ----
#' Creating subset of ADNI data to use with mixture GLMM-JM fit.
#' #####

rm(list=ls())
library(tidyverse)
load("ADNI_JM.RData")

ADNI_JM$month <- as.numeric(as.character(ADNI_JM$Month))
ADNI_JM$Year <- ADNI_JM$month/12

set <- ADNI_JM %>% 
  select(RID, age = AGE, PTGENDER, APOE4,
         month, Year, ADAS13, FAQ, MidTemp, EventTime, cens, RAVLT.perc.forgetting) %>% 
  mutate_at('RAVLT.perc.forgetting', ~ abs(.x)) %>% 
  filter(complete.cases(.),
         EventTime >= Year,
         RAVLT.perc.forgetting <= 100) 

set$gender = as.numeric(set$PTGENDER) - 1
set$bin <- ifelse(set$APOE4 >= 1, 1, 0)

set2 <- set %>% 
  mutate(Y.1 = c(scale(MidTemp)),
         Y.2 = ifelse(RAVLT.perc.forgetting >= 80, 1, 0),
         Y.3 = floor(ADAS13),
         time = Year,
         status = cens,
         survtime = EventTime)

survdata <- set2 %>% distinct(RID, age, bin, survtime, status)
survdata$cont <- c(scale(survdata$age))
survdata$id <- 1:nrow(survdata)
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)

set3 <- left_join(set2, survdata %>% select(RID, id, cont), 'RID')

data <- set3 %>% select(id, cont, bin, survtime, status, time, Y.1, Y.2, Y.3) 

# Make sure there's a period of time between final t and T.
data <- data %>% group_by(id) %>% mutate(tempflag = any(time == survtime)) %>% 
  ungroup %>% mutate(survtime = ifelse(tempflag, survtime +  1e-2, survtime)) %>% 
  select(-tempflag)
# Ensuring no times are duplicated: Staggering by some small time 1e-1.
data <- data %>% group_by(id) %>% mutate(ltime = lag(time)) %>% 
  mutate(time2 = ifelse(time == ltime, time + 1e-1, time),
         time = ifelse(!is.na(time2), time2, time)) %>% ungroup %>% select(-time2, -ltime)

survdata <- data %>% distinct(id, cont, bin, survtime, status)  
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)

ADNI <- list(
  data = data,
  survdata = survdata,
  ph = ph
)
save(ADNI, file = 'ADNI80.RData')
