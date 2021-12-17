library(JMbayes2)
library(dplyr)
pbc <- pbc2 %>% as_tibble %>% 
  select(-status) %>% 
  rename(survtime = years, time = year, status = status2)
pbc$bin <- as.numeric(pbc$drug) - 1

survdata <- pbc %>%
  distinct(id, survtime, status, age, bin) %>% 
  mutate(`cont` = as.numeric(scale(age)))

pbc <- left_join(
  pbc %>% select(-age, -drug),
  survdata %>% select(id, cont), 
  'id'
)

changes <- pbc %>% 
  group_by(id) %>% 
  mutate(
    across(ascites:spiders, ~ ifelse(.x != lag(.x), 1, 0)),
    r = row_number(),
  ) %>% select(id, time, ascites, hepatomegaly, spiders) %>% 
  filter(!is.na(ascites), !is.na(hepatomegaly), !is.na(spiders)) %>% 
  mutate(ascites = cumsum(ascites),
         hepatomegaly = cumsum(hepatomegaly),
         spiders = cumsum(spiders)) %>% 
  mutate(across(ascites:spiders, ~ max(.x))) %>% ungroup %>% select(-time) %>% 
  distinct 

# Ascites changes the least, spiders second least and hepatomegaly the most.
changes %>% 
  apply(.,2,function(x)x==0) %>% 
  colSums()

changes %>% # backed up here
  apply(.,2,function(x) x > 0) %>% 
  colSums()

changes %>% # and here
  apply(.,2,function(x) x > 2) %>% 
  colSums()

# Make our set to export
pbc <- pbc %>% 
  select(id, survtime, status, cont, bin, time, Y.1 = albumin, Y.2 = hepatomegaly, Y.3 = platelets) %>% 
  mutate_at('Y.2', ~ as.numeric(.x) - 1) %>% 
  mutate_at('id', ~ as.numeric(.x)) %>% 
  filter(complete.cases(.))

survdata <- distinct(pbc, id, survtime, status, cont, bin)

pbcdata <- list(pbc = as.data.frame(pbc), survdata = as.data.frame(survdata))
save(pbcdata, file = 'pbc2.RData')


# Fits
source('EM.R')
ph <- coxph(Surv(survtime, status) ~ cont + bin, pbcdata$survdata)
my.fit3 <- EM(pbcdata$pbc, ph, pbcdata$survdata, gh = 3)
save(my.fit3, file = 'myfit2.RData')

# JMBayes
library(glmmTMB)
fitglmmtmb <- glmmTMB(Y.3 ~ time + cont + bin + (1 + time|id), data =pbcdata$pbc, family = poisson())

m1 <- lme(fixed = Y.1 ~ time + cont + bin, random = ~ time | id, data = pbcdata$pbc)
m2 <- mixed_model(Y.2 ~ time + cont + bin, random = ~ time | id, data = pbcdata$pbc, family = binomial())
m3 <- mixed_model(Y.3 ~ time + cont + bin, random = ~ time | id, data = pbcdata$pbc, family = poisson(),
                  initial_values = list(
                    betas = c(fixef(fitglmmtmb)$cond),
                    D = matrix(VarCorr(fitglmmtmb)$cond$id, 2, 2)
                  ))
M <- list(m1, m2, m3)

jmb.fit <- jm(ph, M, time_var = 'time', 
              n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)

# Comparing different GLMMs for Y.3
library(bbmle)
fit1 <- glmmTMB(Y.3 ~ time + cont + bin + (1 + time|id), data =pbcdata$pbc, family = poisson())
fit2a <- glmmTMB(Y.3 ~ time + cont + bin + (1 + time|id), 
                 dispformula = ~ 1 + time, data =pbcdata$pbc, family = glmmTMB::nbinom1())
fit2b <- glmmTMB(Y.3 ~ time + cont + bin + (1 + time|id), 
                 dispformula = ~ 1 + time, data =pbcdata$pbc, family = glmmTMB::nbinom2())
fit3 <- glmmTMB(Y.3 ~ time + cont + bin + (1 + time|id), data =pbcdata$pbc, family = glmmTMB::compois())

AICtab(fit1, fit2a, fit2b)


AIC(fit2b)

# run

# MCMC summary:
#   chains: 1 
# iterations per chain: 11000 
# burn-in per chain: 1000 
# thinning: 1 
# 
