library(dplyr)
pbc2 <- joineRML::pbc2
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

# Number of NA's per column

apply(pbc, 2, function(x) length(which(is.na(x))))
# Ascites, Spiders, Hepatomegaly (three binary), serum cholesterol and platelets all have missing values.



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
  colSums() %>% .[-1]

changes %>% # backed up here
  apply(.,2,function(x) x > 0) %>% 
  colSums() %>% .[-1]

changes %>% # backed up here
  apply(.,2,function(x) x > 1) %>% 
  colSums() %>% .[-1]

changes %>% # and here
  apply(.,2,function(x) x > 2) %>% 
  colSums() %>% .[-1]

pbc <- pbc %>% 
  mutate_at(vars(hepatomegaly, ascites, spiders), ~ as.numeric(.x) - 1)

# Make our set to export
pbc <- pbc %>% 
  select(id, survtime, status, cont, bin, time, Y.1 = albumin, Y.2 = spiders, Y.3 = platelets) %>% 
  mutate_at('Y.2', ~ as.numeric(.x)) %>% 
  mutate_at('id', ~ as.numeric(.x)) %>% 
  filter(complete.cases(.))

if(any(diff(pbc$id) > 1)) stop('IDs non-constant, add line of code!')

survdata <- distinct(pbc, id, survtime, status, cont, bin)

pbcdata <- list(pbc = as.data.frame(pbc), survdata = as.data.frame(survdata))
save(pbcdata, file = 'PBC/pbc-hepatomegaly.RData')



source()
pbc <- as.data.frame(pbc)

fit <- EM(pbc, ph, survdata, quad = F)
fit <- EM(pbc, ph, survdata, quad = T, tol = 1e-3)


glmmTMB(Y.3 ~ time + cont + bin + (1+time|id), pbc, family = poisson) -> temp

library(JMbayes2)
m1 <- lme(Y.1 ~ time + cont + bin, random = ~ time|id, pbc, method = 'ML')
m2 <- mixed_model(Y.2 ~ time + cont + bin, random = ~ 1|id, pbc, family = binomial)
m3 <- mixed_model(Y.3 ~ time + cont + bin, random = ~ time|id, pbc, family = poisson,
                  initial_values = list(
                    betas = fixef(temp)$con,
                    D = matrix(VarCorr(temp)$cond$id,2,2)
                  ))
M <- list(m1,m2,m3)

jmb <- jm(ph, M, 'time',n_chains = 1L,
          n_iter = 11000L, n_burnin = 1000L)

