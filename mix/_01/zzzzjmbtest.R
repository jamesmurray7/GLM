load('~/Downloads/jmfit.RData')
JMbayes2::traceplot(jmb.fit)

source('EM.R')

set.seed(123)
data <- simData_joint()
ph <- coxph(Surv(survtime, status) ~ cont + bin, data$surv.data)

my.fit <- EM(data$data, ph, data$surv.data, gh = 3)
my.fit9 <- EM(data$data, ph, data$surv.data, gh = 9)

library(JMbayes2)
m1 <- lme(Y.1 ~ time + cont + bin, random = ~ time | id, data = data$data)
m2 <- mixed_model(Y.2 ~ time + cont + bin, random = ~ time | id, data = data$data, family = binomial())
m3 <- mixed_model(Y.3 ~ time + cont + bin, random = ~ time | id, data = data$data, family = poisson())


M <- list(m1, m2, m3)
jmb.fit <- jm(ph, M, 'time', n_chains = 1L, n_burnin = 1e3) # 22 minutes !!!! D:
summary(jmb.fit)

my.fit$coeffs$D 
summary(jmb.fit) %>% .$D
true.D

my.fit$coeffs$beta %>% t %>% c %>% matrix(., nc = 4, nr = 3, byrow = T)
c(summary(jmb.fit) %>% .$Outcome1 %>% .$Mean %>% .[1:4],
  summary(jmb.fit) %>% .$Outcome2 %>% .$Mean %>% .[1:4],
  summary(jmb.fit) %>% .$Outcome3 %>% .$Mean %>% .[1:4]) %>% matrix(., nc = 4, nr = 3, byr = T)
true.beta

c(my.fit$coeffs$eta, my.fit$coeffs$gamma)
summary(jmb.fit) %>% .$Survival %>% .$Mean %>% c
c(true.eta, true.gamma)
