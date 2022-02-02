library(nlme)
library(splines)
library(tidyverse)
theme_set(theme_bw())
source('simData.R')

# Quadratic ---------------------------------------------------------------
beta <- rbind(c(1, -0.2, 0.01, 0.33, -0.50),
              c(0, -0.5, 0.05, -0.33, 0.50),
              c(3, 0.1, -0.05, 0.5, 0.1))

D <- as.matrix(Matrix::bdiag(replicate(3, diag(c(0.5^2, .2^2, .05^2)), simplify = F)))

data <- poly.simData(n = 250, ntms = 10, beta = beta, D = D, order = 2, theta0 = -5)

# Fit using raw poly
testpoly <- lme(Y.1 ~ poly(time, 2, raw = TRUE) + cont + bin, 
                random = ~ poly(time, 2, raw = TRUE)|id,
                data$dat,
                method = 'ML',
                control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

# Fit using bs(degree = 2)
testbs <- lme(Y.1 ~ bs(time, degree = 2) + cont + bin, 
                random = ~ bs(time, degree = 2)|id,
                data$dat,
                method = 'ML')

# Fit using ns(df = degree)
testns <- lme(Y.1 ~ ns(time, df = 2) + cont + bin, 
                random = ~ ns(time, df = 2)|id,
                data$dat,
                method = 'ML')

data$dat %>%
  select(id, time, Y.1:Y.3) %>%
  pivot_longer(Y.1:Y.3) %>%
  ggplot(aes(x = time, y = value, group = id)) +
  geom_line(alpha = .33) +
  facet_wrap(~name, scales = 'free') + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ poly(x, 2, raw = TRUE), colour = 'red', se = F) + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ bs(x, degree = 7), colour = 'green', lty = 3, se = F) + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ ns(x, df = 8), colour = 'blue', lty = 5, se = F)


# Cubic -------------------------------------------------------------------
rm(list=ls())
source('simData.R')
#               i    t     t^2   t^3    cont   bin
beta <- rbind(c(1, -0.2,  0.01, -0.01,  0.33, -0.50),
              c(0, -0.5,  0.05, -0.02, -0.33,  0.50),
              c(3,  0.1, -0.05,  0.02,  0.50,  0.1))

D <- as.matrix(Matrix::bdiag(replicate(3, diag(c(0.5^2, .2^2, .05^2, .015^2)), simplify = F)))

data <- poly.simData(n = 250, ntms = 10, beta = beta, D = D, order = 3, theta0 = -5)

# Fit using raw poly
testpoly <- lme(Y.1 ~ poly(time, 3, raw = TRUE) + cont + bin, 
                random = ~ poly(time, 3, raw = TRUE)|id,
                data$dat,
                method = 'ML',
                control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

# Fit using bs(degree = degree)
testbs <- lme(Y.1 ~ bs(time, degree = 3) + cont + bin, 
              random = ~ bs(time, degree = 3)|id,
              data$dat,
              method = 'ML',
              control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

# Fit using ns(df = degree)
testns <- lme(Y.1 ~ ns(time, df = 3) + cont + bin, 
              random = ~ ns(time, df = 3)|id,
              data$dat,
              method = 'ML',
              control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

data$dat %>%
  select(id, time, Y.1:Y.3) %>%
  pivot_longer(Y.1:Y.3) %>%
  ggplot(aes(x = time, y = value, group = id)) +
  geom_line(alpha = .33) +
  facet_wrap(~name, scales = 'free') + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ poly(x, 3, raw = TRUE), colour = 'red', se = F) + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ bs(x, degree = 3), colour = 'green', lty = 3, se = F) + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ ns(x, df = 3), colour = 'blue', lty = 5, se = F)


# Quartic -----------------------------------------------------------------
rm(list=ls())

source('simData.R')
#               i    t     t^2   t^3    t^4    cont   bin
beta <- rbind(c(1, -0.25,  0.05, -0.01,  0.025,  0.33, -0.50),
              c(0, -0.50,  0.10, -0.02, -0.030, -0.33,  0.50),
              c(3,  0.33, -0.15,  0.02,  0.010,  0.50,  0.1))

D <- as.matrix(Matrix::bdiag(replicate(3, diag(c(0.5^2, .2^2, .05^2, .015^2, 0.001^2)), simplify = F)))

data <- poly.simData(n = 250, ntms = 10, beta = beta, D = D, order = 4, theta0 = -5)

# Fit using raw poly
testpoly <- lme(Y.1 ~ poly(time, 4, raw = TRUE) + cont + bin, 
                random = ~ poly(time, 4, raw = TRUE)|id,
                data$dat,
                method = 'ML',
                control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

# Fit using bs(degree = degree)
testbs <- lme(Y.1 ~ bs(time, degree = 4) + cont + bin, 
              random = ~ bs(time, degree = 4)|id,
              data$dat,
              method = 'ML',
              control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

# Fit using ns(df = degree)
testns <- lme(Y.1 ~ ns(time, df = 4) + cont + bin, 
              random = ~ ns(time, df = 4)|id,
              data$dat,
              method = 'ML',
              control = nlme::lmeControl(opt = "optim", msTol = 1e-3))

data$dat %>%
  select(id, time, Y.1:Y.3) %>%
  pivot_longer(Y.1:Y.3) %>%
  ggplot(aes(x = time, y = value, group = id)) +
  geom_line(alpha = .33) +
  facet_wrap(~name, scales = 'free') + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ poly(x, 4, raw = TRUE), colour = 'red', se = F) + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ bs(x, degree = 4), colour = 'green', lty = 3, se = F) + 
  geom_smooth(aes(group = NULL),
              method = 'lm', formula = y ~ ns(x, df = 4), colour = 'skyblue', lty = 5, se = F)

