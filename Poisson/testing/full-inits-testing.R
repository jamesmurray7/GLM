# glmer -----
library(lme4)
glmer.fit <- glmer(Y.1 ~ time + cont + bin + (1 + time|id),
                   data = data, 
                   family = 'poisson',
                   control = glmerControl(
                     optimizer = 'Nelder_Mead', boundary.tol = 1e-3,
                     optCtrl = list(FtolRel = 1e-3, XtolRel = 1e-3)
                   )
)

# MASS -----
library(MASS)
library(nlme)
pql.fit <- glmmPQL(fixed = Y.1 ~ time + cont + bin,
        random = ~ 1 + time | id, family = poisson(link = 'log'), data = data, niter = 25,
        control = lmeControl(tolerance = 1e-2), verbose = F)


# GLMMadaptive -----
library(GLMMadaptive)

microbenchmark::microbenchmark(
  `vanilla` = { # No Control arguments
    glmma.fit1 <- mixed_model(fixed = Y.1 ~ time + cont + bin, 
                              random = ~ 1 + time | id, data = data,
                              family = poisson())
  },
  `penalized` = { # Penalized around {0, 1}
    glmma.fit2 <- mixed_model(fixed = Y.1 ~ time + cont + bin, 
                              random = ~ 1 + time | id, data = data,
                              family = poisson(), penalized = T)
  },
  `tolx` = { # Changing tolerance criterion across fits
    glmma.fit3 <- mixed_model(fixed = Y.1 ~ time + cont + bin, 
                              random = ~ 1 + time | id, data = data,
                              family = poisson(),
                              control = list(tol3 = 1e-3))
  },
  `nAGQ=7` = {  # 4 less quadrature points
    glmma.fit4 <- mixed_model(fixed = Y.1 ~ time + cont + bin, 
                              random = ~ 1 + time | id, data = data,
                              family = poisson(),
                              control = list(nAGQ = 7))
  },
  times = 10
) -> bench

plot(bench)

# Benchmarking all --------------------------------------------------------

microbenchmark::microbenchmark(
  `lme4::glmer()` = {
    glmer.fit <- glmer(Y.1 ~ time + cont + bin + (1 + time|id),
                       data = data, 
                       family = 'poisson',
                       control = glmerControl(
                         optimizer = 'Nelder_Mead', boundary.tol = 1e-3,
                         optCtrl = list(FtolRel = 1e-3, XtolRel = 1e-3)
                       ))
  },
  `MASS::glmmPLQ()` = {
    pql.fit <- glmmPQL(fixed = Y.1 ~ time + cont + bin,
                       random = ~ 1 + time | id, family = poisson(link = 'log'), data = data, niter = 25,
                       control = lmeControl(tolerance = 1e-3), verbose = F)
  },
  `GLMMadaptive, nAGQ=7` = {
    glmma.fit4 <- mixed_model(fixed = Y.1 ~ time + cont + bin, 
                              random = ~ 1 + time | id, data = data,
                              family = poisson(),
                              control = list(nAGQ = 7))
  },
  times = 10 
) -> bench2

bench2

#' #####################
#' MASS::glmmPQL and lme4::glmer are comparable in speed
#' GLMMadaptive -- which exists for some reason -- is incomparably slow
#' #####################
