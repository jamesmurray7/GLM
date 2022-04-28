rm(list=ls())
source('EM.R')
load('../PBC-case-study/PBC.RData')
pbc$serBilir <- log(pbc$serBilir)
long.formulas <- list(
  serBilir ~ drug * time + (1 + time|id),
  albumin ~ drug * time + (1 + time|id)
)
family <- list(gaussian, gaussian)
surv.formula <- Surv(survtime, status) ~ drug

my.fit <- EM(long.formulas, surv.formula, pbc, family, control = list(hessian='auto'))

mjoint(
  formLongFixed = list(
    '1' = serBilir ~ drug * time,
    '2' = albumin ~ drug * time
  ),
  formLongRandom = list(
    '1' = ~ time|id,
    '2' = ~ time|id
  ),
  formSurv = surv.formula, timeVar = 'time',
  data = pbc, control = list(type = 'sobol', convCrit = 'rel', tol2 = 1e-2)
) ->mj

mj$comp.time

my.fit2$coeffs
summary(mj)


# Some nutty longitudinal formula -----------------------------------------

long.formulas <- list(
  serBilir ~ drug * splines::ns(time, knots = c(1, 4)) + 
    (1 + splines::ns(time, knots = c(1, 4)) | id),
  albumin ~ drug * splines::ns(time, knots = c(1, 4)) + 
    (1 + splines::ns(time, knots = c(1, 4)) | id)
)

my.fit.spline <- EM(long.formulas, surv.formula, pbc, family, control = list(verbose=T, hessian='auto'))

mjoint(
  formLongFixed = list(
    '1' = serBilir ~ drug * splines::ns(time, knots = c(1, 4)),
    '2' = albumin ~ drug * splines::ns(time, knots = c(1, 4))
  ),
  formLongRandom = list(
    '1' = ~ 1 + splines::ns(time, knots = c(1, 4)) | id,
    '2' = ~ 1 + splines::ns(time, knots = c(1, 4)) | id
  ),
  formSurv = surv.formula, timeVar = 'time',
  data = pbc, control = list(type = 'sobol', convCrit = 'rel', tol2 = 1e-2)
) ->mj.spline


# Trivariate mixture ------------------------------------------------------

pbc.tri <- pbc[!is.na(pbc$serBilir) & !is.na(pbc$albumin) & !is.na(pbc$platelets), ]
# Two Gaussian, one poisson (Rustand formulation for each is done later)...
long.formulas <- list(
  serBilir ~  drug * time + (1 + time | id),
  albumin ~ drug * time + (1 + time | id),
  platelets ~ drug * time + (1 + time | id)
)
family <- list(gaussian, gaussian, poisson)

my.fit.tri1 <- EM(long.formulas, surv.formula, pbc.tri, family, control = list(verbose = T, hessian = 'auto'))

# Rustand formulation for each.
long.formulas <- list(
  serBilir ~  drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1,4)) | id),
  albumin ~ drug * time + (1 + time | id),
  platelets ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1,4)) | id)
)
my.fit.tri.spline <- EM(long.formulas, surv.formula, pbc.tri, family, control = list(verbose = T, hessian = 'manual', optimiser= 'ucminf'))
my.fit.tri.spline2 <- EM(long.formulas, surv.formula, pbc.tri, family, control = list(verbose = T, hessian = 'auto', optimiser= 'optim'))



# Five-variate mixture as described in Rustand et al. ---------------------
library(splines)
pbc.quin <- pbc[!is.na(pbc$serBilir) & !is.na(pbc$SGOT) & !is.na(pbc$albumin) &
                    !is.na(pbc$platelets) & !is.na(pbc$spiders),]
pbc.quin$SGOT <- log(pbc.quin$SGOT)
any(diff(unique(as.numeric(pbc.quin$id))) > 1) # check continuous ids

long.formulas <- list(
  serBilir ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1,4)) | id),
  SGOT ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1,4)) | id),
  albumin ~ drug * time + (1 + time|id),
  platelets ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1,4)) | id),
  spiders ~ drug * time + (1 + time|id)
)

family <- list(gaussian, gaussian, gaussian, poisson, binomial)

my.rustand <- EM(long.formulas, surv.formula, pbc.quin, family,
                 control = list(hessian = 'auto'))
source('multivariate-compare.R')
# load('RUSTANDFIT.RData')
my.summary(my.rustand)


# JMbayes2 Rustand fit ----------------------------------------------------
library(JMbayes2)
M1 <- lme(serBilir ~ drug * ns(time, knots = c(1,4)), random = ~ns(time, knots = c(1, 4))|id, data = pbc.quin,
          method = "ML",
          control = nlme::lmeControl(opt = "optim", msTol = 1e-3))
M2 <- lme(SGOT ~ drug * ns(time, knots = c(1,4)), random = ~ns(time, knots = c(1, 4))|id, data = pbc.quin,
          method = "ML",
          control = nlme::lmeControl(opt = "optim", msTol = 1e-3))
M3 <- lme(albumin ~ drug * time, random = ~time|id, data = pbc.quin)
temp <- glmmTMB(platelets ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1,4)) | id),
                data = pbc.quin, family = poisson)
M4 <- mixed_model(platelets ~ drug * ns(time, knots = c(1, 4)), random = ~ns(time, knots = c(1, 4)) | id,
                  data = pbc.quin, family = poisson,
                  initial_values = list(
                    betas = glmmTMB::fixef(temp)$cond,
                    D = matrix(glmmTMB::VarCorr(temp)$cond$id, 4, 4)
                  ))
M5 <- mixed_model(spiders ~ drug * time, random = ~ time|id, data = pbc.quin, family = binomial)
Ms <- list(M1, M2, M3, M4, M5)
survdata <- parseCoxph(surv.formula, pbc.quin)$survdata
phfit <- coxph(Surv(survtime, status) ~ drug, survdata)

JMbayes2.rustand <- jm(
  Surv_object = phfit,
  Mixed_objects = Ms,
  time_var = 'time',
  data_Surv = survdata
)




