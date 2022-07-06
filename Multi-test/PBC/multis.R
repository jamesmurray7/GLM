#' ###########
#' multis.R
#' ---
#' Multivariate + backwards elimination for PBC, starting with:
#' Continuous //
#'    serBilir, SGOT, albumin, prothrombin
#' Count //
#'    platelets, 
#' Binomial //
#'    spiders, ascites, hepatomegaly.
#' ###########
rm(list=ls())
source('EM.R')
load('../PBC-case-study/PBC.RData')
data <- na.omit(pbc[, c("id", "survtime", "status","drug","age", 
                        "sex","time","serBilir","SGOT", "albumin", 'prothrombin',
                        "platelets", "ascites")])
data$id <- as.numeric(as.character(data$id))

uid <- unique(data$id)
diff(uid) # check no ID differences > 1, or program wont work.

data$serBilir <- log(data$serBilir)
data$SGOT <- log(data$SGOT)
data$prothrombin <- (.1*data$prothrombin)^(-4)
surv.formula <- Surv(survtime, status) ~ drug # Global


# All gaussians -----------------------------------------------------------
gaussians <- list(
  serBilir ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  SGOT ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  albumin ~ drug * time + (1 + time|id),
  prothrombin ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id)
)
families <- as.list(rep('gaussian', 4))
gauss.fit1 <- EM(gaussians, surv.formula, data, families, control = list(hessian = 'auto'))
gauss.fit2 <- EM(gaussians, surv.formula, data, families, control = list(hessian = 'auto', correlated = F))

# 
# Event-time sub-model: CORRELATED
#   Estimate    SE   2.5%  97.5% p-value
# zeta_drug           -0.193 0.382 -0.942  0.556   0.613
# gamma_serBilir       0.813 0.278  0.269  1.358   0.003
# gamma_SGOT          -0.040 0.557 -1.132  1.051   0.942
# gamma_albumin       -2.020 0.778 -3.545 -0.494   0.009
# gamma_prothrombin   -3.247 1.431 -6.052 -0.442   0.023

# Count -------------------------------------------------------------------
# ((Univariate))
pois <- EM(list(platelets ~  drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id)),
           surv.formula, data, list(poisson), control = list(SEs = 'appx'))

# Rustand -----------------------------------------------------------------
long.formulas <- list(
  serBilir ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  SGOT ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  albumin ~ drug * time + (1 + time|id),
  platelets ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  spiders ~ drug * time + (1 + time|id)
)

rustand.fit <- EM(
  long.formulas, surv.formula, data, list('gaussian', 'gaussian', 'gaussian', 'poisson', 'binomial'),
  control = list(SEs = 'appx')
)

rustand.fit2 <- EM(
  long.formulas, surv.formula, data, list('gaussian', 'gaussian', 'gaussian', 'poisson', 'binomial'),
  control = list(SEs = 'exact')
)

# Full 6-variate model ----------------------------------------------------
long.formulas <- list(
  serBilir ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  SGOT ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  albumin ~ drug * time + (1 + time|id),
  prothrombin ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  platelets ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  # spiders ~ drug * time + (1|id), # Remove these --> makes it take longer and doesn't change outcome.
  ascites ~ drug * time + (1+time|id)
  # hepatomegaly ~ drug * time + (1|id)
)

long.formulas <- list(
  serBilir ~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id),
  SGOT ~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id),
  albumin ~ drug * time + (1 + time|id),
  prothrombin ~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id),
  platelets ~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id),
  # spiders ~ drug * time + (1|id), # Remove these --> makes it take longer and doesn't change outcome.
  ascites ~ drug * time + (1|id)
  # hepatomegaly ~ drug * time + (1|id)
)

long.formulas <- list(
  serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
  SGOT ~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id),
  albumin ~ drug * time + (1 + time|id),
  prothrombin ~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id),
  platelets ~ drug * time + (1 + time|id),
  # spiders ~ drug * time + (1|id), # Remove these --> makes it take longer and doesn't change outcome.
  ascites ~ drug * time + (1|id)
  # hepatomegaly ~ drug * time + (1|id)
)


families <- list(gaussian, gaussian, gaussian, gaussian, poisson, binomial)#, binomial, binomial)

fullfit <- EM(long.formulas, surv.formula, data, families, control = list(correlated = F,
                                                                          SEs = 'exact',
                                                                          verbose = T,
                                                                          maxit = 500))

fullfit2 <- EM(long.formulas, surv.formula, data, families, control = list(correlated = T,
                                                                           SEs = 'exact',
                                                                           verbose = T,
                                                                           maxit = 500))
fit2xtab(fullfit, 15)
fit2xtab(fullfit2, 15)
# Taking trivariate forward: ----------------------------------------------
tri.formula <- list(
  serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
  albumin ~ drug * time + (1 + time|id),
  prothrombin ~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id)
)

trifam <- list('gaussian', 'gaussian', 'gaussian')
trifit <- EM(tri.formula, surv.formula, data, trifam, control = list(correlated = F,
                                                                          SEs = 'exact',
                                                                          verbose = F,
                                                                          maxit = 500))


my.summary(trifit)
fit2xtab(trifit)
trifit2 <- EM(tri.formula, surv.formula, data, trifam, control = list(correlated = T,
                                                                       SEs = 'exact',
                                                                       verbose = T,
                                                                       maxit = 500))
my.summary(trifit2)
fit2xtab(trifit2)
