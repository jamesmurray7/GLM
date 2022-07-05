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
data <- pbc[!is.na(pbc$serBilir) & !is.na(pbc$SGOT) & !is.na(pbc$albumin) & !is.na(pbc$platelets) & !is.na(pbc$spiders) & !is.na(pbc$ascites) & !is.na(pbc$hepatomegaly),]
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

# All binary --------------------------------------------------------------
binoms <- list(
  ascites ~ drug * time + (1|id),
  spiders ~ drug * time + (1|id),
  hepatomegaly ~ drug * time + (1|id)
)

bin.fit1 <- EM(binoms, surv.formula, data, as.list(rep('binomial',3)),
               control = list(hessian = 'manual', gamma.SE = 'exact'))


# Count -------------------------------------------------------------------
# ((Univariate))
pois <- EM(list(platelets ~  drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id)),
           surv.formula, data, list(poisson), control = list(SE.D = F, SEs = 'score'))



# Full 8-variate model ----------------------------------------------------
long.formulas <- list(
  serBilir ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  SGOT ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  albumin ~ drug * time + (1 + time|id),
  prothrombin ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  platelets ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  # spiders ~ drug * time + (1|id), # Remove these --> makes it take longer and doesn't change outcome.
  ascites ~ drug * time + (1|id)
  # hepatomegaly ~ drug * time + (1|id)
)

families <- list(gaussian, gaussian, gaussian, gaussian, poisson, binomial)#, binomial, binomial)

fullfit <- EM(long.formulas, surv.formula, data, families, control = list(correlated = F,
                                                                          SEs = 'exact.gamma',
                                                                          SE.D = F,
                                                                          verbose = T,
                                                                          maxit = 500))

fullfit2 <- EM(long.formulas, surv.formula, data, families, control = list(correlated = T,
                                                                           gamma.SE = 'exact',
                                                                           verbose = T,
                                                                           maxit = 500))
# corr
my.summary(fullfit2)
Event-time sub-model: 
  Estimate    SE   2.5% 97.5% p-value
zeta_drug           -0.302 0.677 -1.629 1.024   0.655
gamma_serBilir       0.640 0.177  0.293 0.987   0.000
gamma_SGOT           0.180 0.359 -0.525 0.884   0.617
gamma_albumin       -1.750 0.925 -3.564 0.064   0.059
gamma_prothrombin   -2.025 1.645 -5.249 1.200   0.218
gamma_platelets     -0.224 0.256 -0.726 0.277   0.381
gamma_ascites        0.304 0.171 -0.031 0.639   0.075

# Non-corr
my.summary(fullfit)
Event-time sub-model: 
                  Estimate    SE   2.5%  97.5% p-value
zeta_drug           -0.192 0.255 -0.692  0.308   0.451
gamma_serBilir       1.136 0.114  0.912  1.360   0.000
gamma_SGOT          -0.436 0.265 -0.955  0.083   0.100
gamma_albumin       -1.868 0.324 -2.502 -1.234   0.000
gamma_prothrombin   -0.783 0.653 -2.064  0.497   0.231
gamma_platelets     -0.458 0.205 -0.859 -0.056   0.025
gamma_ascites        0.091 0.075 -0.057  0.238   0.228
