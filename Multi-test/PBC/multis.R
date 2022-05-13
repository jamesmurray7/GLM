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
  ascites ~ drug * time + (1 + time|id),
  spiders ~ drug * time + (1 + time|id),
  hepatomegaly ~ drug * time + (1 + time|id)
)

bin.fit1 <- EM(binoms, surv.formula, data, as.list(rep('binomial',3)),
               control = list(hessian = 'auto'))


# Count -------------------------------------------------------------------
# ((Univariate))
pois <- EM(list(platelets ~  drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id)),
           surv.formula, data, list(poisson), control = list(hessian='auto'))



# Full 8-variate model ----------------------------------------------------
long.formulas <- list(
  serBilir ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  SGOT ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  albumin ~ drug * time + (1 + time|id),
  prothrombin ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  platelets ~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
  spiders ~ drug * time + (1 + time|id),
  ascites ~ drug * time + (1 + time|id),
  hepatomegaly ~ drug * time + (1 + time|id)
)

families <- list(gaussian, gaussian, gaussian, gaussian, poisson, binomial, binomial, binomial)

fullfit <- EM(long.formulas, surv.formula, data, families, control = list(hessian = 'auto',
                                                                          correlated = F,
                                                                          maxit = 500))
fullfit2 <- EM(long.formulas, surv.formula, data, families, control = list(hessian = 'auto',
                                                                           correlated = T,
                                                                           maxit = 500))
# corr
Event-time sub-model: 
  Estimate    SE   2.5%  97.5% p-value
zeta_drug            -0.353 0.434 -1.203  0.497   0.416
gamma_serBilir        0.529 0.345 -0.147  1.204   0.125
gamma_SGOT            0.213 0.246 -0.268  0.695   0.386
gamma_albumin        -1.470 0.130 -1.724 -1.216   0.000
gamma_prothrombin    -2.257 0.157 -2.564 -1.950   0.000
gamma_platelets      -0.211 0.419 -1.033  0.611   0.615
gamma_spiders        -0.023 0.230 -0.474  0.428   0.921
gamma_ascites         0.299 0.202 -0.096  0.694   0.138
gamma_hepatomegaly    0.092 0.295 -0.485  0.669   0.755

