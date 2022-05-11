#' ###########
#' univs.R
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

fullfit <- EM(long.formula8, surv.formula, data, families, control = list(hessian = 'auto'))



long.formula <- list(spiders ~ drug * time + (1|id))
EM(long.formula, surv.formula ,data,list(binomial), control = list(hessian = 'auto'))
