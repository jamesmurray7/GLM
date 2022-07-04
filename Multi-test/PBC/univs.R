#' ###########
#' univs.R
#' ---
#' Sequence of univariate fits for PBC data, with longit. specification determined by longs.R
#' serBilir         SGOT      albumin  prothrombin    platelets     alkaline      spiders      ascites 
#  "spline"     "spline"     "linear"     "spline"     "spline"  "quadratic"  "quadratic"     "spline" 
#  hepatomegaly 
#  "quadratic" 
#' ###########
rm(list=ls())
source('EM.R')
load('../PBC-case-study/PBC.RData')
pbc$serBilir <- log(pbc$serBilir)
pbc$SGOT <- log(pbc$SGOT)
pbc$prothrombin <- (.1*pbc$prothrombin)^(-4)
surv.formula <- Surv(survtime, status) ~ drug # Global

# define a function for fitting aEM object
EMfit <- function(response, type, family, gamma.SE = 'appx'){
  
  data.subset <- pbc[!is.na(pbc[,response]),]
  data.subset$id <- as.numeric(as.character(data.subset$id))
  ids <- unique(data.subset$id)
  if(any(diff(ids)>1)) stop('need to code unique ID assigment for', response)
  
  
  rhs <- switch(type,
                linear = ' ~ drug * time + (1 + time|id)',
                intercept = '~ drug * time + (1|id)',
                quadratic = '~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id)',
                spline = '~ drug * splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id)')
  long.formula <- list(as.formula(paste0(response, rhs)))
  
  fit <- EM(long.formula, surv.formula, data = data.subset, family = list(family), control = list(hessian = 'manual',
                                                                                                  gamma.SE = gamma.SE))
  print(my.summary(fit))
  return(invisible(fit))
}

# Serum bilirubin ---------------------------------------------------------
fit <- EM(list(serBilir ~ drug * splines::ns(time, knots = c(1,4)) + (1 + splines::ns(time, knots = c(1,4))|id)),
          surv.formula, pbc, family = list(gaussian), control = list(hessian = 'auto'))
my.summary(fit) # serBilir soundly sig. associated...

a <- EMfit('serBilir', 'spline', gaussian)

# SGOT --------------------------------------------------------------------
EMfit('SGOT', 'spline', gaussian)

# Albumin -----------------------------------------------------------------
EMfit('albumin', 'linear', gaussian, 'appx')

# Prothrombin time --------------------------------------------------------
EMfit('prothrombin', 'spline', gaussian)

# Platelets ---------------------------------------------------------------
EMfit('platelets', 'spline', poisson)

# Alkaline ----------------------------------------------------------------
EMfit('alkaline', 'quadratic', poisson) # p >> 0.1

# spiders -----------------------------------------------------------------
EMfit('spiders', 'quadratic', binomial)
EMfit('spiders', 'linear', binomial) # probably just take linear for parsimony
EMfit('spiders', 'intercept', binomial, 'exact')

# Ascites -----------------------------------------------------------------
EMfit('ascites', 'spline', binomial)
EMfit('ascites', 'quadratic', binomial)
EMfit('ascites', 'linear', binomial)
EMfit('ascites', 'intercept', binomial)


# Hepatomegaly ------------------------------------------------------------
EMfit('hepatomegaly', 'quadratic', binomial)
EMfit('hepatomegaly', 'linear', binomial)
EMfit('hepatomegaly', 'intercept', binomial)

