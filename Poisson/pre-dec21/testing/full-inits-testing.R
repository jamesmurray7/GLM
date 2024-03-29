library(MASS)
library(lme4)
library(nlme)
library(GLMMadaptive)
library(glmmTMB)
# Define a function -----
glm.fit <- function(data, pkg = ''){
  if(pkg == 'lme4'){
    res <- glmer(Y.1 ~ time + cont + bin + (1 + time|id),
                 data = data, 
                 family = 'poisson',
                 control = glmerControl(
                   optimizer = 'Nelder_Mead', boundary.tol = 1e-3,
                   optCtrl = list(FtolRel = 1e-3, XtolRel = 1e-3)
                 ))
  }else if(pkg == 'MASS'){
    res <- glmmPQL(fixed = Y.1 ~ time + cont + bin,
                  random = ~ 1 + time | id, family = poisson(link = 'log'), data = data, niter = 25,
                  verbose = F)
  }else if(pkg == 'GLMMadaptive'){
    res <- mixed_model(fixed = Y.1 ~ time + cont + bin, 
                              random = ~ 1 + time | id, data = data,
                              family = poisson())
  }else if(pkg == 'glmmTMB'){
    res <- glmmTMB(Y.1 ~ time + cont + bin + (1+time|id), data = data, family = poisson())
  }else{
    stop('Misspecified `pkg` argument')
  }
  res
}

# Load some data
load('~/Downloads/sim1.RData')
tmb <- glm.fit(data[[1]], 'glmmTMB')
lme <- glm.fit(data[[1]], 'lme4')

# Extracting VarCorr
matrix(glmmTMB::VarCorr(tmb)$cond$id,2,2)
matrix(lme4::VarCorr(lme)$id, 2, 2)





#' #####
#' Benchmarking
#' #####
microbenchmark::microbenchmark(
  `glmer` = {
    res1 <- glmer(Y.1 ~ time + cont + bin + (1 + time|id),
                 data = data[[1]], 
                 family = 'poisson',
                 control = glmerControl(
                   optimizer = 'Nelder_Mead', boundary.tol = 1e-3,
                   optCtrl = list(FtolRel = 1e-3, XtolRel = 1e-3)
                 ))
  },
  `glmmTMB` = {
    res2 <- glmmTMB(Y.1 ~ time + cont + bin + (1+time|id), data = data[[1]], family = poisson())
  },
  `glmmTMB-with-control` = {
    res3 <- glmmTMB(Y.1 ~ time + cont + bin + (1+time|id), data = data[[1]], family = poisson(),
                   control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))
  },
  times = 20
)
# glmmTMB one to beat.

# Try again with short profile
load('~/Downloads/sim3.RData')
data[[1]]

dd <- lapply(data, function(x){
  glmmTMB(Y.1 ~ time + cont + bin + (1+time|id), data = x, family = poisson(),
          control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))
})

dd2 <- lapply(data, function(x){
  glmmTMB(Y.1 ~ time + cont + bin + (1+time|id), data = x, family = poisson())
})

apply(do.call(rbind, lapply(dd, function(x) fixef(x)$cond)),2,mean)
apply(do.call(rbind, lapply(dd2, function(x) fixef(x)$cond)),2,mean)

#' ####
#' TO DO:
#' ----
#' Find that vignette and compare other methods.
#' ####

