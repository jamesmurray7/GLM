#' #####
#' Prereqs for negbinom/EM.R
#' Assumes WD set to ~/.../negbinom
#' #####

library(glmmTMB)
library(survival)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)

# Load cpp source files ---------------------------------------------------
message('Loading ', getwd(), '/ll.cpp')
sourceCpp('./ll.cpp')
message('Loading ', getwd(), '/../gammaCalc.cpp')
sourceCpp('../gammaCalc.cpp')
message('Loading ', getwd(), '/../lambdaUpdate.cpp')
sourceCpp('../lambdaUpdate.cpp')

# Set source on other prerequisites ---------------------------------------
source('./inits/inits.R') # initial conditions
source('./inits/MVLME.R')
source('../DataFunctions/longFns.R')
source('../DataFunctions/survFns.R')
vech <- function(x) x[lower.tri(x, diag = T)]

message('\n============================\n=== Prerequisites loaded ===\n============================')