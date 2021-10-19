library(survival)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)

rm(list=ls())   
vech <- function(x) x[lower.tri(x, diag = T)]
sourceCpp('Poisson/ll.cpp')
source('Simulations/simData.R')
source('./DataFunctions/longFns.R')
source('./DataFunctions/survFns.R')

# Simulate some Data
args(simData)
beta <- rbind(c(0.75, -0.05, -0.2, 0.2),
              c(1, 0.05,  0.5, -0.80))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5

gamma <- c(-0.5, 0.8)
eta <- c(-0.3, 0.5)

data <- simData(250, 10, beta, D, gamma, eta, theta = c(-3, 0.2)) # appx. 50%
# Check all ids actually here, if not just rerun line above
length(unique(data$id))
ph <- coxph(Surv(survtime, status) ~ cont + bin, data = dplyr::distinct(data, id, cont, bin, survtime, status))

X <- getXi(data, 2)
Y <- getYi(data, 2)
Z <- getZi(data, 2)
mi <- getmi(data, 2) # Dont actually think this will be needed(!)

#' Initial conditions -----
source('inits/inits.R')
inits.long <- Longit.inits(2, data)
inits.surv <- TimeVarCox(data, Ranefs(inits.long))

sv <- surv.mod(ph, data, inits.surv$l0.init)
# Extract all survival-related objects
ft <- sv$ft; nev <- sv$nev
surv.ids <- sv$surv.ids; surv.times <- sv$surv.times
Di <- sv$Di
l0 <- sv$l0; l0i <- sv$l0i; l0u <- sv$l0u
Fi <- sv$Fi; Fu <- sv$Fu
K <- getKi(data, ph)

#' MVLME Step ----
source('./inits/MVLME.R')
mvlme.fit <- mvlme(data, Y, X, Z, inits.long, nK = 2, q = 4,  
                   mvlme.tol = 5e-3)

beta <- mvlme.fit$beta.init
D <- mvlme.fit$D.init
b <- mvlme.fit$b # currently a list
gamma <- inits.surv$inits[3:4]
eta <- inits.surv$inits[1:2]


gh <- statmod::gauss.quad.prob(9, 'normal')
w <- gh$w; v <- gh$n


