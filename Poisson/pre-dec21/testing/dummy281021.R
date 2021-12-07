setwd('~/Documents/GLMM/Poisson/')
source('EM.R')
sourceCpp('../beta_test.cpp')
sourceCpp('../temp-gammaCalc.cpp')
source('../Simulations/simData.R')

beta <- rbind(c(0.75, -0.05, -0.2, 0.2),
              c(1, 0.05,  0.5, -0.80))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5

gamma <- c(-0.5, 0.8)
eta <- c(-0.3, 0.5)

data <- simData(250, 15, beta, D, gamma, eta, theta = c(-5, 0.2)) # appx. 50%
# Check all ids actually here, if not just rerun line above
length(unique(data$id))
ph <- coxph(Surv(survtime, status) ~ cont + bin, data = dplyr::distinct(data, id, cont, bin, survtime, status))

test2 <- em(data, ph, 3, nK = 2)
