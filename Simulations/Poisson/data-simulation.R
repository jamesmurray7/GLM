#' ###
#' Simulation Study
#' ###

source('./Simulations/simData.R')

#' ######
#' Simulate some data
#' ######
n <- 250; ntms = 15
beta <- rbind(c(0.75, -0.05, -0.2, 0.2),
              c(1, 0.05,  0.5, -0.80))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5

gamma <- c(-0.5, 1)
eta <- c(0.05, -0.3)
theta <- c(-6, 0.25)

data <- replicate(50, simData(n, ntms, beta = beta, D = D, gamma = gamma, eta = eta, theta = theta), simplify = F)

save(data, file = '~/Downloads/sim1.RData')


theta <- c(-3.75, 0.2)
ntms <- 10

data <- replicate(50, simData(n, ntms, beta = beta, D = D, gamma = gamma, eta = eta, theta = theta), simplify = F)

mean(do.call(c, lapply(data, function(x)
  sum(dplyr::distinct(x, id, status)$status)/n
)))

save(data, file = '~/Downloads/sim2.RData')

theta <- c(-2.6, 0.33)
ntms <- 6

data <- replicate(50, simData(n, ntms=ntms, beta = beta, D = D, gamma = gamma, eta = eta, theta = theta), simplify = F)

save(data, file = '~/Downloads/sim3.RData')

#' ####
#' 25 10 21
#' --
#' sim1: 30% failure rate on 15ntms; 
#' sim2: 40% failure rate on 10ntms.
#' sim3: 50% failure rate on 6ntms
#' ####
