#' ###
#' Simulation Study
#' ###

source('./Simulations/simData.R')
sourceCpp('./beta_test.cpp')
source('EM.R')

#' ######
#' Simulate some data
#' ######
n <- 250; ntms = 10
beta <- rbind(c(0.75, -0.05, -0.2, 0.2),
              c(1, 0.05,  0.5, -0.80))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5

gamma <- c(-0.25, 0.6)
eta <- c(-0.1, 0.25)
theta <- c(-5, 0.18)

num.sets <- 0
full.data <- list(); p <- 1
while(num.sets < 50){
  data <- simData(n, ntms, beta = beta, D = D, gamma = gamma, eta = eta, theta = theta)
  if(length(unique(data$id)) == n){
    num.sets <- num.sets + 1
    full.data[[p]] <- data
    p <- p + 1
  }
}

save(full.data, file = '~/Downloads/sim1.RData')

# Change theta
theta <- c(-3, 0.2)
num.sets <- 0
full.data <- list(); p <- 1
while(num.sets < 50){
  data <- simData(n, ntms, beta = beta, D = D, gamma = gamma, eta = eta, theta = theta)
  if(length(unique(data$id)) == n){
    num.sets <- num.sets + 1
    full.data[[p]] <- data
    p <- p + 1
  }
}

save(full.data, file = '~/Downloads/sim2.RData')
