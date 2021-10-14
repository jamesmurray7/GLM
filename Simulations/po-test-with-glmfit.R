#' ############
#' Short program which simulates a poisson model and then fits
#' using lme4::glmer and GLMMadaptive::mixed_model
#' ############
rm(list=ls())

# Simulate a Poisson Model ----
n <- 100
D <- matrix(c(0.5, 0, 0, 0.1), 2, 2)
b <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = D)
time <- 0:5
cont <- rnorm(n)
bin <- rbinom(n, 1, .5)
beta <- c(3.25, 0.15, -0.15, 0.6)

dfs <- list(); lambda_is <- numeric(n)
for(i in 1:n){
    Yi <- numeric(6)
    lambda_i <- cbind(1, time, cont[i], bin[i]) %*% beta + cbind(1, time) %*% b[i, ]
    for(j in seq_along(lambda_i)){
      Yi[j] <- rpois(1, lambda_i[j])
    }
    dfs[[i]] <- data.frame(id = i, time = time, cont = cont[i], bin = bin[i], Y = Yi)
}

df <- do.call(rbind, dfs)

# Lets compare some model fits ----
# lme4
library(lme4)
glm1 <- glmer(Y ~ time + cont + bin + (1 + time | id), data = df, 
              family = 'poisson')
library(GLMMadaptive)
glm2 <- mixed_model(fixed = Y ~ time + cont + bin,
                    random = ~ 1 + time | id,
                    family = poisson(),
                    data = df)

summary(glm1)$coef
glm2$coefficients
glm2$D
