## simulate some data ##
# set seed
library(glmmTMB)
set.seed(10)
# random effect
u <- rep(rnorm(500, 0, .1), each = 25)
# predictor
x <- rnorm(500 * 25, 0, .1)

# outcome
y <- rgamma(500 * 25, shape = exp(.5), scale = exp(1 + u + x)/exp(.5)) # log link...
#y <- rgamma(500 * 25, 1/(1 + u + x)) # inverse link

plot(y)

# id
id <- factor(rep(1:500, each = 25))

(m <- glmmTMB::glmmTMB(y ~ x + (1 | id), family = Gamma(link='log')))
summary(m)
exp(fixef(m)$dis)

# Do my implementations of log-likelihood match-up? -----------------------
library(Rcpp)
library(RcppArmadillo)
beta <- fixef(m)$cond
b <- ranef(m)$cond$id$`(Intercept)`
mu <- exp(x * beta + b)
shape <- exp(fixef(m)$disp)
glmmTMB:::logLik.glmmTMB(m)
dgamma_(y, shape, mu/shape)
