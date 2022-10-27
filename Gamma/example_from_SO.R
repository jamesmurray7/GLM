## simulate some data ##
# set seed
set.seed(10)
# random effect
u <- rep(rnorm(500, 0, .1), each = 25)
# predictor
x <- rnorm(500 * 25, 0, .1)

# outcome
y <- rgamma(500 * 25, shape = 5.00, scale = exp(1 + u + x)/5.00) # log link...

#y <- rgamma(500 * 25, 1/(1 + u + x)) # inverse link
plot(y)

# id
id <- factor(rep(1:500, each = 25))

(m <- glmmTMB::glmmTMB(y ~ x + (1 | id), family = Gamma(link='log')))
summary(m)
exp(fixef(m)$dis)
