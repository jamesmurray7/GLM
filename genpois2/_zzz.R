long.formula = Y ~ time + cont + bin + (1|id)
surv.formula = Surv(survtime, status) ~ bin
disp.formula = ~ 1
control = list(verbose = T, delta.update.quad = T, beta.update.quad = F, include.all = T)
genpois.inits = T



#' temp <- function(delta, X, Z, b, Y, tau){
#' out <- numeric(length(w))
#' for(l in 1:length(w)){
#  mu <- exp(X %*% beta + Z %*% b + tau * v[l])
#   phi <- exp(delta)
#   out[l] <- w[l] * ll_genpois(mu, phi, Y)
# }
# sum(out)
# }