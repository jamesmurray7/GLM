testfn <- function(beta, Y, X, Z, b){
  crossprod(Y, X %*% beta + Z %*% b) - exp(sum(X %*% beta + Z %*% b))
}

testfn2 <- function(beta, Y, X, Z, b){
  crossprod(Y, X %*% beta + Z %*% b) - sum(exp(X %*% beta + Z %*% b))
}

testfn(beta, Y[[2]], X[[2]], Z[[2]], b.hat[[2]])
testfn2(beta, Y[[1]], X[[1]], Z[[1]], b.hat[[1]])

numDeriv::grad(testfn, beta, method = 'Richardson', side = NULL, method.args = list(),
               Y = Y[[1]], X = X[[1]], Z = Z[[1]], b = b.hat[[1]])

numDeriv::hessian(testfn2, beta, method = 'Richardson', method.args = list(),
                  Y = Y[[1]], X = X[[1]], Z = Z[[1]], b = b.hat[[1]])

sourceCpp('./beta_test.cpp')

Sbeta.ND <- mapply(function(X, Z, Y, b){
  numDeriv::grad(beta_ll, beta, method = 'Richardson', side = NULL, method.args = list(),
                 Y = Y, X = X, Z = Z, b = b)
}, X = X, Z = Z, Y = Y, b = b.hat, SIMPLIFY = F)

Ibeta.ND <- mapply(function(X, Z, Y, b){
  -1 * numDeriv::hessian(beta_ll, beta, method = 'Richardson', method.args = list(),
                 Y = Y, X = X, Z = Z, b = b)
}, X = X, Z = Z, Y = Y, b = b.hat, SIMPLIFY = F)


# eta ---------------------------------------------------------------------
Fi.sumgb <- mapply(function(Fi, b){
  out <- 0;
  for(k in 1:nK) out <- out + gamma[k] * (Fi %*% b[[k]])
  out
}, Fi = Fi.list, b = b.hat.split, SIMPLIFY = F)

Fu.sumgb <- mapply(function(Fu, b){
  out <- 0;
  for(k in 1:nK) out <- out + gamma[k] * (Fu %*% b[[k]])
  out
}, Fu = Fu, b = b.hat.split, SIMPLIFY = F)

eta_ll <- function(eta, Di, K, Kr, Fi.sum, Fu.sum, l0u, tau.surv, w, v){
  lhs <- Di * (K %*% eta + Fi.sum)
  rhs <- numeric(length(l0u))
  for(l in 1:length(w)) rhs <- rhs + w[l] * exp(Kr %*% eta + Fu.sum + tau.surv * v[l])
  lhs - l0u %*% rhs
}

eta_ll(eta, Di[1], K[[1]], Krep[[1]], Fi.sumgb[[1]], Fu.sumgb[[1]], l0u[[1]], tau.surv[[1]], w, v)

Seta.ND <- numDeriv::grad(eta_ll, eta, method = 'Richardson', side = NULL, method.args = list(),
                          Di = Di[5], K = K[[5]], Kr = Krep[[5]], Fi.sum = Fi.sumgb[[5]], Fu.sum = Fu.sumgb[[5]], 
                          l0u = l0u[[5]], tau.surv = tau.surv[[5]], w=w, v=v)

Ieta.ND <- numDeriv::hessian(eta_ll, eta, method = 'Richardson',method.args = list(),
                          Di = Di[1], K = K[[1]], Kr = Krep[[1]], Fi.sum = Fi.sumgb[[1]], Fu.sum = Fu.sumgb[[1]], 
                          l0u = l0u[[1]], tau.surv = tau.surv[[1]], w=w, v=v)


gamma_ll <- function(gamma, Di, K, eta, Kr, Fi, Fu, l0u, bmat){
  gb <- as.numeric(gamma %*% bmat)
  Di * (K %*% eta + Fi %*% gb) - l0u %*% exp(Kr %*% eta + Fu %*% gb)
}
gamma_ll(gamma, Di[1], K[[1]], eta, Krep[[1]], Fi[1, ], Fu[[1]], l0u[[1]], temp[[1]])

numDeriv::grad(gamma_ll, gamma,  method = 'Richardson', side = NULL, method.args = list(),
               Di = Di[1], K = K[[1]], eta = eta, Kr = Krep[[1]], Fi = Fi[1, ], 
               Fu = Fu[[1]], l0u = l0u[[1]], bmat = temp[[1]])

numDeriv::hessian(gamma_ll, gamma,  method = 'Richardson', method.args = list(),
               Di = Di[1], K = K[[1]], eta = eta, Kr = Krep[[1]], Fi = Fi[1, ], 
               Fu = Fu[[1]], l0u = l0u[[1]], bmat = temp[[1]])



mu <- Fu[[1]] %*% as.numeric(gamma %*% temp[[1]]) 
dim(mu)
gb <- as.numeric(gamma %*% temp[[1]])


Di[1] * crossprod(temp[[1]], Fi[1,])

sourceCpp('../beta_test.cpp')
beta_ll(beta, X[[1]], Z[[1]], Y[[1]], b.hat[[1]])
beta_ll_quadrature(beta, X[[1]], Z[[1]], Y[[1]], b.hat[[1]], tau.long[[1]], w, v, 9)

numDeriv::grad(beta_ll, beta, method = 'Richardson', side = NULL, method.args = list(),
               Y = Y[[1]], X = X[[1]], Z = Z[[1]], b = b.hat[[1]])

numDeriv::grad(beta_ll_quadrature, beta, method = 'Richardson', side = NULL, method.args = list(),
               Y = Y[[1]], X = X[[1]], Z = Z[[1]], b = b.hat[[1]], tau = tau.long[[1]], w = w, v = v, gh = 9)
Sbeta.ND[[1]]

numDeriv::hessian(beta_ll_quadrature, beta, method = 'Richardson',method.args = list(),
                  Y = Y[[1]], X = X[[1]], Z = Z[[1]], b = b.hat[[1]], tau = tau.long[[1]], w = w, v = v, gh = 9)
