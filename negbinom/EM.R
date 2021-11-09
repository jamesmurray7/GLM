#' ########
#' EM.R ----
#' nK Poisson-Gamma deviate responses (i.e. negative binomial) with a survival sub-model.
#' Assumes wd is ~/.../GLMM/negbinom
#' ########

source('./prerequisites.R')

em <- function(data, ph, gh.nodes, collect.hist = T, max.iter = 200, 
               tol = 0.01, diff.type = "abs.rel",
               b.tol = .1, numOptims = 4, nK = 3, post.process = T, verbose = F){
  start.time <- proc.time()[3]
  #' Algorithm-specific ----
  q <- nK * 2
  b.diff <- diff <- 100; iter <- 0; 
  
  #' Data Matrices ----
  X <- getXi(data, nK); Y <- getYi(data, nK); Z <- getZi(data, nK)
  K <- getKi(data, ph)
  n <- length(Y)
  
  #' Initial Conditions ----
  inits.long <- Longit.inits(nK, data)
  inits.surv <- TimeVarCox(data, Ranefs(inits.long))
  # MVLME for more optimal starting values in longit part
  mvlme.fit <- mvlme(data, Y, X, Z, inits.long, nK = nK, q = q,  
                     mvlme.tol = 5e-3, verbose = verbose)
  
  #' Survival-related objects ----
  sv <- surv.mod(ph, data, inits.surv$l0.init)
  ft <- sv$ft; nev <- sv$nev
  surv.ids <- sv$surv.ids; surv.times <- sv$surv.times
  Di <- sv$Di
  l0 <- sv$l0; l0i <- sv$l0i; l0u <- sv$l0u
  Fi <- sv$Fi; Fu <- sv$Fu
  # List-versions of several data objects for use of mapply()
  Fi.list <- lapply(1:nrow(Fi), function(i) Fi[i, ])
  rvFi.list <- lapply(1:nrow(Fi), function(i) do.call(c, replicate(nK, Fi[i,], simplify = F)))
  Krep <- sapply(1:n, function(x){  # For updates to \eta
    x <- apply(K[[x]], 2, rep, nrow(Fu[[x]]))
    if("numeric" %in% class(x)) x <- t(as.matrix(x))
    x
  }) 
  
  #' Collect all inits ----
  beta <- mvlme.fit$beta
  D <- mvlme.fit$D
  b <- mvlme.fit$b 
  theta <- mvlme.fit$theta
  gamma <- inits.surv$inits[3:length(inits.surv$inits)]; gr <- rep(gamma, each = 2)   # 2 for intercept/slope
  eta <- inits.surv$inits[1:2]
  theta.list <- lapply(1:n, function(i){ # list version of theta.
    out <- c();
    for(k in 1:nK) out <- c(out, rep(theta[k], nrow(X[[i]])/2))
    out
  })
  
  # And cast to parameter vector
  params <- c(vech(D), beta, theta, gamma, eta)
  names(params) <- c(rep('D', length(vech(D))), # vech(D)
                     names(beta), names(theta),
                     names(gamma), 
                     paste0('eta_', names(eta))
                     )
  
  #' Quadrature ----
  gh <- statmod::gauss.quad.prob(gh.nodes, 'normal')
  w <- gh$w; v <- gh$n
  
  #' Collect data matrices for post-processing and the iteration history
  if(post.process) dmats <- list(Y = Y, X = X, Z = Z, K = K, sv=sv)
  if(collect.hist) iter.hist <- data.frame(iter = iter, t(params))
  inds <- split(seq(nK * 2), rep(1:nK, each = 2))
  
  #' Start EM ----
  EMtime <- proc.time()[3]
  message('Starting EM Algorithm')
  while(diff > tol & iter <= max.iter){
     start.estep <- proc.time()[3]
     #' #################
     #' E-step =========#
     #' #################
     
     #' bhat and \Sigma_i ---------
     b.hat <- mapply(function(b, Y, lfactY, X, Z, theta, K, Delta, l0i, Fi, l0u, Fu, rvFi){
       ucminf::ucminf(
         b, ll, gradll, Y, lfactY, X, Z, D, theta, K, Delta, l0i, Fi, l0u, Fu, 
         beta, eta, gr, rvFi, nK, q
       )$par}, 
       b = b, Y = Y, lfactY = lapply(Y, lfactorial), X = X, Z = Z, theta = theta.list, K = K, Delta = as.list(Di),
       l0i = as.list(l0i), Fi = Fi.list, l0u = l0u, Fu = Fu, rvFi = rvFi.list, SIMPLIFY = F)
     b.hat.split <- lapply(b.hat, function(y) lapply(inds, function(x) y[x]))
     
     Sigmai <- mapply(function(b, X, Z, Y, theta, K, l0u, Fu){
       solve(-1 * hessll(b, X, Z, Y, D, theta, K, l0u, Fu, beta, eta, gr, nK))
     }, b = b.hat, Y = Y, X = X, Z = Z, K = K, theta=theta.list,
     l0u = l0u, Fu = Fu, SIMPLIFY = F)
     S <- lapply(Sigmai, function(y) lapply(inds, function(x) y[x,x]))
     
     #' Steps to update Longitudinal and RE parameters --------
     Sbeta <- mapply(function(b, Y, X, Z, theta){
       crossprod(X, Y) - crossprod(X, (Y + theta) * exp(X %*% beta + Z %*% b) / (exp(X %*% beta + Z %*% b) + theta))
     }, b = b.hat, Y = Y, X = X, Z = Z, theta = theta.list, SIMPLIFY = F)
     
     Ibeta <- mapply(function(b, X, Z, Y, theta){ # NB this -I(\beta)
       mu <- exp(X %*% beta + Z %*% b) # This may well be slower than doing an 'outside' mapply for my first.
       term1 <- (Y + theta) * mu
       term2 <- mu + theta
       -1 * (crossprod(diag(x = c(term1 / term2)) %*% X, X) - crossprod(diag(x = c(term1 * mu / (term2 * term2))) %*% X, X))
     }, b = b.hat, X = X, Z = Z, Y = Y, theta = theta.list, SIMPLIFY = F)
     
     #' --- Steps to update theta ----
     #'           ((NYI))
     #' ------------------------------
     
     #' Steps to update D -----------
     D.news <- mapply(function(S, b){
       S + tcrossprod(b)
     }, S = Sigmai, b = b.hat, SIMPLIFY = F)
     
     #' Steps to update survival parameters ------------
     # Define mu.surv
     mu.surv <- mapply(function(K, Fu, b){
       rhs <- 0
       for(k in 1:nK) rhs <- rhs + gamma[k] * b[inds[[k]]]
       exp(K %*% eta + Fu %*% rhs)
     }, K = Krep, Fu = Fu, b = b.hat, SIMPLIFY = F)
     
     # Define tau objects
     tau <- mapply(function(Fu, S){
       out <- numeric(nrow(Fu))
       for(k in 1:nK) out <- out + diag(gamma[k]^2 * tcrossprod(Fu %*% S[[k]], Fu))
       out
     }, Fu = Fu, S = S, SIMPLIFY = F)
     
     tau.tilde <- mapply(function(Fu, S){
       mat <- matrix(0, nr = nK, nc = nrow(Fu))
       for(k in 1:nK) mat[k, ] <- diag(tcrossprod(Fu %*% S[[k]], Fu))
       mat
     }, Fu = Fu, S = S, SIMPLIFY = F)
     
     tau.surv <- lapply(tau, sqrt)
     tau2.surv <- lapply(tau, function(x){
       x <- x^(-0.5)
       if(any(is.infinite(x))) x[which(is.infinite(x))] <- 0   # Avoid NaN
       x
     })
     
     #' S(\gamma) ----
     Sgamma <- mapply(function(Delta, tau.surv, mu.surv, l0u, Fu, Fi, b){
       t(Sgammacalc(gamma, Delta, tau.surv, mu.surv, l0u, Fu, Fi, w, v, b, nK, gh.nodes))
     }, Delta = as.list(Di), tau.surv = tau.surv, mu.surv = mu.surv, l0u = l0u,
     Fu = Fu, Fi = Fi.list, b = b.hat.split, SIMPLIFY = F)
     
     #' I(\gamma)
     Igamma <- mapply(function(tau.tilde, tau.surv, tau2.surv, mu.surv, Fu, l0u, b){
       gamma2Calc(gamma, tau.tilde, tau.surv, tau2.surv, mu.surv, w, v, Fu, l0u, b, nK, gh.nodes)
     }, tau.tilde = tau.tilde, tau.surv = tau.surv, tau2.surv = tau2.surv,
     mu.surv = mu.surv, Fu = Fu, l0u = l0u, b = b.hat.split, SIMPLIFY = F)
     
     #' S(\eta) ----
     Seta <- mapply(function(K, KK, Delta, l0u, mu.surv, tau.surv){
       cen <- Delta %*% K
       rhs <- c(0, 0)
       for(l in 1:gh.nodes) rhs <- rhs + w[l] * t(KK) %*% (l0u * (mu.surv * exp(tau.surv * v[l])))
       cen-t(rhs)
     }, K = K, KK = Krep, Delta = as.list(Di), l0u = l0u, mu.surv = mu.surv, tau.surv = tau.surv, SIMPLIFY = F)
     
     #' I(\eta) ----
     Ieta <- mapply(function(K, KK, tau.surv, mu.surv, l0u){
       Ietacalc(2, K, KK, tau.surv, mu.surv, l0u, w, v, gh.nodes)
     }, K = K, KK = Krep, tau.surv = tau.surv, mu.surv = mu.surv, l0u = l0u, SIMPLIFY = F)
     
     #' Second derivates of \gamma and \eta ('cross-terms') -----
     Igammaeta <- list() # for some reason I can't get mapply() to work with this 
     # Old
     for(i in 1:n){
       Igammaeta[[i]] <- Igammaetacalc(2, Krep[[i]], tau.surv[[i]], tau.tilde[[i]],
                                       mu.surv[[i]], l0u[[i]], Fu[[i]], b.hat.split[[i]], gamma, w, v, nK, gh.nodes)
     }

     #' #################
     #' M-step =========#
     #' #################
     
     #' \beta -----
     beta.new <- beta + solve(-Reduce('+', Ibeta), 
                              rowSums(do.call(cbind, Sbeta)))
     
     #' D -----
     D.new <- Reduce('+', D.news)/n
     
     #' The baseline hazard, \lambda_0
     lambda <- lambdaUpdate(surv.times, ft, gamma, eta, K, S,
                            b.hat.split, n, w, v, gh.nodes)
     
     l0.new <- nev/rowSums(lambda)
     l0u.new <- lapply(l0u, function(x){
       ll <- length(x); l0.new[1:ll]
     })
     l0i.new <- c()
     l0i.new[which(Di == 0)] <- 0 
     l0i.new[which(Di == 1)] <- l0.new[match(Fi[which(Di==1), 2], ft)]
     
     #' (\gamma, \eta) ----
     # Set up score vector and information matrix
     Sge <- c(colSums(do.call(rbind, Sgamma)), colSums(do.call(rbind, Seta)))
     Imat <- as.matrix(Matrix::bdiag(Reduce('+', Igamma),
                                     Reduce('+', Ieta)))
     # Fill in off-block diagonal with Igammaeta
     eta.inds <- (nK+1):ncol(Imat)
     for(k in 1:nK){
       Imat[k, eta.inds] <- Imat[eta.inds, k] <- rowSums(do.call(cbind, lapply(Igammaeta, '[[', k)))
     }
     
     gamma.eta.new <- c(gamma, eta) + solve(Imat, Sge)
     
     gamma.new <- gamma.eta.new[1:nK]
     eta.new <- gamma.eta.new[(nK + 1):length(gamma.eta.new)]
     
     EMtime[iter + 1] <- proc.time()[3] - start.estep # M step finishes here ---
     
     #' Update parameters and print ----
     params.new <- c(vech(D.new), beta.new, theta, gamma.new, eta.new); names(params.new) <- names(params)
     if(verbose) print(sapply(params.new, round, 4))
     # Take differences (user input)
     if(diff.type == "abs"){
       diffs <- abs(params.new-params)
       b.diff <- max(abs(do.call(rbind, b.hat) - do.call(rbind, b)))
     }else if(diff.type == "abs.rel"){
       diffs <- abs(params.new-params)/(abs(params) + 1e-3)
       b.diff <- max(abs(do.call(rbind, b.hat) - do.call(rbind, b))/(abs(do.call(rbind, b)) + 1e-3))
     }
     diff <- max(diffs)
     # Message output (max relative diff)
     message("\nIteration ", iter + 1, " maximum difference: ", round(diff, 5))
     message("Largest change: ", names(params)[which(diffs==diff)])
     message("--> old: ", params[which(diffs==diff)], " new: ", params.new[which(diffs==diff)])
     message("Largest change in random effects: ", round(b.diff, 3))
     
     #' Update
     params <- params.new
     b <- b.hat
     D <- D.new; beta <- beta.new; gamma <- gamma.new; eta <- eta.new; # theta <- theta.new (NYI)
     gr <- rep(gamma, each = nK)
     l0 <- l0.new; l0u <- l0u.new; l0i <- l0i.new
     iter <- iter + 1
     if(collect.hist) iter.hist <- rbind(iter.hist, c(iter = iter, t(params)))
  }
  coeffs <- list(D = D, beta = beta, gamma = gamma, eta = eta, theta = theta, hazard = cbind(ft, l0))
  rtn <- list(REs = do.call(rbind, b.hat), coeffs = coeffs, 
              # Elapsed times //
              EMtime = round(sum(EMtime), 2),
              mvlme.time = mvlme.fit$elapsed.time,
              comp.time = round(proc.time()[3] - start.time, 2))
  if(collect.hist) rtn$history <- iter.hist
  if(post.process) message('\n\nSEs NYI!\n')
  rtn
}
