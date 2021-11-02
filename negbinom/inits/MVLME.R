#' ######
#' Multivariate extension of linear mixed effects model
#' for negative bin (nK responses) 
#' ######


mvlme <- function(d, Y, X, Z,                       # data-specific args
                  inits.long, nK, q,                # parameters etc.
                  mvlme.tol = 5e-3, verbose = F){   # (relative) tolerance.
  
  uids <- unique(d$id)
  n <- length(uids)
  
  # Unpack inits.long
  D <- inits.long$D.init; Dinv <- solve(D)
  beta <- inits.long$beta.init
  theta <- inits.long$theta.init
  # b, the REs
  b0 <- as.matrix(Ranefs(inits.long))
  id_name <- grepl('id', colnames(b0))
  b0 <- lapply(1:nrow(b0), function(x) b0[x, !id_name])
  
  # And make parameter vector
  params <- c(vech(D), beta, theta)
  names(params) <- c(rep("D", length(vech(D))), names(beta), names(theta))
  
  # For ease of splitting b and beta
  b.inds <- split(seq(nK * 2), rep(1:nK, each = 2))
  beta.inds <- split(seq(length(beta)), rep(1:nK, each = 4))
  
  # We need to repeat theta Y.1 + Y.2 times for each id
  theta.list <- lapply(1:n, function(i){
    out <- c();
    for(k in 1:nK) out <- c(out, rep(theta[k], nrow(X[[i]])/2))
    out
  })
  
  # EM Algorithm
  diff <- 100; iter <- 0;
  mvlmetime <- c()
  EMstart <- proc.time()[3]
  while(diff > mvlme.tol){
    Estart <- proc.time()[3]
    
    b <- mapply(function(b0, Y, lfactY, X, Z, theta){
      ucminf::ucminf(b0, nb_ll, nb_grad,
                     Y, lfactY, X, Z, D, beta, theta, q, control = list())$par
    }, b0=b0, Y=Y, lfactY = lapply(Y, lfactorial), X = X, Z = Z, theta = theta.list,
    SIMPLIFY = F)
    
    Sigmai <- mapply(function(b, X, Z, Y, theta){
      solve(-1 * nb_hess(b, X, Z, Y, D, beta, theta))
    }, b = b, X = X, Z = Z, Y = Y, theta = theta.list,
    SIMPLIFY = F)
    
    # Variance V under `nbinom2`    << LIKELY NOT NEEDED >>
    # V <- mapply(function(X, Z, b, mi){
    #   vv <- numeric(nK)
    #   for(k in 1:nK){
    #     muk <- exp(X[[k]] %*% beta[beta.inds[[k]]] + Z[[k]] %*% b[b.inds[[k]]])
    #     vv[k] <- crossprod(muk, 1+muk/theta[k])
    #   }
    #   diag(x = rep(vv, mi), ncol = sum(mi))
    # }, X = Xk, Z = Zk, b = b, mi = mi, SIMPLIFY = F)
    # 
    # VV <- mapply(function(Z, V){
    #   solve(t(Z) %*% solve(V) %*% Z + solve(D))
    # }, V = V, Z = Z, SIMPLIFY = F)
    
    # E[bbT|]
    EbbT <- mapply(function(Sigmai, b){ # E[bbT|\lambda, D]
      Sigmai + tcrossprod(b)
    }, Sigmai = Sigmai, b = b, SIMPLIFY = F)
    
    # Steps to update \beta
    Sbeta <- mapply(function(b, Y, X, Z, theta){
      crossprod(X, Y) - crossprod(X, (Y + theta) * exp(X %*% beta + Z %*% b) / (exp(X %*% beta + Z %*% b) + theta))
    }, b = b, Y = Y, X = X, Z = Z, theta = theta.list, SIMPLIFY = F)
    
    Ibeta <- mapply(function(b, X, Z, Y, theta){ # NB this -I(\beta)
      mu <- exp(X %*% beta + Z %*% b) # This may well be slower than doing an 'outside' mapply for my first.
      term1 <- (Y + theta) * mu
      term2 <- mu + theta
      -1 * (crossprod(diag(x = c(term1 / term2)) %*% X, X) - crossprod(diag(x = c(term1 * mu / (term2 * term2))) %*% X, X))
    }, b = b, X = X, Z = Z, Y = Y, theta = theta.list, SIMPLIFY = F)
    
    #' M-step -------------------
    D.new <- Reduce('+', EbbT)/n
    
    beta.new <- beta + solve(-Reduce('+', Ibeta),
                             c(rowSums(do.call(cbind, Sbeta))))
    
    names(beta.new) <- names(beta)
    
    theta.new <- theta#+ solve(Reduce('+', Itheta), rowSums(Stheta))
    
    # New parameter vector
    params.new <- c(vech(D.new), beta.new, theta.new)
    names(params.new) <- names(params)
    # Take difference & report
    if(verbose) print(sapply(params.new, round, 4))
    diffs <- abs(params.new - params)/(abs(params) + 1e-3)
    diff <- max(diffs)
    message("\nIteration ", iter + 1, " largest relative difference = ", round(diff, 5))
    message("For: ", names(params)[which(diffs == diff)])
    message("---> Old: ", params[which(diffs == diff)], ", new: ", params.new[which(diffs == diff)])
    
    # Update parameters and loop
    params <- params.new
    b0 <- b
    D <- D.new; Dinv <- solve(D)
    theta <- theta.new
    beta <- beta.new; 
    iter <- iter + 1
  }
  message("Converged after ", iter, " iterations ")
  list(
    beta = beta, D = D, b = b, theta = theta,
    elapsed.time = proc.time()[3] - EMstart, EMtime = sum(mvlmetime)
  )
}