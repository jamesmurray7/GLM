mvlme <- function(d, Y, X, Z, inits.long, nK, q,    # data-specific args
                  mvlme.tol = 5e-3){                # (relative) tolerance.
  
  uids <- unique(d$id)
  n <- length(uids)
  
  # Unpack inits.long
  D <- inits.long$D.init; Dinv <- solve(D)
  beta <- inits.long$beta.init
  # b, the REs
  b0 <- as.matrix(Ranefs(inits.long))
  id_name <- grepl('id', colnames(b0))
  b0 <- lapply(1:nrow(b0), function(x) b0[x, !id_name])
  
  # And make parameter vector
  params <- c(vech(D), beta)
  names(params) <- c(rep("D", length(vech(D))), names(beta))
  
  gh <- statmod::gauss.quad.prob(9, 'normal')
  w <- gh$w; v <- gh$n
  
  # EM Algorithm
  diff <- 100; iter <- 0;
  mvlmetime <- c()
  EMstart <- proc.time()[3]
  while(diff > mvlme.tol){
    Estart <- proc.time()[3]
    #' E-step -------------------
    # E[b]
    # Sb <- mapply(function(b, Y, X, Z){
    #   crossprod(Z, Y) - crossprod(Z, exp(X %*% beta - Z %*% b)) - Dinv %*% b
    # }, b = b, Y = Y, X = X, Z = Z, SIMPLIFY = F)
    # 
    # Ib <- mapply(function(b, X, Z){
    #   -crossprod(diag(c(exp(X %*% beta + Z %*% b))) %*% Z,
    #              Z) - Dinv
    # }, b = b, X = X, Z = Z, SIMPLIFY = F)
    
    b <- mapply(function(b0, Y, lfactY, X, Z){
      ucminf::ucminf(b0, po_ll, po_grad,
                     Y, lfactY, X, Z, D, beta, q, control = list())$par
    }, b0=b0, Y=Y, lfactY = lapply(Y, lfactorial), X = X, Z = Z,
    SIMPLIFY = F)
    
    Sigmai <- mapply(function(b, X, Z){
      solve(-1 * po_hess(b, X, Z, D, beta))
    }, b = b, X = X, Z = Z,
    SIMPLIFY = F)
    
    # Old -v- (?)
    # pois.vars <- mapply(function(b, X, Z){ # Assume dispersion = 1
    #   diag(c(exp(X %*% beta + Z %*% b)))
    # }, b = bb, X = X, Z = Z, SIMPLIFY = F)
    
    # Cov(b)
    # Sigmai <- mapply(function(Z, V){  # Total variation
    #   solve(crossprod(Z, solve(V) %*% Z) + Dinv)
    # }, Z = Z, V = pois.vars, SIMPLIFY = F)

    # E[bbT|\lambda, D]
    EbbT <- mapply(function(Sigmai, b){ # E[bbT|\lambda, D]
      Sigmai + tcrossprod(b)
    }, Sigmai = Sigmai, b = b, SIMPLIFY = F)
    
    # Steps to update \beta
    Sbeta <- mapply(function(b, Y, X, Z){
      crossprod(X, Y) - crossprod(X, exp(X %*% beta + Z %*% b))
    }, b = b, Y = Y, X = X, Z = Z, SIMPLIFY = F)
    
    Ibeta <- mapply(function(b, X, Z){
      -crossprod(diag(c(exp(X %*% beta + Z %*% b))) %*% X,
                 X)
    }, b = b, X = X, Z = Z, SIMPLIFY = F)
    
    #' M-step -------------------
    D.new <- Reduce('+', EbbT)/n
    
    beta.new <- beta + solve(-Reduce('+', Ibeta),
                             c(rowSums(do.call(cbind, Sbeta))))
    
    names(beta.new) <- names(beta)
    
    # New parameter vector
    params.new <- c(vech(D.new), beta.new)
    names(params.new) <- names(params)
    # Take difference & report
    print(sapply(params.new, round, 4))
    diffs <- abs(params.new - params)/(abs(params) + 1e-3)
    diff <- max(diffs)
    message("\nIteration ", iter + 1, " largest relative difference = ", round(diff, 5))
    message("For: ", names(params)[which(diffs == diff)])
    message("---> Old: ", params[which(diffs == diff)], ", new: ", params.new[which(diffs == diff)])
    
    # Update parameters and loop
    params <- params.new
    b0 <- b
    D <- D.new; Dinv <- solve(D)
    beta <- beta.new; 
    iter <- iter + 1
  }
  message("Converged after ", iter, " iterations ")
  list(
    beta = beta, D = D, b = b,
    elapsed.time = proc.time()[3] - EMstart, EMtime = sum(mvlmetime)
  )
}
