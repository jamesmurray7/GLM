mvlme <- function(d, Y, X, Z, inits.long, mvlme.tol = 5e-3){
  diff <- 100; iter <- 0;
  uids <- unique(d$id)
  n <- length(uids)
  
  # Unpack inits.long
  D <- inits.long$D.init; Dinv <- solve(D)
  beta <- inits.long$beta.init
  # b, the REs
  b <- as.matrix(Ranefs(inits.long))
  id_name <- grepl('id', colnames(b))
  b <- lapply(1:nrow(b), function(x) b[x, !id_name])
  
  # And make parameter vector
  params <- c(vech(D), beta)
  names(params) <- c(rep("D", length(vech(D))), names(beta))
  
  # EM Algorithm
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
    
    pois.vars <- mapply(function(b, X, Z){ # Assume dispersion = 1
      diag(c(exp(X %*% beta + Z %*% b)))
    }, b = b, X = X, Z = Z, SIMPLIFY = F)
    
    b <- mapply(function(Z, V, X, Y){
      ZD <- Z %*% D
      ZDZtV <- Z %*% D %*% t(Z) + V
      t(ZD) %*% solve(ZDZtV) %*% (Y - exp(X %*% beta))
    }, Z = Z, V = pois.vars, X = X, Y = Y, SIMPLIFY = F)
    
    # Cov(b)
    Sigmai <- mapply(function(Z, V){  # Total variation
      solve(crossprod(Z, solve(V) %*% Z) + Dinv)
    }, Z = Z, V = pois.vars, SIMPLIFY = F)  
    
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
    
    # b.new <- mapply(function(b, Sb, Ib){
    #   c(b + solve(-Ib, Sb))
    # }, b = b, Sb = Sb, Ib = Ib, SIMPLIFY = F)
    
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
    D <- D.new; Dinv <- solve(D)
    beta <- beta.new; 
    iter <- iter + 1
  }
  message("Converged after ", iter, " iterations ")
  list(
    beta = beta, D = D, b = b,
    V = pois.vars,
    elapsed.time = proc.time()[3] - EMstart, EMtime = sum(mvlmetime)
  )
}

