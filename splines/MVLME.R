#' ###
#' MVLME.R // Multivariate linear mixed effects
#' ---
#' Uses longit initial conditions from inits.R
#' and fits MVLME to find optimal starting values for joint model with quadratic RE structure on all responses
#' ###

# Prerequisites -----------------------------------------------------------
# sourceCpp("./mvlmecpp.cpp") # Loading Rcpp Source files...

if(!"vech"%in%ls()) vech <- function(x) x[lower.tri(x, diag = T)]

splitbbT <- function(bbT, nK, degree){
  uids <- length(bbT)
  out <- list()
  for(i in 1:uids){
    out[[i]] <- lapply(split(seq(nK * (degree + 1)), rep(1:nK, each = (degree + 1))), function(x) bbT[[i]][x, x])
  }
  out
}

mvlme <- function(d, Y, X, Z, # BLOCK MATRICES
                  Yk, Xk, Zk, m, # Single matrices
                  inits.long, nK, degree, 
                  tol.mvlme = 5e-3){
  # Data
  diff <- 100; iter <- 0;
  uids <- unique(d$id)
  n <- length(uids)
  XtX <- lapply(X, crossprod)
  
  # Unpack inits.long
  D <- inits.long$D.init
  var.e <- inits.long$var.e.init
  beta <- inits.long$beta.init; beta.mat <- matrix(beta, nr = nK, byrow = T)
  V <- lapply(m, function(i) {
    diag(x = rep(var.e, i), ncol = sum(i))
  })
  
  # And make parameter vector
  params <- c(vech(D), var.e, beta)
  
  # EM Algorithm
  mvlmetime <- c()
  EMstart <- proc.time()[3]
  
  while(diff > tol.mvlme){
    Estart <- proc.time()[3]
    #' E-step -------------------
    # E[b]
    b <- Eb(Y, X, Z, V, D, beta, n)
    bmat <- lapply(b, function(x) matrix(x, nc = (degree + 1), byrow = T))
    
    # Sigmai
    Sigmai <- covb(Z, V, solve(D), n)
    
    # E[bbT]
    bbT <- EbbT(b, Sigmai, n)
    bbTk <- splitbbT(bbT, nK, degree)
    
    #' M-step -------------------
    
    # D //
    D.new <- Reduce('+', bbT) / n
    
    # \beta //
    beta.rhs <- betaRHS(X, Y, Z, b, n)
    beta.new <- solve(Reduce('+', XtX)) %*% Reduce('+', beta.rhs)
    
    # var.e
    var.e.new <- Ee(Yk, Xk, Zk, beta.mat, bmat, bbTk, n, nK)/colSums(do.call(rbind, m))
    
    mvlmetime[iter + 1] <- proc.time()[3] - Estart
    
    # New parameter vector
    params.new <- c(vech(D.new), var.e.new, beta.new)
    # Take difference & report
    diffs <- abs(params.new - params)/(abs(params) + 1e-3)
    diff <- max(diffs)
    message("\nIteration ", iter + 1, " largest relative difference = ", round(diff, 5))
    
    # Update parameters and loop
    params <- params.new
    D <- D.new; var.e <- var.e.new; 
    beta <- beta.new; beta.mat <- matrix(c(beta), nr = nK, byrow = T)
    V <- lapply(m, function(i) {
      diag(x = rep(var.e.new, i), ncol = sum(i))
    })
    iter <- iter + 1
  }
  message("Converged after ", iter, " iterations ")
  list(
    beta = beta, D = D, var.e = var.e, b = b,
    V = V, XtX = XtX,
    elapsed.time = proc.time()[3] - EMstart, EMtime = sum(mvlmetime)
  )
}