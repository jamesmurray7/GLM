#' ######
#' Multivariate extension of linear mixed effects model
#' for negative bin (nK responses) 
#' ######

# Define functions we need for update to theta
SthetaFun <- function(theta, Y, mu){ 
  digamma(theta+Y) - digamma(theta) + 1 + log(theta) - (log(mu + theta) + (theta + Y)/(mu + theta))
}
IthetaFun <- function(theta, Y, mu){
  ones <- rep(1, length(theta)); l1 <- length(ones)
  term1 <- diag(x = c(trigamma(theta + Y)), nc = l1) - diag(x = c(trigamma(theta)), nc = l1) + diag(x = c(ones / theta), nc = l1)
  term2 <- diag(x = c(ones / (mu + theta)), nc = l1)
  term3 <- diag(x = c((theta + Y) / ((mu + theta) * (mu + theta))), nc = l1)
  term1 - (2 * term2 - term3)
}



mvlme <- function(d, Y, X, Z,                       # data-specific args
                  Yk, Xk, Zk, mi,                   # Split-out data matrices for the K responses
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
  
  # We need to repeat theta_k m_{ik} times for each id...
  theta.list <- lapply(1:n, function(i){
    out <- list();
    for(k in 1:nK) out[[k]] <- rep(theta[k], nrow(Xk[[i]][[k]]))
    out
  })
  
  LL <- apply(do.call(rbind, mi), 2, max)  # Max profile length.
  
  # EM Algorithm
  diff <- 100; iter <- 0;
  mvlmetime <- c()
  EMstart <- proc.time()[3]
  while(diff > mvlme.tol){
    Estart <- proc.time()[3]
    
    b <- mapply(function(b0, Y, lfactY, X, Z, theta){
      ucminf::ucminf(b0, nb_ll, nb_grad,
                     Y, lfactY, X, Z, D, beta, theta, q)$par
    }, b0=b0, Y=Y, lfactY = lapply(Y, lfactorial), X = X, Z = Z, theta = lapply(theta.list, function(x) do.call(c, x)),
    SIMPLIFY = F)
    
    Sigmai <- mapply(function(b, X, Z, Y, theta){
      solve(-1 * nb_hess(b, X, Z, Y, D, beta, theta))
    }, b = b, X = X, Z = Z, Y = Y, theta = lapply(theta.list, function(x) do.call(c, x)),
    SIMPLIFY = F)
    
    # E[bbT|]
    EbbT <- mapply(function(Sigmai, b){ # E[bbT|\lambda, D]
      Sigmai + tcrossprod(b)
    }, Sigmai = Sigmai, b = b, SIMPLIFY = F)
    
    # Steps to update \beta
    Sbeta <- mapply(function(b, Y, X, Z, theta){
      crossprod(X, Y) - crossprod(X, (Y + theta) * exp(X %*% beta + Z %*% b) / (exp(X %*% beta + Z %*% b) + theta))
    }, b = b, Y = Y, X = X, Z = Z, theta = lapply(theta.list, function(x) do.call(c, x)), SIMPLIFY = F)
    
    Ibeta <- mapply(function(b, X, Z, Y, theta){ # NB this -I(\beta)
      mu <- exp(X %*% beta + Z %*% b) # This may well be slower than doing an 'outside' mapply for my first.
      term1 <- (Y + theta) * mu
      term2 <- mu + theta
      -1 * (crossprod(diag(x = c(term1 / term2)) %*% X, X) - crossprod(diag(x = c(term1 * mu / (term2 * term2))) %*% X, X))
    }, b = b, X = X, Z = Z, Y = Y, theta = lapply(theta.list, function(x) do.call(c, x)), SIMPLIFY = F)
    
    # Steps to update \theta
    muk <- mapply(function(X, Z, b){
      out <- list()
      for(k in 1:nK){
        out[[k]] <- exp(X[[k]] %*% beta[beta.inds[[k]]] + Z[[k]] %*% b[b.inds[[k]]])
      }
      out
    }, X = Xk, Z = Zk, b = b, SIMPLIFY = F)
    
    Sthetas <- Ithetas <- list()
    for(i in 1:n){
      thisS <- thisI <- list()
      for(k in 1:nK){
        tempSk <- SthetaFun(theta.list[[i]][[k]], Yk[[i]][, k], muk[[i]][[k]])
        tempSk <- c(tempSk, rep(0, LL[k] - length(tempSk)))           # Concatenate with a load of zeroes so we can rbind later
        thisS[[k]] <- tempSk
        tempIk <- -1 * IthetaFund(theta.list[[i]][[k]], Yk[[i]][, k], muk[[i]][[k]])
        tempIk <- diag(tempIk)
        tempIk <- diag(x = c(tempIk, rep(0, LL[k] - length(tempIk)))) # Concatenate with a load of zeroes so we can rbind later
        thisI[[k]] <- tempIk
      }
      Sthetas[[i]] <- thisS
      Ithetas[[i]] <- thisI
    }
    
    #' M-step -------------------
    D.new <- Reduce('+', EbbT)/n
    
    beta.new <- beta + solve(-Reduce('+', Ibeta),
                             c(rowSums(do.call(cbind, Sbeta))))
    
    Stheta <- lapply(1:nK, function(k) colSums(do.call(rbind, lapply(Sthetas, '[[', k))))
    Itheta <- lapply(1:nK, function(k) Reduce('+', lapply(Ithetas, '[[', k)))
    
    theta.new <- lapply(1:nK, function(k){
      rep(theta[k], LL[k]) + solve(Itheta[[k]], Stheta[[k]])
    })
    theta.new <- sapply(1:nK, function(k) mean(theta.new[[k]]))
    
    names(beta.new) <- names(beta)
    names(theta.new) <- names(theta)
    
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
    theta.list <- lapply(1:n, function(i){
      out <- list();
      for(k in 1:nK) out[[k]] <- rep(theta[k], nrow(Xk[[i]][[k]]))
      out
    })
    beta <- beta.new; 
    iter <- iter + 1
  }
  message("Converged after ", iter, " iterations ")
  
  list(
    beta = beta, D = D, b = b, theta = theta,
    elapsed.time = proc.time()[3] - EMstart, EMtime = sum(mvlmetime)
  )
}
