library(lme4)
df <- simData() # from ./Simulation.R

# Get X, Y, Z...
# Functions to extract subject-specific data
getXi <- function(data){
  uids <- unique(data$id)
  xx <- list()
  for(i in uids){
    xx[[i]] <- as.matrix(cbind(1, data[data$id == i, c("x1", "x2", 'x3')]))
  }
  xx
}

getYi <- function(data){
  uids <- unique(data$id)
  xx <- list()
  for(i in uids){
    xx[[i]] <- c(data[data$id == i, c("Y")])
  }
  xx
}

getZi <- function(data){
  uids <- unique(data$id)
  xx <- list()
  for(i in uids){
    xx[[i]] <- as.matrix(cbind(1, data[data$id == i, "x1"]))
  }
  xx
}

# Set out complete data loglikelihood -------------------------------------
X <- getXi(df); Y <- getYi(df); Z <- getZi(df)
# FOR SUBJECT i
ll <- function(b, Y, X, Z, beta, D){
  rtn <- -sum(exp(X %*% beta + Z %*% b)) + Y %*% (X %*% beta + Z %*% b) - sum(lfactorial(Y)) - 
    1/2 * log(2 * pi) - 1/2 * log(det(D)) - 1/2 * crossprod(b, solve(D) %*% b)
  -rtn
}

# FOR SUBJECT i, wrt b
gradll <- function(b, Y, X, Z, beta, D){
  grad <- -1 * c(crossprod(Z, exp(X %*% beta + Z %*% b))) + Y %*% Z - c(solve(D) %*% b)
  -grad
}

ll2 <- function(b, Y, X, Z, beta, D){
  rtn <- -sum(exp(X %*% beta + Z %*% b)) + Y %*% (X %*% beta + Z %*% b) - sum(lfactorial(Y)) - 
    1/2 * log(2 * pi) - 1/2 * log(det(D)) - 1/2 * crossprod(b, solve(D) %*% b)
  attr(rtn, 'gradient') <- -1 * (-1 * c(crossprod(Z, exp(X %*% beta + Z %*% b))) + Y %*% Z - c(solve(D) %*% b))
  -rtn
}

# Sigmai deriv step
d2b.ll <- function(b, X, Z, beta, D){
  crossprod(-1 * diag(as.numeric(exp(X %*% beta + Z %*% b))) %*% Z, Z) - solve(D)
}

#  EM ----------------------------------------------------------------------
# Quadrature points
gh.nodes <- 15
gh <- statmod::gauss.quad.prob(gh.nodes, dist = 'normal')
w <- gh$weights; v <- gh$nodes
# While/for loop-specific items
diff <- maxiter <- 100
tol <- 5e-3; iter <- 1
uids <- unique(df$id); n <- length(uids)
# Initial conditions - bad
b0 <- matrix(0.001, n, 2)
D <- diag(2)
beta <- c(0.1, 0.05, -0.1, 0.1)
# initial conditions - better
glmfit <- glmer(Y ~ x1 + x2 + x3 + (1 + x1|id), data = df, family = poisson)
beta <- fixef(glmfit)
D <- matrix(as.numeric(VarCorr(glmfit)$id), 2, 2)
b0 <- as.matrix(ranef(glmfit)$id)
# param vector
vech <- function(x) x[lower.tri(x, diag = T)]
params <- c(c(beta), vech(D))
while(diff > tol && iter <= maxiter){
  # Estep
  Sbetai <- matrix(NA, nr = length(uids), nc = 4)
  Ibetai <- D.newi <- list()
  bi.mat <- matrix(NA, nr = length(uids), nc = 2)
  for(i in uids){
    # bi <- ucminf(c(b0[i,]), ll, gradll, Y[[i]], X[[i]], Z[[i]], beta, D)$par
    # bi <- nlm(ll2, b0[i,], Y[[i]], X[[i]], Z[[i]], beta, D)$estimate
    bi <- ucminf::ucminf(b0[i,], ll, gradll, 
                         Y[[i]], X[[i]], Z[[i]], beta, D)$par
    Sigmai <- solve(-1  * d2b.ll(bi, X[[i]], Z[[i]], beta, D))
    
    bi.mat[i,] <- bi
    # Update for D
    D.newi[[i]] <- Sigmai + tcrossprod(bi)
    
    # Update for beta
    tau <- sqrt(diag(tcrossprod(Z[[i]] %*% Sigmai, Z[[i]])))
    Sbetai.store <- matrix(NA, nr = gh.nodes, nc = 5)
    for(k in 1:gh.nodes){
      Sbetai.store[k,] <- w[k] * exp(X[[i]] %*% beta + Z[[i]] %*% bi + v[k] * tau)
    }
    Sbetai[i, ] <- crossprod(X[[i]], Y[[i]]) - crossprod(X[[i]], colSums(Sbetai.store))
    
    Ibetai.store <- list()
    for(k in 1:gh.nodes){
      Ibetai.store[[k]] <-crossprod(diag(as.numeric(w[k] * exp(X[[i]] %*% beta + Z[[i]] %*% bi + v[k] * tau))) %*% X[[i]], X[[i]])
    }
    Ibetai[[i]] <- Reduce('+', Ibetai.store)
  }
  
  D.new <- Reduce('+', D.newi)/n
  
  # beta update
  Sbeta <- colSums(Sbetai)
  Ibeta <- Reduce('+', Ibetai)
  beta.new <- beta + solve(Ibeta, Sbeta)
  
  # New params
  params.new <- c(c(beta.new), vech(D.new))
  diff <- max(abs(params.new - params) / (abs(params) + 1e-3))
  message("\nIteration ", iter, " Relative difference ", round(diff, 5),
          "\nNew beta: ", paste0(sapply(beta.new, round, 4), collapse =' '),
          "\nNew vech(D): ", paste0(sapply(vech(D.new), round, 4), collapse = ' '))
  beta <- beta.new; D <- D.new
  b0 <- bi.mat
  params <- params.new
  iter <- iter + 1
}
