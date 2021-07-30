#' ###
#' Bernhardt poisson
#' ###

# simulate some data
n <- 200
mi <- 5

# outcome
Y <- rpois(n*mi,8)

# X
id <- 1:200
time <- rep(0:(mi-1), n)
age <- rnorm(n)
X <- cbind(id = rep(id, each = mi), time, age = rep(age, each = mi))

# Data and fit 
data <- data.frame(Y, X)
glmerfit <- lme4::glmer(Y ~ time + age + (1 + time|id), data, family = 'poisson')
true.res <- lme4::ranef(glmerfit)$id

betaTrue <- c(2, 0.01, -0.002)

# Functions to extract subject-specific data
getXi <- function(data){
  uids <- unique(data$id)
  xx <- list()
  for(i in uids){
    xx[[i]] <- as.matrix(cbind(1, data[data$id == i, c("time", "age")]))
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
    xx[[i]] <- cbind(1, 0:(mi-1))
  }
  xx
}

# Set out complete data loglikelihood -------------------------------------
X <- getXi(data); Y <- getYi(data); Z <- getZi(data)
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

nlmfit <- nlm(ll2, b0[1,], Y[[1]], X[[1]], Z[[1]], beta, D)
Sigmai <- solve(-1 * d2b.ll(nlmfit$estimate, X[[1]], Z[[1]], betaTrue, D))

# EM ----------------------------------------------------------------------
gh <- statmod::gauss.quad.prob(9, dist = 'normal')
w <- gh$weights; v <- gh$nodes
# Construct a while loop
diff <- maxiter <- 100
tol <- 1e-4; iter <- 1
uids <- unique(data$id)
# Inits
b0 <- matrix(0.001, n, 2)
D <- matrix(lme4::VarCorr(glmerfit)$id, 2, 2)
beta <- c(summary(glmerfit)$coef[,1])
# param vector
vech <- function(x) x[lower.tri(x, diag = T)]
params <- c(c(beta), vech(D))
while(diff > tol && iter < maxiter){
  # Estep
  Sbetai <- matrix(NA, nr = length(uids), nc = 3)
  Ibetai <- D.newi <- list()
  bi.mat <- matrix(NA, nr = length(uids), nc = 2)
  for(i in uids){
    # bi <- ucminf(c(b0[i,]), ll, gradll, Y[[i]], X[[i]], Z[[i]], beta, D)$par
    bi <- nlm(ll2, b0[i,], Y[[i]], X[[i]], Z[[i]], beta, D)$estimate
    Sigmai <- solve(-1  * d2b.ll(bi, X[[i]], Z[[i]], beta, D))
    
    bi.mat[i,] <- bi
    # Update for D
    D.newi[[i]] <- Sigmai + tcrossprod(bi)
    
    # Update for beta
    tau <- sqrt(diag(tcrossprod(Z[[i]] %*% Sigmai, Z[[i]])))
    Sbetai.store <- matrix(NA, nr = 9, nc = 3)
    for(k in 1:9){
      Sbetai.store[k,] <- c(crossprod(-X[[i]], w[k] * exp(X[[i]] %*% beta + Z[[i]] %*% bi + v[k] * tau))) + Y[[i]] %*% X[[i]]
    }
    Sbetai[i, ] <- colSums(Sbetai.store)
    
    Ibetai.store <- list()
    for(k in 1:9){
      Ibetai.store[[k]] <- w[k] * crossprod(diag(as.numeric(exp(X[[i]] %*% beta + Z[[i]] %*% bi + v[k] * tau))) %*% X[[i]], X[[i]])
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


