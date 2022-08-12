#' 1)

data <- test$data
n <- 250

# Get data matrices
X <- Y <- Z <- G <- vector('list', n)
for(i in 1:n){
  i.dat <- subset(data, id == i)
  Y[[i]] <- i.dat$Y
  X[[i]] <- model.matrix(~time + cont + bin, i.dat)
  Z[[i]] <- G[[i]] <- model.matrix(~ 1, i.dat)
}
 

fit <- glmmTMB(Y ~ time + cont + bin + (1|id),
               data, family = 'poisson')
beta <- fixef(fit)$cond
b <- ranef(fit)$cond$id
b <- lapply(1:n, function(i) b[i,])

eta <- mapply(function(X, Z, b) X %*% beta + Z %*% b, X = X, Z = Z, b = b)
mu <- lapply(eta, exp)

vars <- lapply(Y, var)
use.var <- sapply(vars, is.na)

# 'double optim' ----------------------------------------------------------
# We want to minimise error between `actual` Var[Y] and its expression.
summax <- max(sapply(Y, max)) * 2 # can change this around ... 
f <- function(delta, beta, b, Y, X, Z, G){
  nu <- exp(G %*% delta)
  mu <- exp(X %*% beta + Z %*% b)
  lambda <- lambda_appx(mu, nu, summax)
  EY <- lambda^(1/nu) - (nu - 1)/(2 * nu)
  var(Y) - EY/nu
}

# Function to find either side of interval to look for root
# otherwise uniroot will just error, and no conclusions can be drawn.
determine.interval <- function(mu){
  lower <- -2; upper <- 2 # Nominally choose lower and upper values of interval to check over.
  test.lower <- lambda_appx(mu, exp(lower), summax)
  test.upper <- lambda_appx(mu, exp(upper), summax)
  
  flag.lower <- is.nan(test.lower); flag.upper <- is.nan(test.upper)
  while(flag.lower){
    lower <- lower + 0.2
    flag.lower <- is.nan(lambda_appx(mu, exp(lower), summax))
  }
  while(flag.upper){
    upper <- upper - 0.2
    flag.upper <- is.nan(lambda_appx(mu, exp(upper), summax))
  }
  c(lower = lower, upper = upper)
}

# UNIROOTS ----
# Only use Ys with non-zero variance and profile length > 1
Y.to.use <- sapply(1:n, function(i){
  y <- Y[[i]]
  length(y) > 1 && var(y) != 0
})
uniroot.deltas <- vector('list', sum(Y.to.use))
inds <- which(Y.to.use)
p <- 1
for(i in inds){
  y <- Y[[i]]; l <- length(y)
  x <- X[[i]]; z <- Z[[i]]; g <- G[[i]]; bb <- b[[i]]
  mu.i <- exp(x %*% beta + z %*% bb)
  out <- numeric(l)
  if(var(y)/mean(y) != 0 && !is.nan(var(y)/mean(y))){ # if dispersion detected solve for nu, else return zeros
    for(j in 1:l){
      xj <- x[j,,drop=F]; zj <- z[j,,drop=F]; gj <- g[j,,drop=F]
      .interval <- determine.interval(mu.i[j])
      out[j] <- tryCatch(uniroot(f, interval = c(.interval[1], .interval[2]), 
                        beta = beta, b = bb, Y = y, X = xj, Z = zj, G = gj)$root,
                        error = function(e) NA)
    }
  }
  uniroot.deltas[[p]] <- out
  cat(sprintf('%d/%d\r', p, sum(Y.to.use)))
  p <- p + 1
}

# Median/mean values
uniroot.deltas


# Different parameterisation of E/Var[Y] ----------------------------------

Ey <- function(mu, nu, summax){
  lam <- lambda_appx(mu, nu, summax)
  logZ <- logZ_c_scalar(log(lam), nu, summax)
  rhs <- numeric(summax)
  for(j in 1:summax){
    rhs[j] <- exp(log(j-1) + (j-1) * log(lam) - nu * lgamma(j) - logZ)
  }
  sum(rhs)
}

VarY <- function(mu, nu, summax, Ey){
  lam <- lambda_appx(mu, nu, summax)
  logZ <- logZ_c_scalar(log(lam), nu, summax)
  rhs <- numeric(summax)
  for(j in 1:summax){
    rhs[j] <- exp(
      2 * log(j - 1) + (j - 1) * log(lam) - nu * lgamma(j) - logZ
    )
  }
  sum(rhs) - Ey^2
}

#  Seek to minimise difference between real and expected variance.
ff <- function(delta, Y, G, mu, summax){
  nu <- exp(G %*% delta)
  YY <- Ey(mu, nu, summax)
  vY <- VarY(mu, nu, summax, YY)
  var(Y) - vY
}


uniroot.deltas2 <- vector('list', length(inds))
p <- 1
for(i in inds){
  y <- Y[[i]]; l <- length(y)
  x <- X[[i]]; z <- Z[[i]]; g <- G[[i]]; bb <- b[[i]]
  mu.i <- exp(x %*% beta + z %*% bb)
  out <- numeric(l)
  if(var(y)/mean(y) != 0 && !is.nan(var(y)/mean(y))){ # if dispersion detected solve for nu, else return zeros
    for(j in 1:l){
      xj <- x[j,,drop=F]; zj <- z[j,,drop=F]; gj <- g[j,,drop=F]
      # .interval <- determine.interval(mu.i[j])
      out[j] <- tryCatch(uniroot(ff, interval = c(-1, 3), 
                                 Y = y, G = gj, mu = mu.i[j], summax = 50)$root,
                         error = function(e) NA)
    }
  }
  uniroot.deltas2[[p]] <- out
  cat(sprintf('%d/%d\r', p, length(inds)))
  p <- p + 1
}

# Compare these two methods ----
# medians
med1 <- sapply(uniroot.deltas, median, na.rm = T)
med2 <- sapply(uniroot.deltas2, median, na.rm = T)
mean(med2, na.rm = T)
mean(med1, na.rm = T)

# Means
mean1 <- sapply(uniroot.deltas, mean, na.rm = T)
mean2 <- sapply(uniroot.deltas2, mean, na.rm = T)
plot(mean1, mean2)
abline(0, 1)
median(mean1)
median(mean2, na.rm = T) # From this, average median delta root...
