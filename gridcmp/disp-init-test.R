test.sets <- replicate(100,
                       simData_joint2(n = 250, delta = c(.8,-.2), 
                                      ntms = 10, theta = c(-2, .1), fup = 3,
                                      beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                      D = matrix(c(0.25, 0, 0, 0.00), 2, 2)),
                       simplify = F)

test.sets <- lapply(test.sets, el)
# At each time point and then lm'd ----------------------------------------
time.disps <- lapply(test.sets, function(x){
  time.disps <- with(x, tapply(Y, time, var))/with(x, tapply(Y, time, mean))
  df <- data.frame(time = as.numeric(names(time.disps)), disp = unname(time.disps))
  df
})

#' Okay, so correctly captures initial \delta value (0.8), and captures fact that dispersion
#' INCREASES over time. But this doesn't directly translate into true value -0.2 (value into \nu).

plot(sapply(time.disps, function(x) x$disp[1]))
abline(h=0.8, col = 'red', lty = 5)

# Linear models
time.disp.lm <- lapply(time.disps, function(x) lm(disp ~ time, data = x)$coeff)
points(sapply(time.disp.lm, function(x) x[1]), pch = 3) # Intercept values a little different.

collect.disps.lms <- do.call(rbind, time.disp.lm)

apply(collect.disps.lms, 2, quantile) # Log (intercept) for larger Y values?


# Turn this into a function -----------------------------------------------
disp.init <- function(formulas, disp.formula, data){
  
  data$Y <- data[,formulas$response]
  
  # Dispersion across time
  time.disps <- with(data, tapply(Y, time, var))/with(data, tapply(Y, time, mean))
  
  # Setup a dataframe
  df <- data.frame(time = as.numeric(names(time.disps)), disp = unname(time.disps))
  
  # Return linear model
  f <- paste0('disp', paste0(as.character(disp.formula), collapse = ''))
  lm(f, data = df)$coeff
  
}
formulas <- list(response='Y')

# For time specification
disp.formula1 <- ~ time

inits1 <- lapply(test.sets, function(x) disp.init(formulas, disp.formula1, x))
apply(do.call(rbind, inits1),2,mean) # Pretty good.


# Intercept only...
test.setsI <- replicate(100,
                       simData_joint2(n = 250, delta = c(.8, 0.0), 
                                      ntms = 10, theta = c(-2, .1), fup = 3,
                                      beta = c(0.0, 0.5, 0.05, -0.1), gamma = 0.6, zeta= c(0.0, -0.2),
                                      D = matrix(c(0.25, 0, 0, 0.00), 2, 2)),
                       simplify = F)
test.setsI <- lapply(test.setsI, el)
disp.formula2 <- ~ 1
inits2 <- lapply(test.setsI, function(x) disp.init(formulas, disp.formula2, x)) # Seems to overshoot a little



# Straight var(x)/mean(x) //
1/sapply(test.sets, function(x){
  y <- x$Y
  var(y)/mean(y)
}) ## Not good at all.

# var(x)/mean(x) across first Y terms only //
first.only <- sapply(test.sets, function(x){
  ys <- sapply(1:250, function(i){
    y <- x$data[x$data$id == i, 'Y']
    y[1]
  })
  var(ys)/mean(ys)
})

hist(first.only)
abline(v=.8,col='red') ## Quite bad, but target seems to be in one of the tails at least.

# var(x)/mean(x) as function of time, at start and end of each individual profile //
time.fun <- lapply(test.sets, function(x){
  ys <- lapply(1:250, function(i){
    dat <- x$data[x$data$id == i, ]
    y <- dat[, 'Y']; t <- dat[, 'time']; tl <- length(t)
    data.frame(t = c(t[1], t[tl]), Y = c(y[1], y[tl]))
  })
  do.call(rbind, ys)
})

# Carry out lm on each
estimates <- lapply(time.fun, function(x){
  time.disps <- with(x, tapply(Y, t, var))/with(x, tapply(Y, t, mean))
  df <- data.frame(disp = unname(time.disps), t = as.numeric(names(time.disps)))
  lm(disp ~ t,df)$coeff
})

estimates <- do.call(rbind, estimates)
apply(estimates,2,quantile) ## This quite good!

# What about just taking Y values from start and end of profile (at __tmax__)
time.fun2 <- lapply(test.sets, function(x){
  ys <- lapply(1:250, function(i){
    dat <- x$data[x$data$id == i, ]
    if(max(dat$time == 3)){
      y <- dat[,'Y']; t <- dat[,'time']
      return (data.frame(t = c(t[1], t[nrow(dat)]), Y = c(y[1], y[nrow(dat)])))
    }else{
      return(NULL)
    }
  })
  do.call(rbind, ys)
})

estimates2 <- lapply(time.fun2, function(x){
  time.disps <- with(x, tapply(Y, t, var))/with(x, tapply(Y, t, mean))
  df <- data.frame(disp = unname(time.disps), t = as.numeric(names(time.disps)))
  df
  lm(disp ~ t,df)$coeff
})

estimates2 <- do.call(rbind, estimates2)
apply(estimates2,2,quantile) ## One with first/last measurements better than 
apply(estimates,2,quantile)  ## just taking first/tmax measurements

# Complete profiles, finally -->
time.fun3 <- lapply(test.sets, function(x){
  yt <- x$data[,c('time', 'Y')]
  time.disps <- with(yt, tapply(Y, time, var))/with(yt, tapply(Y, time, mean))
  df <- data.frame(disp = unname(time.disps), t = as.numeric(names(time.disps)))
  df
  lm(disp ~ t,df)$coeff
})

estimates3 <- do.call(rbind, time.fun3)
apply(estimates2,2,quantile)  ## One with 'full' picture on dispersion better than 
apply(estimates,2,quantile)   ## just taking first/tmax measurements
apply(estimates3,2,quantile)  ## as well as just taking first/last measurements

par(mfrow=c(1,2))
apply(estimates3, 2, hist, breaks = 20) # Seems to over-estimate the time component a little;
par(mfrow=c(1,1))                       # All methods seemed to over-estimate/get this wrong; so maybe just set as zero.