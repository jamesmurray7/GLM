rm(list=ls())

# Load data
load('/data/c0061461/disp-fits-3.RData')
load('/data/c0061461/disp-fits-8.RData')
load('/data/c0061461/disp-fits-minus5.RData')
 
# Function to plot median values against true
plot.medians.on.target <- function(x, target){
  # x a list
  meds <- do.call(c, lapply(x, function(xx) xx$median.estimate))
  plot(meds, pch = 20, ylab = bquote("Median estimate for " ~ delta==.(target)), xlab = '', xaxt = 'n')
  abline(h = target, lty = 5, col = 'red')
}

# Function to plot median with IQRs on target
plot.medIQR.on.target <- function(x, target){
  # x a list
  quan <- as.data.frame(do.call(rbind, lapply(lapply(x, el), quantile, probs = c(.25, .50, .75), na.rm = T)))
  ylims <- c(min(pmin(quan$`25%`, quan$`75%`)), max(pmax(quan$`25%`, quan$`75%`)))
  n <- 1:nrow(quan)
  plot(quan$`50%` ~ n, pch = 20, xaxt = 'n', ylim = ylims, ylab = bquote("Median [IQR] estimate for " ~ delta==.(target)), xlab = '')
  arrows(x0=n, y0 = quan$`25%`, y1 = quan$`75%`, code = 3, angle = 90, length = 0.02)
  abline(h = target, lty = 5, col = 'red')
}

# Function to plot median with user-define quantiles of values on target
plot.medquantile.on.target <- function(x, target, lower.quantile, upper.quantile){
  # x a list
  quan <- as.data.frame(do.call(rbind, lapply(lapply(x, el), quantile, probs = c(lower.quantile, .50, upper.quantile), na.rm = T)))
  ylims <- c(min(pmin(quan[, 1], quan[, 3])), max(pmax(quan[, 1], quan[, 3])))
  n <- 1:nrow(quan)
  
  # Work out coverage
  is.covered <- quan[, 1] <= target & quan[, 3] >= target
  cp <- paste0(round(100 * sum(is.covered)/nrow(quan), 2), '%')
  
  # Plot
  plot(quan$`50%` ~ n, pch = 20, xaxt = 'n', ylim = ylims, 
       xlab = paste0('Lower precentile: ', round(100 * lower.quantile, 2), 
                     ', upper precentile: ', round(100 * upper.quantile, 2),
                     '; coverage: ', cp),
       ylab = bquote("Median [IQR] estimate for " ~ delta==.(target)))
  arrows(x0=n, y0 = quan[, 1], y1 = quan[, 3], code = 3, angle = 90, length = 0.02)
  abline(h = target, lty = 5, col = 'red')
}

# Function to plot median with standard deviation on target
plot.medSD.on.target <- function(x, target, scale = qnorm(.975)){
  # x a list
  meds <- do.call(c, lapply(x, function(xx) xx$median.estimate))
  sds <- do.call(c, lapply(x, function(xx) xx$sd.estimates))
  df <- data.frame(med = meds, lb = meds - scale * sds, ub = meds + scale * sds)
  
  # Work out coverage
  is.covered <- df$lb <= target & df$ub >= target
  cp <- paste0(round(100 * sum(is.covered)/nrow(df), 2), '%')
  
  # Plot
  ylims <- c(min(df$lb), max(df$ub))
  n <- 1:nrow(df)
  plot(df$med ~ n, pch = 20, xaxt = 'n', ylim = ylims, 
       ylab = expression(paste("Median ", phantom(.)%+-%phantom(.), "sd estimate for ", delta)), 
       xlab = paste0('Scale: ', round(scale, 2), '; Coverage: ', cp))
  arrows(x0=n, y0 = df$lb, y1 = df$ub, code = 3, angle = 90, length = 0.02)
  abline(h = target, lty = 5, col = 'red')
}

# Function to overlay many densities
plot.many.densities <- function(x, target){
  # x a list
  ests <- lapply(x, el) # extract the i=1,...,valid(n) estimates
  # get densities and ranges for x/ylims
  dens <- lapply(ests, density, na.rm = T)
  ranges.x <- do.call(rbind, lapply(dens, function(x) range(x$x, na.rm = T)))
  ranges.y <- do.call(rbind, lapply(dens, function(x) range(x$y, na.rm = T)))
  xlims <- c(min(ranges.x[,1]), max(ranges.x[, 2]))
  ylims <- c(min(ranges.y[,1]), max(ranges.y[, 2]))
  
  # Plot first in list
  plot(dens[[1]], main = paste0('Densities of esimates for true dispersion: ', target),
       xlim = xlims, ylim = ylims, xlab = expression(delta))
  for(i in 2:length(dens)){
    lines(dens[[i]], col = i)
  }
  
}

# Function to plot elapsed times
plot.elapsed.time <- function(x){
  invisible(2+2)
}

# Tabulating coverages
tab.iqr.cp <- function(x, target, lower.quantile, upper.quantile){
  # x a list
  quan <- as.data.frame(do.call(rbind, lapply(lapply(x, el), quantile, probs = c(lower.quantile, .50, upper.quantile), na.rm = T)))
  ylims <- c(min(pmin(quan[, 1], quan[, 3])), max(pmax(quan[, 1], quan[, 3])))
  n <- 1:nrow(quan)
  
  # Work out coverage
  is.covered <- quan[, 1] <= target & quan[, 3] >= target
  cp <- round(100 * sum(is.covered)/nrow(quan), 2)
  
  cp
}

tab.sd.cp <- function(x, target, scale){
  # x a list
  # x a list
  meds <- do.call(c, lapply(x, function(xx) xx$median.estimate))
  sds <- do.call(c, lapply(x, function(xx) xx$sd.estimates))
  df <- data.frame(med = meds, lb = meds - scale * sds, ub = meds + scale * sds)
  
  # Work out coverage
  is.covered <- df$lb <= target & df$ub >= target
  cp <- round(100 * sum(is.covered)/nrow(df), 2)
  
  cp
}


