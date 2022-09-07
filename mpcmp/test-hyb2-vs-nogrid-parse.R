# Load --------------------------------------------------------------------
rm(list=ls())
load('/data/c0061461/hybrid2_fits.RData')
load('/data/c0061461/no-grids_fits.RData')

# Check for null fits -----------------------------------------------------
to.remove <- unique(c(which(unlist(lapply(fits1, is.null))), which(unlist(lapply(fits2, is.null)))))
# Can't compare if one fit doesn't exist!

# Function to extract estimates of interest -------------------------------
targets <- c(.25, 2, -0.1, 0.1, -0.2, 0.8, 0.6, -0.2)
extract.model.estimates <- function(fit){
  if(is.null(fit)) return(NA)
  
  #' Coefficients & SEs.
  SE <- fit$SE; nSE <- names(SE)
  Omega <- setNames(c(c(fit$coeffs$D), c(fit$coeffs$beta), c(fit$coeffs$delta), 
                      c(fit$coeffs$gamma), c(fit$coeffs$zeta)),
                    nSE)
  
  #' Working out coverage.
  lb <- Omega - qnorm(.975) * SE; ub <- Omega + qnorm(.975) * SE
  CP <- lb <= targets & ub >= targets
  
  return(list(
    coeffs = Omega,
    CP = CP
  ))
}

f1 <- lapply(fits1, extract.model.estimates)
f2 <- lapply(fits2, extract.model.estimates)

# Compare estimates -------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
all.ests <- lapply(setdiff(1:50, to.remove), function(i){
  a <- data.frame(estimate = f1[[i]]$coeffs, target = targets, method = 'Hybrid grid', i = i)
  a <- cbind(parameter = row.names(a), a)
  row.names(a) <- NULL
  b <- data.frame(estimate = f2[[i]]$coeffs, target = targets, method = 'No grid', i = i)
  b <- cbind(parameter = row.names(b), b)
  row.names(b) <- NULL
  
  rbind(a, b)
})


all.ests <- do.call(rbind, all.ests)

ggplot(all.ests, aes(x = estimate, colour = method)) + 
  geom_vline(aes(xintercept = target), col = 'black', lty = 5) + 
  geom_density(alpha = .5) + 
  facet_wrap(~ parameter, scales = 'free') + 
  theme_light()

# These look quite similar then (at least at face-value), 

# Coverages...
all.CPs <- lapply(setdiff(1:50, to.remove), function(i){
  a <- data.frame(estimate = as.numeric(f1[[i]]$CP), target = targets, method = 'Hybrid grid', i = i)
  a$parameter <- names(f1[[i]]$coeffs)
  b <- data.frame(estimate = as.numeric(f2[[i]]$CP), target = targets, method = 'No grid', i = i)
  b$parameter <- names(f1[[i]]$coeffs)
  rbind(a, b)
})

all.CPs <- do.call(rbind, all.CPs)

all.CPs %>% 
  group_by(method, parameter) %>% 
  summarise(.groups = 'keep',
            CP = sum(estimate)/length(setdiff(1:50, to.remove))) # Exactly the same!


# Times taken... ----------------------------------------------------------
# Start with total comp. times -->
comp.times <- lapply(setdiff(1:50, to.remove), function(i){
  f1 <- fits1[[i]]$elapsed.time; f2 <- fits2[[i]]$elapsed.time
  t1 <- unname(f1[grepl('Total', names(f1))]); t2 <- unname(f2[grepl('Total', names(f2))])
  c(`No grid` = t2, `Hybrid grid` = t1)
})
comp.times <- do.call(rbind, comp.times)
plot(comp.times[,1], comp.times[,2], pch = 20, ylab = 'Hybrid', xlab = 'No grids',
     main = "Total computation time for underdispersed data")
abline(0, 1, col = 'red', ) # So hybrid nearly always slower.

# What about EM times?
EM.times <- lapply(setdiff(1:50, to.remove), function(i){
  f1 <- fits1[[i]]$elapsed.time; f2 <- fits2[[i]]$elapsed.time
  t1 <- unname(f1[grepl('EM ', names(f1))]); t2 <- unname(f2[grepl('EM ', names(f2))])
  c(`No grid` = t2, `Hybrid grid` = t1)
})
EM.times <- do.call(rbind, EM.times)
plot(EM.times[,1], EM.times[,2], pch = 20, ylab = 'Hybrid', xlab = 'No grids',
     main = "EM time for light underdispersed data")
abline(0, 1, col = 'red', ) # Again, hybrid nearly always slower.

# Is this simply because too few iterations occur?
iters <- lapply(setdiff(1:50, to.remove), function(i){
  f1 <- fits1[[i]]$iter; f2 <- fits2[[i]]$iter
  c(`No grid` = f2, `Hybrid` = f1)
})
iters <- do.call(rbind, iters)
iters[,1] > iters[,2]
EM.times[iters[,1] > iters[,2],]
EM.times[EM.times[,2] < EM.times[,1],]

same.iters <- iters[,1] == iters[,2]
longer.hybr <- iters[,2] > iters[,1]
longer.nogr <- iters[,1] > iters[,2]

EM.times[same.iters,]
EM.same <- as.data.frame(EM.times[same.iters,])
EM.lh <- as.data.frame(EM.times[longer.hybr,])
EM.ng <- as.data.frame(EM.times[longer.nogr,])
plot(EM.same); abline(0,1)
plot(EM.lh);   abline(0,1)
plot(EM.ng);   abline(0,1)

# Is there some `crossover` for non-identical iters?
EM.diffs <- as.data.frame(EM.times[!same.iters,])
EM.diffs <- cbind(diff.iters = iters[!same.iters, 'Hybrid'] - iters[!same.iters, 'No grid'], EM.diffs)
plot(EM.diffs$diff.iters, EM.diffs$`Hybrid grid` - EM.diffs$`No grid`, pch= 20,
     xlab = 'Hybrid - Nogrid iteration',
     ylab = 'Hybrid - Nogrid elapsed time (s)')

EM.diffs <- as.data.frame(EM.times)
EM.diffs <- cbind(diff.iters = iters[, 'Hybrid'] - iters[, 'No grid'], EM.diffs)
plot(EM.diffs$diff.iters, EM.diffs$`Hybrid grid` - EM.diffs$`No grid`, pch= 20,
     xlab = 'Hybrid - Nogrid iteration',
     ylab = 'Hybrid - Nogrid elapsed time (s)', main = 'EM times')
abline(h = 0, lty = 5, col = 'lightblue')

# Contrast w/ total computation times...
comp.times <- as.data.frame(comp.times)
comp.times <- cbind(diff.iters = iters[, 'Hybrid'] - iters[, 'No grid'], comp.times)
plot(comp.times$diff.iters, comp.times$`Hybrid grid` - comp.times$`No grid`, pch= 20,
     xlab = 'Hybrid - Nogrid iteration',
     ylab = 'Hybrid - Nogrid elapsed time (s)', main = 'Total computation time')
abline(h = 0, lty = 5, col = 'lightblue')
