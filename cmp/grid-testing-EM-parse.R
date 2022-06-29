load('/data/c0061461/cmp-fits-intonly-28-06-22-with-QUAD.RData')
sum(unlist(lapply(fits, is.null))) # no null

coeffs <- function(f){
  if(!is.null(f)){
    x <- f$coeffs
    D <- x$D[lower.tri(x$D, T)]
    b <- c(x$beta)
    d <- c(x$delta)
    g <- x$gamma
    z <- x$zeta
    setNames(c(D,b,d,g,z), names(f$SE))
  }
}
ests <- as.data.frame(do.call(rbind, lapply(fits, coeffs)))
# targets <- setNames(c(.2, 0, 0.05, 0.0, -0.1, 0.05, -0.1, 0.5, -0.1, .6, -.2),
#                     colnames(ests)) 

targets <- setNames(c(.16, 0.0, -0.1, 0.05, -0.1, 0.8, .6, -.2),
                    colnames(ests))

# Graphs ------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

ests.long <- ests %>% 
  pivot_longer(everything(), names_to = 'parameter', values_to = 'estimate')

ests.long <- left_join(ests.long, data.frame(parameter = names(targets), target = targets), 'parameter')

ests.long %>% 
  ggplot(aes(x = estimate)) + 
  geom_vline(aes(xintercept = target), lty = 5, colour = 'blue') + 
  geom_density() + 
  facet_wrap(~ parameter, scales = 'free') + 
  theme_light()

ests.long %>% 
  ggplot(aes(x = estimate)) + 
  geom_boxplot() + 
  geom_vline(aes(xintercept = target), lty = 5, colour = 'blue') + 
  facet_wrap(~ parameter, scales = 'free') + 
  theme_light()


# Tabulate ----------------------------------------------------------------
SE <- as.data.frame(do.call(rbind, lapply(fits, function(x) x$SE)))
lb <- ests - 1.96*SE;ub <- ests+1.96*SE
targets.mat <- apply(t(targets),2,rep,100)
CP <- colSums(lb <= targets.mat & ub >= targets.mat)/nrow(SE)
