#' #####
#' longs.R
#' ---
#' Determining best longitudinal fits for each biomarker in PBC data
#' #####

rm(list=ls())
library(glmmTMB)
library(tidyverse)
library(ggplot2)
load('../PBC-case-study/PBC.RData')

markers <- data.frame(
  marker = c('serBilir', 'SGOT', 'albumin', 'prothrombin', 'platelets', 'alkaline', 'spiders', 'ascites', 'hepatomegaly'),
  family = c(rep('gaussian', 4), rep('poisson', 2), rep('binomial', 3)),
  stringsAsFactors = F
)

pbc$serBilir <- log(pbc$serBilir)
pbc$SGOT <- log(pbc$SGOT)
pbc$prothrombin <- (.1*pbc$prothrombin)^(-4)

# Plotting trajectories ---------------------------------------------------
plotter <- function(val, group){
  ggplot(pbc[pbc$status==1,], aes(x = time, y = {{val}}, group = id)) + 
    geom_line(colour = 'grey', alpha = .25) + 
    geom_smooth(aes(group={{group}}), method = 'loess', colour = 'red', formula=y~x) +
    # geom_smooth(aes(group = {{group}}), colour = 'red', method = 'lm',
    #             formula = formula) + 
    theme_light()
}

plotter(serBilir, NULL)    # curvilinear
plotter(albumin, NULL)     # approximately linear
plotter(SGOT, NULL)        # approximately linear
plotter(prothrombin, NULL) # quadratic?
plotter(platelets, drug)   # curvilinear
plotter(alkaline, sex)     # curvilinear

# Functions to obtain glmmTMB fit -----------------------------------------
makeFormula <- function(y, type = 'linear'){
  rhs <- switch(type,
                linear = ' ~ drug * time + (1 + time|id)',
                intercept = '~ drug * time + (1|id)',
                quadratic = '~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id)',
                spline = '~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id)')
  as.formula(paste0(y, rhs))
}

TMBfit <- function(y, type){
  form <- makeFormula(y, type)
  family <- markers[markers$marker == y, 'family']
  fit <- glmmTMB(form, data = pbc, family = family)
  fit
}


out <- setNames(vector('list', length = length(markers$marker)),
                markers$marker)
for(i in seq_along(out)){
  this.marker <- names(out)[i]
  aics <- setNames(numeric(4), c('intercept', 'linear', 'quadratic', 'spline'))
  p <- 1
  for(type in c('intercept', 'linear', 'quadratic', 'spline')){
    aics[p] <- summary(TMBfit(this.marker, type))$AIC[1]
    p <- p + 1
  }
  out[[i]] <- aics
  message(this.marker, ' done')
}

lapply(out, which.min)


types <- c('intercept', 'linear', 'quadratic', 'spline')
out2 <- setNames(vector('list', length = length(markers$marker)),
                 markers$marker)
for(i in seq_along(out2)){
  this.marker <- names(out2)[i]
  p.val <- 0
  p <- 1
  last.fit <- TMBfit(this.marker, types[p])
  cat(paste0('Beginning model selection by ANOVA for ', this.marker, '.\n'))
  while(p.val < 0.05 && p <= 3){
    p <- p + 1
    next.fit <- tryCatch(TMBfit(this.marker, types[p]), warning = function(w) w)
    if(!inherits(next.fit, 'warning')){
      p.val <- anova(last.fit, next.fit)$`Pr(>Chisq)`[2]
      last.fit <- next.fit
      cat(paste0('Old type: ', types[p-1], ', new type: ', types[p], ', sig: ', round(p.val, 4),'.\n'))
    }else{
      cat(paste0('Model selection stopping for ', this.marker, '. TMB fit failed for ', types[p],'.\n'))
      p.val <- 1;p <- p-1
    }
  }
  cat('Recommending ', types[p], 'random + fixed specification for ', this.marker, '.\n\n')
  out2[[i]] <- types[p]
}


# Agreement?
unlist(lapply(out, function(x) names(x)[which.min(x)]))
unlist(out2)

fits <- setNames(vector('list', length = length(markers$marker)),
                 markers$marker)
for(i in seq_along(fits)){
  this.marker <- names(fits)[i]
  fits[[i]] <- TMBfit(this.marker, out2[[i]])
  message(this.marker)
}

lapply(fits, summary) # Interesting here there doesn't seem to always be sig. association with time (This backed-up somewhat in arXiv application).

checktimeSig <- function(f){ # Can verify this systematically
  cond <- summary(f)$coeff$cond
  time.vars <- grepl('time', row.names(cond))
  ps <- cond[time.vars, 'Pr(>|z|)']
  check <- ps < 5e-2
  length(which(check))
}

lapply(fits, checktimeSig)

