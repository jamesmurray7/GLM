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
# makeFormula <- function(y, type = 'linear'){
#   rhs <- switch(type,
#                 linear = ' ~ drug * + (1 + time|id)',
#                 intercept = '~ drug * time + (1|id)',
#                 quadratic = '~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id)',
#                 spline = '~ drug * splines::ns(time, df = 3) + (1 + splines::ns(time, df = 3)|id)')
#   as.formula(paste0(y, rhs))
# }

makeFormula <- function(y, type = 'linear', drug = T){
  RE <- '+ (1|id)' # Model selection based on the fixed effects part only. I think REs are influencing heavily the conclusion one draws from AIC
  mid <- switch(type,
                linear = ' ~ drug * time',
                quadratic = '~ drug * (time + I(time^2))',
                spline = '~ drug * splines::ns(time, df = 3)')
  if(!drug) mid <- gsub('drug\\s\\*\\s', '', mid)
  as.formula(paste0(y, mid, RE))
}

TMBfit <- function(y, type, drug){
  form <- makeFormula(y, type, drug)
  family <- markers[markers$marker == y, 'family']
  fit <- glmmTMB(form, data = pbc, family = family)
  fit
}


out <- setNames(vector('list', length = length(markers$marker)),
                markers$marker)
for(i in seq_along(out)){
  this.marker <- names(out)[i]
  aics <- setNames(numeric(3), c('linear', 'quadratic', 'spline'))
  p <- 1
  for(type in c('linear', 'quadratic', 'spline')){
    aics[p] <- summary(TMBfit(this.marker, type))$AIC[1]
    p <- p + 1
  }
  out[[i]] <- aics
  message(this.marker, ' done')
}

lapply(out, which.min)


types <- c('linear', 'quadratic', 'spline')
out2 <- setNames(vector('list', length = length(markers$marker)),
                 markers$marker)

drug <- T #set T/F
for(i in seq_along(out2)){
  this.marker <- names(out2)[i]
  p.val <- 0
  message(this.marker)
  linr <- TMBfit(this.marker, 'linear', drug)
  quad <- TMBfit(this.marker, 'quadratic', drug)
  p <- anova(linr, quad)$`Pr(>Chisq)`[2]
  cat('linear -> quadratic p.val =', p, '\n')
  if(p < 0.05){
    spli <- TMBfit(this.marker, 'spline', drug)
    p <- anova(quad, spli)$`Pr(>Chisq)`[2]
    cat('quadratic -> spline p.val = ', p, '\n')
    if(p < 0.05){
      cat('Recommending spline fit for ', this.marker, '.\n', sep='')
    }else{
      cat('Recommending quadratic fit for ', this.marker, '.\n', sep='')
    }
  }else{
    cat('Recommending linear fit for ', this.marker, '.\n', sep='')
  }
}


# Plot competing for cont/counts ------------------------------------------
plotCompeting <- function(marker, type1, type2, drug = T){
  fit1 <- TMBfit(marker, type1, drug)
  fit2 <- TMBfit(marker, type2, drug)
  newData <- pbc[!is.na(pbc[, marker]), ]
  newData$pred1 <- fitted(fit1, 'response')
  newData$pred2 <- fitted(fit2, 'response')
  newData$y <- newData[,marker]
  
  ggplot(newData, aes(x = time, y = y, group = id)) + 
    geom_line(alpha = .25) + 
    geom_smooth(aes(y = pred1, group = NULL), method = 'loess', formula = y ~ x, colour = 'red') + 
    geom_smooth(aes(y = pred2, group = NULL), method = 'loess', formula = y ~ x, colour ='blue') + 
    theme_light()
}

plotCompeting('serBilir', 'quadratic', 'spline')
plotCompeting('SGOT', 'linear', 'quadratic')
plotCompeting('albumin', 'quadratic', 'spline')
plotCompeting('prothrombin', 'linear', 'quadratic')
plotCompeting('platelets', 'quadratic', 'spline')
plotCompeting('alkaline', 'quadratic', 'spline')
