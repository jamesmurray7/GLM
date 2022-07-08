#' #####
#' longs.R
#' ---
#' Determining best longitudinal fits for each biomarker in PBC data
#' #####

rm(list=ls())
library(glmmTMB)
# library(tidyverse)
library(dplyr)
library(tidyr)
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
    geom_smooth(aes(group={{group}}), method = 'lm', colour = 'red', formula=y~splines::ns(x,knots=c(0.987077, 3.942723 ))) +
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
                spline = '~ drug * splines::ns(time, knots = c(0.987077, 3.942723))')
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
    aics[p] <- summary(TMBfit(this.marker, type, T))$AIC[1]
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


# Plot for paper ----------------------------------------------------------
library(forcats)
theme_csda <- function(base_family = "Arial", base_size = 12){
  theme_light(base_family = base_family, base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 8.5, colour = "black"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}


# Version 1: Failures only ------------------------------------------------
pbc %>% 
  filter(status == 1) %>% 
  pivot_longer(cols=serBilir:prothrombin, names_to='biomarker') %>% 
  mutate(biomarker = case_when(
    biomarker == 'serBilir' ~ "log(Serum~bilirubin)",
    biomarker == 'albumin' ~ "Albumin",
    biomarker == 'prothrombin' ~ '(0.1 ~ x ~ Prothrombin~time)^{-4}',
    biomarker == 'SGOT' ~ "log(AST)",
    biomarker == "platelets" ~ "Platelet~count",
    biomarker == 'alkaline' ~ "Alkaline~phosphatase",
    T ~ 'AA'
  ),
    tt = -1 * (survtime-time)
  ) %>% 
  mutate(f.biomarker = factor(biomarker, levels = c('log(Serum~bilirubin)',
                                                    'log(AST)', 'Albumin', '(0.1 ~ x ~ Prothrombin~time)^{-4}',
                                                    'Platelet~count', "Alkaline~phosphatase"))) %>% 
  filter(biomarker != 'AA') %>% 
  ggplot(aes(x=tt, y = value, group = id)) + 
  geom_vline(xintercept = 0, colour = 'black', alpha = .25) + 
  geom_line(alpha = .10) + 
  geom_smooth(aes(group=NULL), colour = 'black', method = 'loess', formula = y~x) +
  # geom_smooth(aes(group=NULL), colour = 'black', method = 'lm', formula = y~splines::ns(x, 3))+
  # geom_smooth(aes(group=NULL), colour = 'red', method = 'lm', formula = y~x) + 
  # geom_smooth(aes(group=NULL), colour = 'blue', method = 'lm', formula = y~x+I(x^2)) + 
  facet_wrap(~f.biomarker, scales = 'free', strip.position = 'left', labeller = label_parsed) + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  labs(y = NULL,
       x = 'Time (years) from death (0: time of death)') + 
  theme_csda() + 
  theme(strip.placement = 'outside',
        strip.text = element_text(vjust = 1))
#ggsave('~/Downloads/PBCtrajectories.png', width = 140, height = 90, units = 'mm')

pbc %>% 
  pivot_longer(cols=serBilir:prothrombin, names_to='biomarker') %>% 
  mutate(biomarker = case_when(
    biomarker == 'serBilir' ~ "log(Serum~bilirubin)",
    biomarker == 'albumin' ~ "Albumin",
    biomarker == 'prothrombin' ~ '(0.1 ~ x ~ Prothrombin~time)^{-4}',
    biomarker == 'SGOT' ~ "log(AST)",
    biomarker == "platelets" ~ "Platelet~count",
    biomarker == 'alkaline' ~ "Alkaline~phosphatase",
    T ~ 'AA'
  )) %>% 
  mutate(f.biomarker = factor(biomarker, levels = c('log(Serum~bilirubin)',
                                                    'log(AST)', 'Albumin', '(0.1 ~ x ~ Prothrombin~time)^{-4}',
                                                    'Platelet~count', "Alkaline~phosphatase"))) %>% 
  filter(biomarker != 'AA') %>% 
  ggplot(aes(x=time, y = value, group = id,  as.factor(drug))) + 
  geom_vline(xintercept = 0, colour = 'black', alpha = .25) + 
  geom_line(alpha = .10) + 
  # geom_smooth(aes(group=NULL), colour = 'black', method = 'loess', formula = y~x) +
  geom_smooth(aes(group=NULL), colour = 'black', method = 'lm', formula = y~splines::ns(x, df = 3))+
  # geom_smooth(aes(group=NULL), colour = 'red', method = 'lm', formula = y~x) + 
  # geom_smooth(aes(group=NULL), colour = 'blue', method = 'lm', formula = y~x+I(x^2)) + 
  facet_wrap(~f.biomarker, scales = 'free', strip.position = 'left', labeller = label_parsed) + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  labs(y = NULL,
       x = 'Time (years) from death (0: time of death)') + 
  theme_csda() + 
  theme(strip.placement = 'outside',
        strip.text = element_text(vjust = 1))
