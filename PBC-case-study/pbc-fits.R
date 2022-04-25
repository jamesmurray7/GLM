#' ###
#' Program which fits all univariate GLMM joint models to PBC data.
#' ###
rm(list=ls())
library(dplyr)
source('EM.R')
source('univariate-compare.R')
load('PBC.RData')


# Define a function that 'resets' the data --------------------------------
newpbc <- function(Y){
  y <- pbc[,Y]
  inds.to.remove <- is.na(y)
  if(length(inds.to.remove) > 0){
    rtn <- pbc[!inds.to.remove, ]
    rtn <- rtn %>% 
      group_by(id) %>% 
      mutate(new_id = cur_group_id()) %>% 
      ungroup %>% 
      select(-id) %>% 
      rename(id = new_id) %>% 
      as.data.frame()
  }else{
    rtn <- pbc
  }
  rtn
}

#' ####
#' GAUSSIANS
#' i.e. (log-)normal or Box-Cox transformed normals
#' ####

# log(serum bilirubin) ----------------------------------------------------
bil.pbc <- newpbc('serBilir')
bil.pbc$serBilir <- log(bil.pbc$serBilir)
my.serBilir <- EM(serBilir ~ drug + time + (1 + time|id),
                 Surv(survtime, status) ~ drug,
                 data = bil.pbc, family = gaussian, control = list(optimiser = 'optim', hessian = 'manual'))
joineR.serBilir <- joineR.fit(bil.pbc, 'serBilir', serBilir ~ drug + time, Surv(survtime, status) ~ drug,
                              c('time'), c('drug'))
joineRML.serBilir <- joineRML.fit(bil.pbc, serBilir ~ drug + time, ~ 1 + time|id, Surv(survtime, status) ~ drug)
my.summary(my.serBilir, 'log(Serum bilirubin)')
compare.gaussians(my.serBilir, joineR.serBilir, joineRML.serBilir, 'Serum Bilirubin')

JMbayes2.serBilir <- JMbayes2.fit(bil.pbc, gaussian, serBilir ~ drug + time, ~ 1 + time|id, Surv(survtime, status) ~ drug) # Don't run this on linux! Crashes

# Spline fits
long.formula <- serBilir ~ drug * splines::ns(time, knots = c(1, 4)) + 
                       (1 + splines::ns(time, knots = c(1, 4)) | id)
surv.formula <- Surv(survtime, status) ~ drug

my.serBilir <- EM(long.formula, surv.formula,
                  data = bil.pbc, family = gaussian, control = list(optimiser = 'optim', hessian = 'auto'))

joineR.serBilir <- joineR.fit(bil.pbc, 'serBilir', serBilir ~ drug + splines::ns(time,knots=c(1,4)), Surv(survtime, status) ~ drug,
                              c('time'), c('drug'))

joineRML.serBilir <- joineRML.fit(bil.pbc, serBilir ~ drug * splines::ns(time,knots=c(1,4)), ~ 1 + splines::ns(time,knots=c(1,4))|id, Surv(survtime, status) ~ drug)
summary(joineRML.serBilir)
my.summary(my.serBilir)


# log(serum aspartate aminotransferase) -----------------------------------
# AST, or aspartate aminotransferase, is one of the two liver enzymes. 
# It is also known as serum glutamic-oxaloacetic transaminase, or SGOT. AST is a protein made by liver cells. 
# When liver cells are damaged, AST leaks out into the bloodstream and the level of AST in the blood becomes elevated.

AST.pbc <- newpbc('SGOT')
AST.pbc$AST <- log(AST.pbc$SGOT)
my.AST <- EM(AST ~ drug + time + (1 + time|id),
             Surv(survtime, status) ~ drug,
             data = AST.pbc, family = gaussian)
my.summary(my.AST)
jML.AST <- joineRML.fit(AST.pbc, AST ~ drug + time, ~ 1 + time|id, Surv(survtime, status) ~ drug)
joineR.AST <- joineR.fit(AST.pbc, 'AST', AST ~ drug + time, Surv(survtime, status) ~ drug,
                              c('time'), c('drug'))
compare.gaussians(my.AST, joineR.AST, jML.AST, 'Aspartate aminotransferase')

# Serum albumin -----------------------------------------------------------
alb.pbc <- newpbc('albumin')
my.alb <- EM(albumin ~ drug + time + (1 + time|id),
             Surv(survtime, status) ~ drug,
             data = alb.pbc, family = gaussian)
my.summary(my.alb,  'Serum albumin')

jML.alb <- joineRML.fit(alb.pbc, albumin ~ drug + time, ~ 1 + time|id, Surv(survtime, status) ~ drug)
joineR.alb <- joineR.fit(alb.pbc, 'albumin', albumin ~ drug + time, Surv(survtime, status) ~ drug,
                         c('time'), c('drug'))
compare.gaussians(my.alb, joineR.alb, jML.alb, 'Serum albumin')


# (Transformed) prothrombin time ------------------------------------------
prot.pbc <- newpbc('prothrombin')
prot.pbc$prot <- (.1 * prot.pbc$prothrombin)^(-4)

my.prot <- EM(prot ~ drug + time + (1 + time|id),
              Surv(survtime, status) ~ drug,
              data = prot.pbc, family = gaussian)
my.summary(my.prot, 'Transformed prothrombin time')
jML.prot <- joineRML.fit(prot.pbc, prot ~ drug + time, ~ 1 + time|id, Surv(survtime, status) ~ drug)
joineR.prot <- joineR.fit(prot.pbc, 'prot', prot ~ drug + time, Surv(survtime, status) ~ drug,
                         c('time'), c('drug'))
compare.gaussians(my.prot, joineR.prot, jML.prot, 'Transformed prothrombin time')

#' ###
#' Counts
#' ###


# Platelet counts ---------------------------------------------------------
pla.pbc <- newpbc('platelets')
my.pla <- EM(platelets ~ drug + time + (1 + time|id),
             Surv(survtime, status) ~ drug,
             data = pla.pbc, family = poisson)
my.pla2 <- EM(platelets ~ drug + time + (1 + time|id),
             Surv(survtime, status) ~ drug, 
             data = pla.pbc, family = 'negative.binomial', control = list(hessian = 'auto'))
my.summary(my.pla, 'Platelet count')
my.summary(my.pla2, 'Platelet count (Negative binomial model)')

# Spline fit 
my.pla.spline <- EM(platelets ~ drug + splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
                    Surv(survtime, status) ~ drug,
                    data = pla.pbc, family = poisson, control = list(hessian = 'auto'))
my.summary(my.pla.spline, 'Platelet count (poisson w/ natural spline)')
my.pla.spline2 <- EM(platelets ~ drug + splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id),
                    Surv(survtime, status) ~ drug,
                    data = pla.pbc, family = 'negative.binomial', control = list(hessian = 'auto'))
my.summary(my.pla.spline2, 'Platelet count (Negative binomial w/ natural spline)')

#' JMbayes2 fit
# No splines
JMb.pla.pois <- JMbayes2.fit(pla.pbc, poisson, # Does not converge
                             platelets ~ drug + splines::ns(time, knots = c(1, 4)), ~ 1 + splines::ns(time, knots = c(1, 4))|id,
                             Surv(survtime, status) ~ drug)

# Alkaline Phosphates -----------------------------------------------------
alk.pbc <- newpbc('alkaline')
my.alk <- EM(alkaline ~ drug + time + (1 + time|id),
             Surv(survtime, status) ~ drug,
             data = alk.pbc, family = poisson, control = list(hessian = 'auto'))
my.alk2 <- EM(alkaline ~ drug + time + (1 + time|id),
              Surv(survtime, status) ~ drug,
              data = alk.pbc, family = 'negative.binomial', control = list(hessian = 'auto'))
my.summary(my.alk,  'Alkaline Phosphates')
my.summary(my.alk2, 'Alkaline Phosphates (Negative binomial model)') # Is sig. with negbin?

#' TO DO: JMbayes2 fit
JMb.alk <- JMbayes2.fit(alk.pbc, poisson,
                        alkaline ~ drug + time, ~ 1+time|id, Surv(survtime, status) ~ drug)

#' ###
#' Binomial
#' ###


# Spiders -----------------------------------------------------------------
# (Presence of blood vessel malformations in the skin.)
spi.pbc <- newpbc('spiders')
my.spi <- EM(spiders ~ drug * time + (1 + time|id),
             Surv(survtime, status) ~ drug,
             data = spi.pbc, family = binomial, control = list(hessian = 'auto'))
my.summary(my.spi, 'Spiders (binomial)')

#' TO DO: JMbayes2 fit
jm.spi <- JMbayes2.fit(spi.pbc,
                       binomial,
                       spiders ~ drug * time,
                       ~ 1 + time | id,
                       Surv(survtime, status) ~ drug)

my.summary(my.spi, 'Spiders (binomial)')
summary(jm.spi)

# Ascites -----------------------------------------------------------------
# Presence of abnormal accumulation of fluid in abdomen.
asc.pbc <- newpbc('ascites')
my.asc <- EM(ascites ~ drug * time + (1 + time|id),
             Surv(survtime, status) ~ drug,
             data = asc.pbc, family = binomial, control = list(hessian = 'auto'))
my.summary(my.asc, 'Ascites (binomial)')

#' JMbayes2 fit
jm.asc <- JMbayes2.fit(asc.pbc,
                       binomial,
                       ascites ~ drug * time,
                       ~ 1 + time | id,
                       Surv(survtime, status) ~ drug)

my.summary(my.asc, 'Ascites (binomial)')
summary(jm.asc)

# Hepatomegaly ------------------------------------------------------------
# Enlarged liver
hep.pbc <- newpbc('hepatomegaly')
my.hep <-  EM(hepatomegaly ~ drug * time + (1 + time|id),
              Surv(survtime, status) ~ drug,
              data = hep.pbc, family = binomial, control = list(hessian = 'auto'))
my.summary(my.hep, 'Hepatomegaly (binomial)')

#' JMbayes2 fit
jm.hep <- JMbayes2.fit(hep.pbc,
                       binomial,
                       hepatomegaly ~ drug * time,
                       ~ 1 + time | id,
                       Surv(survtime, status) ~ drug)

my.summary(my.hep, 'Hepatomegaly (binomial)')
summary(jm.hep)
