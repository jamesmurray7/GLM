rm(list=ls())
source('EM.R')
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



# log(serum bilirubin) ----------------------------------------------------
bil.pbc <- newpbc('serBilir')
pbc$serBilir <- log(pbc$serBilir)
EM(serBilir ~ drug + time + (1 + time|id),
   Surv(survtime, status) ~ drug,
   data = pbc, family = gaussian)
# To debug 21/4...