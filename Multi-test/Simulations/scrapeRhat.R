library(tidyverse)
data.dir <- paste0(getwd(), '/Simulations/fits/')
files <- as.list(dir(data.dir, pattern = 'JMb'))
to.summarise <- lapply(files, function(x){
  assign('fits', get(load(paste0(data.dir,x))))
  out <- lapply(fits, function(y){
    df <- as.data.frame(
      rbind(y$Outcome1,
            y$Outcome2,
            y$Outcome3,
            y$survival)
    )
    df %>% rownames_to_column('parameter') %>% 
      select(parameter, Rhat) %>% 
      mutate(fit = x)
  })
  do.call(rbind, out)
})

lapply(to.summarise, function(x){
  cat(x$fit[1],'\n')
  # with(x, tapply(Rhat, parameter, range))
  print(range(x$Rhat))
})
