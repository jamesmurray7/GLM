setwd('~/Documents/GLMM/mix/_01')
load('pbc.RData') # i think this is ascites?
pbcdata

source('EM.R')

data <- pbcdata$pbc
survdata <- pbcdata$survdata
ph <- coxph(Surv(survtime, status) ~ cont + bin, survdata)

fit <- EM(data, ph, survdata, gh = 3)
fit

source('dynSurv/draw.R')
source('dynSurv/prepdata.R')
source('dynSurv/dynSurv.R')
sourceCpp('dynSurv/helper.cpp')

data[data$id == 1, ]
data %>% filter(status == 1) %>% distinct(survtime) %>% max
args(dynSurv)
ds <- dynSurv(fit, data, id = 1, u = c(2,3,4,5,6,7,8,9))
ds <- do.call(rbind, ds)
ds

plot(ds[,1] ~ c(2:9),
     ylim = c(min(pmin(ds[, 2], ds[, 3])), max(pmax(ds[,2], ds[,3]))),
     type = 'o', col = 'red',
     xlab = 'u', ylab = 'S')
lines(ds[,2] ~ c(2:9), lty = 5)
lines(ds[,3] ~ c(2:9), lty = 5)

# bit of a longer profile

data[data$status == 1, ]
data %>% filter(status == 1) %>% distinct(survtime) %>% max
args(dynSurv)
ds <- dynSurv(fit, data, id = 8, u = c(7,8,9,10,11))
ds <- do.call(rbind, ds)

plot(ds[,1] ~ c(7:11),
     ylim = c(min(pmin(ds[, 2], ds[, 3])), max(pmax(ds[,2], ds[,3]))),
     type = 'o', col = 'red',
     xlab = 'u', ylab = 'S')
lines(ds[,2] ~ c(7:11), lty = 5)
lines(ds[,3] ~ c(7:11), lty = 5)

# Someone who survives

data[data$status == 0, ]
ds <- dynSurv(fit, data, id = 2, u = c(9, 10, 11, 12, 13))
ds <- do.call(rbind, ds)
plot(ds[,1] ~ c(9:13),             # So model would expect this person to fail?
     ylim = c(min(pmin(ds[, 2], ds[, 3])), max(pmax(ds[,2], ds[,3]))),
     type = 'o', col = 'red',
     xlab = 'u', ylab = 'S')
lines(ds[,2] ~ c(9:13), lty = 5)
lines(ds[,3] ~ c(9:13), lty = 5)

# Sample of some who survive
data %>% 
  filter(status == 0, id != 2) %>% 
  group_by(id) %>% 
  mutate(r = row_number(),
         r = max(r)) %>% ungroup %>% 
  arrange(-r) %>% distinct(id, r) %>% head(6) %>% pull(id) -> ids

ft <- surv.mod(ph, data)$ft
dynsurvs <- list()
for(i in ids){
  survtime <- unique(data[data$id == i, 'survtime'])
  project.times <- ft[which(ft > survtime)]
  times <- unique(data[data$id == i, 'time'][-c(1,2)])
  times <- times[which(times < max(ft))]
  uu <- union(times, project.times)
  print(uu)
  ds <- dynSurv(fit, data, id = i, u = uu, nsim = 100)
  dynsurvs[[i]] <- list(ds = do.call(rbind, ds),
                        survtime = survtime,
                        u = uu)
}

pdf('~/Downloads/Survived.pdf')
par(mfrow = c(3, 2))
lapply(dynsurvs,
       function(x) if(!is.null(x)){
         S <- x$ds[, 1]; lb <- x$ds[, 2]; ub <- x$ds[, 3]
         plot(S ~ x$u, type = 'o', col = 'red', ylab = 'S', xlab = 'u',
              ylim = c(min(pmin(lb, ub)),
                       max(pmax(lb, ub))),
              xlim = c(0, max(x$u, x$survtime)))
         lines(lb ~ x$u, lty = 5)
         lines(ub ~ x$u, lty = 5)
         abline(v = x$survtime, lty = 3, col = 'grey10')
       })
dev.off()


# Sample of those who did not survive -------------------------------------
data %>% 
  filter(status == 1) %>% 
  group_by(id) %>% 
  mutate(r = row_number(),
         m = max(r)) %>% ungroup %>% 
  filter(between(m, 4, 9)) %>% 
  distinct(id) %>% 
  pull %>% 
  sample(., 6, T) -> ids
  
dynsurvs <- list()
p <- 1
for(i in ids){
  survtime <- unique(data[data$id == i, 'survtime'])
  project.times <- ft[which(ft > survtime)]
  times <- unique(data[data$id == i, 'time'][-c(1,2)])
  times <- times[which(times < max(ft))]
  uu <- union(times, project.times)
  print(p)
  ds <- dynSurv(fit, data, id = i, u = uu, nsim = 100)
  dynsurvs[[p]] <- list(ds = do.call(rbind, ds),
                        survtime = survtime,
                        u = uu,
                        id = i)
  p <- p + 1
}


pdf('~/Downloads/Failed.pdf')
par(mfrow = c(3, 2))
lapply(dynsurvs,
       function(x) if(!is.null(x)){
         S <- x$ds[, 1]; lb <- x$ds[, 2]; ub <- x$ds[, 3]
         plot(S ~ x$u, type = 'o', col = 'red', ylab = 'S', xlab = 'u', pch = 20,
              ylim = c(min(pmin(lb, ub)),
                       max(pmax(lb, ub))),
              xlim = c(0, max(x$u, x$survtime)))
         lines(lb ~ x$u, lty = 3)
         lines(ub ~ x$u, lty = 3)
         abline(v = x$survtime, lty = 5, col = 'grey10')
       })
dev.off()



# Can i work out how to use geom_ribbon?
dynsurvs <- setNames(dynsurvs, paste0('id = ', ids))
all.ds <- as.data.frame(do.call(rbind, lapply(dynsurvs, el, 1)))
all.ds <- all.ds %>% 
  rownames_to_column('u') %>% 
  mutate(u = str_remove_all(u, '^u\\.\\.\\.'),
         u = str_remove_all(u, '\\.\\d$'),
         u = as.numeric(u))

all.ids <- unname(do.call(c, lapply(dynsurvs, function(x) rep(x$id, nrow(x$ds)))))
all.sts <- unname(do.call(c, lapply(dynsurvs, function(x) rep(x$survtime, nrow(x$ds)))))

all.ds$ids <- all.ids; all.ds$survtime <- all.sts

all.ds %>% mutate(id = paste0('id = ', ids)) %>% 
  ggplot(aes(x = u, y = `50%`)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = 'grey80') +
  geom_vline(aes(xintercept = survtime), lty = 5) + 
  geom_line(colour = 'red') + 
  geom_point(pch = 19, colour = 'red') + 
  facet_wrap(~id, scales = 'free') + 
  scale_x_continuous(breaks = 1:14)
  theme_light()


