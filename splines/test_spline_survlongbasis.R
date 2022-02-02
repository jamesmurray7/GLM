dd <- data
sd <- data$survdat
sd$survtime

aa <- bs(sd$survtime, intercept = T, knots = quantile(sd$survtime, probs = c(.33, .66)))
bb <- bs(sd$survtime, intercept = F, knots = quantile(sd$survtime, probs = c(.33, .66)))

head(aa)
head(bb)

quantile(sd$survtime, probs = seq(0,1, length.out = 11))
quantile(sd$survtime %>% unique, probs = seq(0,1, length.out = 11))


ust <- unique(sd$survtime)
cut.ends <- function(vec){
  vl <- length(vec)
  vec[-c(1, vl)]
}

q10 <- quantile(ust, probs = seq(0, 1, length.out = 11))

cc <- cbind(sort(ust), bs(sort(ust), intercept = T, knots = cut.ends(q10)))
bb <- cbind(ust, bs(ust, intercept = F, knots = cut.ends(q10)))


head(round(aa, 5))
head(round(cc, 5))
head(round(bb, 5))

ust
sort(ust)

bs
cbind(0:8, bs(0:8))
cbind(0:8, bs(0:8, degree = 3))
image(cbind(0:8, bs(0:8, knots = c(0,2,4,6,8), degree = 3))[,2:9])
matplot(cbind(0:8, bs(0:8, degree = 3))[,2:4],
        type = 'l', xlab = '0:8')



a <- cbind(sort(ust), 
           bs(sort(ust), knots = cut.ends(q10)))
matplot(a[,-1], type = 'l', xaxt = 'n')
axis(at = seq(1, length(ust), length.out = 11), side = 1,
     labels = round(q10,2),  cex.axis = 1)

b <- cbind(sort(ust), 
           ns(sort(ust), knots = cut.ends(q10)))   #  so natural splines a bit less interpretable from this angle?
matplot(b[,-1], type = 'l', xaxt = 'n') 
axis(at = seq(1, length(ust), length.out = 11), side = 1,
     labels = round(q10,2),  cex.axis = 1)




# Basis on the actual event times -----------------------------------------
library(survival)
ph <- coxph(Surv(survtime, status) ~ cont  +bin, data$survdat)
bh <- basehaz(ph)
plot(bh[,1]~bh[,2], type = 's')
degree <- 3
.Tis<- unique(data$survdat$survtime) # Not sure if best practise to sort this...
.ft <- unique(data$survdat[data$survdat$status == 1, 'survtime'])
survbs <- bs(.Tis, degree = degree) # Not necessary to sort in bs(...)



# This is a matrix ordered by event times.
survbs.ft <- survbs[match(.ft, .Tis), ]    # and using this to obtain basis for failure times only.

