rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(scales)
library(rethinking)
#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))

#figure out how to deal with offsets. Simulate data and confirm with model.
TC <- c(rpois(30, 200),  rpois(30, 20)) #true counts
prop <- c(rep(.05, 15), rep(.2, 15), rep(1,30)) #proportion sampled
counts <- rpois(length(TC), TC*prop) #sampled counts
sp <- c(rep("ACMA", 30), rep("LIDE", 30))
fake <- data.frame(sp, prop, counts, TC)
#summary of fake data. means of sampled counts are the same, but higher sd 
fake %>% group_by(sp) %>% summarise(mean(counts/prop), sd(counts/prop), mean(TC), sd(TC))
ggplot(fake, aes(sp, TC)) +
  geom_boxplot()
ggplot(fake, aes(sp, counts/prop)) +
  geom_boxplot()
ggplot(fake, aes(sp, counts)) +
  geom_boxplot()+
  geom_violin()

#build a model
fake2 <- list(sp=as.integer(fake$sp), prop=fake$prop, count=fake$counts)
m0 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- log(prop) + a[sp], 
    #priors
    a[sp] ~ dnorm(5, 1)
  ), data=fake2, chains=3
)

#not the right one. you do want to just log(prop) as the offset. wierdly the model seemed to fit the data really well, which is concerning. how do I adequately check the veracity of the model?
m0.1 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- log(1/prop) + a[sp], 
    #priors
    a[sp] ~ dnorm(5, 1)
  ), data=fake2, chains=3
)

#check out model performance. how well does it fit the data? pretty well.
model=m0
simd <- sim(model)
sp1 <- 1:30
sp2 <- 31:60
par(mfrow=c(1,2))
#sp1. half samples have different offsets.

dens(fake2$count[sp1], col=NA, ylim=c(0,.2))
for(i in 1:200){
  dens(simd[i,sp1], add=T, col=alpha('black', .2))
}
dens(fake2$count[sp1], add=T, col='blue')
#sp2
dens(fake2$count[sp2], col=NA, ylim=c(0,.2))
for(i in 1:200){
  dens(simd[i,sp2], add=T, col=alpha('black', .2))
}
dens(fake2$count[sp2], add=T, col='blue')

#check out model. did the offset do what you want? yeah it did.
traceplot(model)
dev.off()
precis(model, depth=2)
pd <- data.frame(sp=c(1,2), prop=1)
postraw <- link(model, pd)
postm <- apply(postraw, 2, mean)
postsd <- apply(postraw, 2, sd)
#make summary table with obs and pred
fake %>% group_by(sp) %>% 
  summarise(obsm=mean(TC), obssd=sd(TC)) %>% 
  mutate(pm=postm, psd=postsd)
