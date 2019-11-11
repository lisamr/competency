rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(scales)
library(rethinking)
library(brms)
library(tidybayes)
#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))

#figure out how to deal with offsets. Simulate data and confirm with model.
TC <- c(rpois(30, 200),  rpois(30, 20)) #true counts
aID <- rnorm(60, 1, 1) #random observation level noise
prop <- c(rep(.05, 15), rep(.2, 15), rep(1,30)) #proportion sampled
counts <- rpois(length(TC), TC*prop) #sampled counts
sp <- c(rep("ACMA", 30), rep("LIDE", 30))
fake <- data.frame(sp, leafID=1:length(sp), prop, counts, TC)

#summary of fake data. means of sampled counts are the same, but higher sd 
fake %>% group_by(sp) %>% summarise(mean(counts/prop), sd(counts/prop), mean(TC), sd(TC))
ggplot(fake, aes(sp, TC)) +
  geom_boxplot()
ggplot(fake, aes(sp, counts/prop)) +
  geom_boxplot()
ggplot(fake, aes(sp, counts)) +
  geom_boxplot()+
  geom_violin()

################################
#brms model
#first one for distinguishing species
#model it
m1 <- brm(
  counts ~ sp + offset(log(prop)),
  data = fake, family = poisson(),
  chains = 3, cores = 4, sample_prior = T,
  control = list(adapt_delta = .95, max_treedepth=15))
summary(m1)
plot(m1, pars = "^b_")#looks good

#check out model fit
pp_check(m1, nsamples = 30)
msim <- fake %>% 
  select(sp, perc_lesion, counts, prop) %>% 
  add_predicted_draws(m1) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(sp, .draw) %>% 
  sample_draws(50) 
msim %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.5, color='grey')+
  stat_density(data=fake, aes(x=counts), geom="line", color='slateblue')+
  facet_wrap(~sp, scales = 'free')

#get predicted means. yep recovers alpha well.
pr <- fake %>%
  modelr::data_grid(sp) %>% mutate(prop=1)%>%
  add_fitted_draws(m1,scale = 'response') 
pr %>% median_qi(.width=c(.9, .95))

##################
#add in lesion area data
beta <- rep(c(5, .1), each=30)
dens(rgamma(1000, shape = 1, scale = 1), xlim=c(0,10))#see how gamma distributions work. scale is the mu and shape is a multiplier
beta2 <- rgamma(length(beta), shape = 1, scale=beta)#variation around the beta
fake$perc_lesion <- log(rpois(length(fake$TC), fake$TC))/beta2
fake$perc_lesion[fake$perc_lesion>100] <- 100 #cap it at 100
ggplot(fake, aes(perc_lesion, TC)) +
  geom_point(aes(color=sp)) + xlim(0,100)

#then another one to look at spores vs lesion size
fm2 <- bf(counts ~ perc_lesion + (perc_lesion|sp) + offset(log(prop)), family = poisson())

#need to set priors
#check them out visually
get_prior(fm2, fake)
a <- rnorm(1000, 1, 2.5)#prior for int
dens(exp(a), xlim=c(0, 100))
b <- rnorm(5000, 0, .1)#prior for slope 
sp <- seq(0,100, length.out = 100)
lam <- sapply(1:length(sp), function(x) exp(a + b*sp[x]))
plot(NULL, xlim=c(0,100), ylim=c(0,250), xlab="lesion", ylab="#spores")
lapply(1:200, function(x) lines(sp, lam[x,], col=alpha(1, .3)))
dens(rlkjcorr(1000, 2, 2)[,,1][,2])#make correlation center around 0, not uniform

prior1 <- c(
  set_prior("normal(0,.1)", class = "b", coef='perc_lesion'),
  set_prior("normal(1, 2.5)", class = "Intercept"),
  set_prior("exponential(1)", class = "sd"),
  set_prior("lkj(2)", class='cor'))

#model it
m2 <- brm(
  fm2, data = fake, prior = prior1,
  chains = 3, cores = 4, sample_prior = T,
  control = list(adapt_delta = .95, max_treedepth=15))

summary(m2)
plot(m2)#had to narrow priors
pp_check(m2)#not bad

#see beta
m2 %>% 
  gather_draws(r_sp[sp, term]) %>% 
  ggplot(aes(y = fct_rev(sp), x = .value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  facet_grid(~term, scales = 'free')

#check out model predictions
xx <- seq(0,100,length.out = 100)
newd <- data.frame(perc_lesion=xx, prop=1, sp=rep(unique(fake$sp), each=length(xx)))
pp1 <- fitted(m2, newdata = newd, summary = F, scale='response')
plotmodfit <- function(i, pp, data=fake, ...){
  sp <- as.vector(unique(newd$sp))
  use <- as.integer(factor(newd$sp))==i
  #estimate observed counts with offset
  data$counte <- data$counts/data$prop
  plot(data$perc_lesion[data$sp==sp[i]], data$counte[data$sp==sp[i]], ...,
       xlab='lesion', ylab='# chamydospores', main=sp[i], 
       pch=16, col=alpha('slateblue', .5), xlim=c(0,100), ylim=c(0,300))
  points(data$perc_lesion[data$sp==sp[i]], data$TC[data$sp==sp[i]], col=alpha('red', .2), pch=16)
  #sapply(1:100, function(x) lines(xx, pp[x,use], col=alpha('black', .1)))
  apply(pp[,use], 2, median) %>% lines(xx, .)
  apply(pp[,use], 2, PI, .9) %>% shade(., xx)
}
#plot predictions. blue are estimated counts, red are true counts, which actually line up with model predictions really well.
par(mfrow=c(1, 2))
sapply(1:length(unique(newd$sp)), function(x) plotmodfit(x, pp1))
sapply(1:length(unique(newd$sp)), function(x) plotmodfit(x, pp2))

################################
################################

#rethinking package

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
