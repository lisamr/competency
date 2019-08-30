#contrast chlamydospore counts, similar to what you did for sporangia. will have more zeros and might need to do a zeroinflated model?

setwd("~/Box/Competency project/competency.git/data2019")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)
library(scales)

#read in master file
df <- read.csv(file = 'master_tall.csv')
#analyze just broad leaf assay
df1 <- df %>% filter(leafID<=530) 

###########################################
#DATA VIZ
#map the plants
df1 %>% filter(spore_assay=="S", trt=="T") %>%
  ggplot(aes(easting, northing)) +
  geom_point(aes(color=species))+
  scale_color_viridis_d()
#plot chlamydospore counts on a map for select species
df1 %>% filter(spore_assay=="C", trt=="T", species=="CEOL") %>%
  ggplot(aes(easting, northing)) +
  geom_point(aes(color=count1))+
  scale_color_viridis_c()

#let's try to see how chlamydo counts differ across species
#2 choices: estimate # of spores by dividing counts by proportion sampled OR (zero inflated to handle zeros and gamma distribution?) OR use a count dist with offsets (poisson)

#start with averages
df2<- df1 %>% mutate(count_est=count1/prop) %>% filter(spore_assay=="C", !is.na(count_est), trt=="T")

#plot it. well they aren't all that different!
df2 %>% 
  ggplot( aes(species, (count_est)))+
  geom_boxplot()+
  geom_point(alpha=.4)

#there's a fair amount of zeros. I wonder if it's warranted to do a zero inflated model. anyways, viz omitting <1.
df2 %>% filter( !count_est<1) %>% 
  ggplot( aes(species, log(count_est)))+
  geom_boxplot()+
  geom_point(alpha=.4)


#chlamydos vs lesion size?
df2 %>% 
  ggplot( aes(count_est, perc_lesion))+
  geom_point(alpha=.4)+
  facet_wrap(~species, scales = 'free_x')

####################################
#statistical model to compare chlamydos across species.
#need a poisson glm with sampling intensity included as offset 
#might need a mixture model to deal with zero inflation (mostly from todi and lide)

library(rethinking)

head(df1)
#create data for stan
dat <- df1 %>% 
  filter(spore_assay=="C", trt=="T", !is.na(count1)) %>% 
  mutate(logprop=log(prop), 
         spID=as.integer(factor(species))) %>% 
  select(spID, count=count1, logprop) %>% 
  as.list()

#figure out prior for the species intercept
curve( dlnorm( x, 7, 1.5) , from=0 , to=5000 , n=200, xlab = "mean # spores (lambda)")

#model it
m1 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- logprop + a[spID], 
    #priors
    a[spID] ~ dnorm(7, 1.5)
  ), data=dat, chains=3
)

traceplot(m1) #looks good
par(mfrow=c(1,1))
stancode(m1)#check out stan code

precis(m1, depth = 2) #it sampled very well

#check out model predictions
fplot <- function(model, dat){
  #check out model prediction vs observed
  dpred <- data.frame(spID=1:6, logprop=0)
  #CI for the model fit
  fit <- link(model, dpred)
  fmean <- apply(fit, 2, mean)
  fpi <- apply(fit, 2, PI, .95) %>% round(1) %>%  t %>% as.data.frame() %>% rename(mulow='3%', muhigh='98%')
  #CI for predicted values
  p <- sim(model, dpred)
  ppi <- apply(p, 2, PI, .95) %>% round(1) %>%  t %>% as.data.frame() %>% rename(prlow='3%', prhigh='98%')
  #put into df
  spp <- df1 %>% filter(spore_assay=="C", trt=="T", !is.na(count1)) %>% pull(species) %>% unique()
  pred <- data.frame(species=spp, fmean, fpi, ppi)
  #get observed df
  c.est <- dat$count/(exp(dat$logprop))
  obs <- data.frame(dat$spID, species=spp[dat$spID], c.est) 
  
  #plot
  ggplot(pred, aes(species, fmean))+
    geom_boxplot(data=obs, aes(species, c.est))+
    geom_point(alpha=.5) +
    geom_errorbar(aes(ymin=mulow, ymax=muhigh), width=.1, color="red") +
    geom_errorbar(aes(ymin=prlow, ymax=prhigh), width=.1, color="red")
}

fplot(m1, dat)

################################
################################
postcheck(m1)
################################
################################
#predict using the same data frame? that way dont need to estimate counts.
simm1 <- sim(m1)

viz <- function(i, fit, data, ...){
  dens(data$count[data$spID==i], main=paste("species",i), col="blue", ...)
  fit <- fit[,data$spID==i]
  sapply(1:200, function(x) dens(fit[x,], add=T, col=alpha("black", .05)))
}
#ugggghhh theyre so bad.
lapply(1:6, function(x) viz(x, simm1, dat)) 

################################
################################
#maybe negative binomial model???? adds phi as a dispersion par
#model it
m2 <- ulam(
  alist(
    count ~ dgampois(lambda, phi), #likelihood
    log(lambda) <- logprop + a[spID], 
    #priors
    a[spID] ~ dnorm(7, 1.5),
    phi ~ dexp(1)
  ), data=dat, chains=3
)
precis(m2, depth = 2)
postcheck(m2) #definitely better, but I think I need to do a zero-infl model

#check out model predictions
simm2 <- sim(m2)
lapply(1:6, function(x) viz(x, simm2, dat)) 



#better than the poisson, but still needs work
fplot(m1, dat)
fplot(m2, dat)


#making phi species specific. nope can't do it.

################################
################################
#how about a zero-inflated negative binomial model???? lets start with zi-pois since there's already a helper function for that distribution. this better work.

#does removing hear help? no, makes no diff. don't bother.
#gonna remove HEAR becuase its screwing with the model. dont know how to keep it in.
rm <- dat$spID==4
dat2 <- lapply(1:3, function(x) dat[[x]][!rm])
names(dat2) <- names(dat)
dat2$spID <- as.integer(factor(dat2$spID))


#check out priors. include the distribution and then backtransform it to see it on the outcome scale.
p <- inv_logit(rnorm(1000, -1, 1.5))
dens(p, xlim=c(0,1)) #prob of no chl present
min(p);max(p)
lam <- exp(rnorm(1000, 6, 1.5))
dens(lam, xlim=c(0,5000)) #mean chl when infected
mean(lam)
min(lam);max(lam)

#model it
m3 <- ulam(
  alist(
    count ~ dzipois(p, lambda), #likelihood
    logit(p) <- ap[spID], #probability of no chl produced
    log(lambda) <- logprop + al[spID], #mean spores prod
    #priors
    ap[spID] ~ dnorm(-1, 1.5),
    al[spID] ~ dnorm(6, 1.5)
  ), data=dat2, chains=3
)

precis(m3, depth = 2)

fplot2 <- function(model){
  #check out model prediction vs observed
  dpred <- data.frame(spID=1:5, logprop=0)
  #CI for the model fit
  fit <- link(model, dpred)
  fmean <- apply(fit$lambda, 2, mean)
  fpi <- apply(fit$lambda, 2, PI, .95) %>% round(1) %>%  t %>% as.data.frame() %>% rename(mulow='3%', muhigh='98%')
  #CI for predicted values
  p <- sim(model, dpred)
  ppi <- apply(p, 2, PI, .95) %>% round(1) %>%  t %>% as.data.frame() %>% rename(prlow='3%', prhigh='98%')
  #put into df
  spp <- c("ACMA", "ARME", "CEOL", "LIDE", "TODI")
  pred <- data.frame(species=spp, fmean, fpi, ppi)
  #get observed df
  c.est <- dat2$count/(exp(dat2$logprop))
  obs <- data.frame(dat2$spID, species=spp[dat2$spID], c.est) 
  
  #plot
  ggplot(pred, aes(species, fmean))+
    geom_boxplot(data=obs, aes(species, c.est))+
    geom_point(alpha=.5) +
    geom_errorbar(aes(ymin=mulow, ymax=muhigh), width=.1, color="red") +
    geom_errorbar(aes(ymin=prlow, ymax=prhigh), width=.1, color="red")
}


#check out model predictions
fplot2(m3) #estimated counts
simm3 <- sim(m3) #same data
lapply(1:5, function(x) viz(x, simm3, dat2)) 
spp
par(mfrow=c(3, 6))
viz(1, simm2, dat, xlim=c(0,500), ylim=c(0,.02))#acma
viz(2, simm2, dat, xlim=c(0,100), ylim=c(0,.03))#arme
viz(3, simm2, dat, xlim=c(0,650), ylim=c(0,.02))#ceol
viz(4, simm2, dat, xlim=c(0,10), ylim=c(0,.1))#hear
viz(5, simm2, dat, xlim=c(0,100), ylim=c(0,.1))#lide
viz(6, simm2, dat, xlim=c(0,200), ylim=c(0,.1))#lide
viz(1, simm1, dat, xlim=c(0,500), ylim=c(0,.02))#acma
viz(2, simm1, dat, xlim=c(0,100), ylim=c(0,.03))#arme
viz(3, simm1, dat, xlim=c(0,650), ylim=c(0,.02))#ceol
viz(4, simm1, dat, xlim=c(0,10), ylim=c(0,.1))#hear
viz(5, simm1, dat, xlim=c(0,100), ylim=c(0,.1))#lide
viz(6, simm1, dat, xlim=c(0,200), ylim=c(0,.1))#lide

viz(1, simm3, dat2, xlim=c(0,500), ylim=c(0,.02))#acma
viz(2, simm3, dat2, xlim=c(0,100), ylim=c(0,.03))#arme
viz(3, simm3, dat2, xlim=c(0,650), ylim=c(0,.02))#ceol
viz(4, simm3, dat2, xlim=c(0,100), ylim=c(0,.1))#lide
viz(4, simm3, dat2, xlim=c(0,100), ylim=c(0,.1))#lide
viz(5, simm3, dat2, xlim=c(0,200), ylim=c(0,.1))#todi

################################
################################
#need to do a zi-negbinom model. ughhhggghghg.

