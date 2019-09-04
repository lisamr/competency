#contrast chlamydospore counts, similar to what you did for sporangia. will have more zeros and might need to do a zeroinflated model?

#use m2.1, which is a poisson model with random intercept for observation. good fit with the data, not great out of sample predictions. just use for inference with your data?

setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)
library(scales)
library(rethinking)
#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
options(scipen = 999) #turns off scientific notation

##################################################
#functions for plotting posterior preds

#check out posterior predictions against observed.
fdens <- function(s, sim, a=.2, ...){
  cols <- dat$spID==s
  dens(dat$count[cols], col=NA, main=spnames[s], ...)
  for(i in 1:200){
    dens(sim[i,cols], add=T, col=alpha('black', a))
  }
  dens(dat$count[cols], add=T, col='blue')
}

dpred <- data.frame(spID=1:6, prop=1)

#check out model predictions
fboxplot <- function(model, dat){
  #check out model prediction vs observed
  dpred <- data.frame(spID=1:6, prop=1)
  #CI for the model fit
  fit <- link(model, dpred)
  fmean <- apply(fit, 2, mean)
  fpi <- apply(fit, 2, PI, .95) %>% round(1) %>%  t %>% as.data.frame() %>% rename(mulow='3%', muhigh='98%')
  #CI for predicted values
  p <- sim(model, dpred)
  ppi <- apply(p, 2, PI, .95) %>% round(1) %>%  t %>% as.data.frame() %>% rename(prlow='3%', prhigh='98%')
  #put into df
  pred <- data.frame(species=spnames, fmean, fpi, ppi)
  #get observed df
  c.est <- dat$count/dat$prop
  obs <- data.frame(dat$spID, species=spnames[dat$spID], c.est)
  #plot
  ggplot(pred, aes(species, fmean))+
    geom_boxplot(data=obs, aes(species, c.est))+
    geom_point(alpha=.5) +
    geom_errorbar(aes(ymin=mulow, ymax=muhigh), width=.1, color="red") +
    geom_errorbar(aes(ymin=prlow, ymax=prhigh), width=.1, color="red") #+ scale_y_continuous(trans='log2')
}
########################################

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
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

head(df1)
#create data for stan
dat <- df1 %>% 
  filter(spore_assay=="C", trt=="T", !is.na(count1)) %>% 
  mutate(spID=as.integer(factor(species))) %>% 
  select(species, spID, count=count1, prop) %>% as.list()
spnames <- as.vector(unique(dat$species))
dat$species <- NULL

#quick try with glm. need something more snazzy.
datdf <- as.data.frame(dat)
datdf$countest <- datdf$count/datdf$prop
m0 <- glm(count ~ as.factor(spID) -1 + offset(log(prop)), data = datdf, family = poisson)
summary(m0)
plot(predict(m0, type='response'), datdf$countest)+abline(0,1)
m0.1 <- glm(count ~ as.factor(spID) -1 + offset(log(prop)), data = datdf, family = quasipoisson)
summary(m0.1)
plot(predict(m0.1, type='response'), datdf$countest)+abline(0,1)
m0.2 <- MASS::glm.nb(count ~ as.factor(spID) -1 + offset(log(prop)), data = datdf)
summary(m0.2)
plot(predict(m0.2, type='response'), datdf$countest)+abline(0,1)

#figure out prior for the species intercept
curve( dlnorm( x, 7, 1.5) , from=0 , to=5000 , n=200, xlab = "mean # spores (lambda)")

#model it
m1 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- log(prop) + a[spID], 
    #priors
    a[spID] ~ dnorm(7, 1.5)
  ), data=dat, chains=3
)

traceplot(m1) #looks good
par(mfrow=c(1,1))
stancode(m1)#check out stan code

precis(m1, depth = 2) #it sampled very well
postcheck(m1) #very poor fit. intervals way too small.
dev.off()

sim1 <- sim(m1)
#pdf('plots/chlamydos/poisson_pred_dens.pdf', width = 8, height = 6)
xmax <- c(.05, .1, .05, 4, .4, .2)
par(mfrow=c(2,3))
lapply(1:6, function(x) fdens(x, sim1, ylim=c(0, xmax[x])))
reset()
#dev.off()

fboxplot(m1, dat)

################################
################################
#maybe negative binomial model???? adds phi as a dispersion par
#model it
dens(rexp(10000, .1), xlim=c(-1, 100)) #prior for phi

m2 <- ulam(
  alist(
    count ~ dgampois(lambda, phi), #likelihood
    log(lambda) <- log(prop) + a[spID], 
    #priors
    a[spID] ~ dnorm(7, 1.5),
    phi ~ dexp(.1)
  ), data=dat, chains=3
)
precis(m2, depth = 2) %>% plot
postcheck(m2) #definitely better, but I think I need to do a zero-infl model
dev.off()

#check out model predictions
#pdf('plots/chlamydos/nbinom_pred_dens.pdf', width = 8, height = 6)
sim2 <- sim(m2)
par(mfrow=c(2,3))
lapply(1:6, function(x) fdens(x, sim2, ylim=c(0, xmax[x])))
reset()
#dev.off()

#better than the poisson, but still needs work
fboxplot(m2, dat)
#make phi species specific. nope can't do it.

################################
#poisson with random observation level variable?
dat2 <- df1 %>% 
  filter(spore_assay=="C", trt=="T", !is.na(count1)) %>% 
  mutate(spID=as.integer(factor(species)), 
         leafID=as.integer(factor(leafID))) %>% 
  select(spID, count=count1, prop, leafID) %>% as.list()

m2.1 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- log(prop) + a[spID] + b[leafID], 
    #priors
    a[spID] ~ dnorm(7, 1.5),
    b[leafID] ~ dnorm(0, sigmab),
    sigmab ~ dexp(1)
  ), data=dat2, chains=3
)
precis(m2.1, pars = c('a', 'sigmab'),  depth=2)#really bad sampling

#estimated species intercepts
precis(m2.1, pars='a', depth = 2) %>% plot

#WAAAYYY better fit! just have to fix the sampling now.
#pdf('plots/chlamydos/randobs_pois_pred_dens.pdf', width = 8, height = 6)
sim2.1 <- sim(m2.1)
par(mfrow=c(2,3))
lapply(1:6, function(x) fdens(x, sim2.1,a = .05, ylim=c(0, xmax[x])))
reset()
#dev.off()

#get posterior fit of lambda
ex <- extract.samples(m2.1)

#simulate posterior
ll <- log(1) + ex$a + 0 #no obs random effect
#ll <- log(1) + ex$a + rnorm(nrow(ex$a), 0, ex$sigmab)
post <- exp(ll)
postm <- apply(post, 2, mean)
postpi <- apply(post, 2, PI,.95)
obs <- dat %>% as.data.frame() %>% mutate(countest=count/prop) %>%  group_by(spID) %>% summarise(obsm=mean(countest))

#pdf('plots/chlamydos/randobs_pois_postmean.pdf', width = 6, height = 4)
postdf <- data.frame(spnames, spID=1:6, postm, t(postpi))
names(postdf)[c(4,5)] <- c('lower', 'upper')
ggplot(postdf, aes(spnames, postm))+
  geom_point()+
  geom_point(data=obs, aes(spID, obsm), color='blue') +
  geom_errorbar(aes(ymin=lower, ymax=upper, width=.1)) +
  scale_y_log10()
#dev.off()
#############################################
#predictions with the same data are good but not out of sample. maybe make random obs coef interact with species?

#there are real differences in the random obs coef depending on species.
f <- function(i){
  use <- dat$spID==i
  list(mean(ex$b[,use]), sd(ex$b[,use]))
}
sapply(1:6, f)


#model
rlkjcorr(10000,2, 2)[,2,1] %>% dens #rho prior

m3 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- log(prop) + a[spID] + b[leafID,spID], 
    #adaptive priors
    vector[6]:b[leafID] ~ multi_normal(0,Rho,sigmab),
    #fixed priors
    a[spID] ~ dnorm(7, 1.5),
    Rho ~ dlkjcorr(4),
    sigmab ~ dexp(1)
  ), data=dat2, chains=4, cores=4
)

#Warning messages:
#1: There were 13 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
#2: There were 4 chains where the estimated Bayesian Fraction of Missing Information was low. See http://mc-stan.org/misc/warnings.html#bfmi-low 
#3: Examine the pairs() plot to diagnose sampling problems

precis(m3, depth=2)
precis(m3, pars=c('a', 'sigmab'), depth=3)

#predictions similar to m2.1
#pdf('plots/chlamydos/randobs_pois_pred_dens.pdf', width = 8, height = 6)
sim3 <- sim(m3)
par(mfrow=c(2,3))
lapply(1:6, function(x) fdens(x, sim3,a = .05, ylim=c(0, xmax[x])))
reset()
#dev.off()

#get posterior fit of lambda. acma, arme and hear are good but others suck.
ex <- extract.samples(m3)
str(ex)
#simulate posterior
fpost <- function(spID){
  ll <- log(1) + ex$a[,spID] + rnorm(nrow(ex$a), 0, ex$sigmab[,spID])
  exp(ll)
}
xmax2 <- c(10000, 1000, 1000, 100, 1000, 1000)
par(mfrow=c(2,3))
post <- sapply(1:6, function(x)fpost(x))
#plot density of post
lapply(1:6, function(x) dens(post[,x], xlim=c(0,xmax2[x])))
#summary table of post
postm <- apply(post, 2, mean)
postpi <- apply(post, 2, PI)
dat %>% as.data.frame() %>% mutate(countest=count/prop) %>%  group_by(spID) %>% summarise(mean(countest))

postdf <- data.frame(spnames, spID=1:6, postm, t(postpi))
names(postdf)[c(4,5)] <- c('lower', 'upper')
ggplot(postdf, aes(spnames, postm))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper, width=.1))
