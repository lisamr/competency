##########################################
# RELATIONSHIP BETWEEN SPORANGIA AND LESION SIZE?

setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)
library(scales)
library(ggridges)

#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf treatment assay
df1 <- df %>% filter(leafID<=530, trt=="T", spore_assay=="S") 

#visualize first
#make data tall in regards to counts
dtall2 <- melt(df1, id.vars = c("species", "leafID", "perc_lesion"), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count)) %>% 
  mutate(perc_lesion=perc_lesion/100) %>% 
  arrange(leafID) 
head(dtall2)

#all counts vs lesions
ggplot(dtall2, aes(perc_lesion, count))+
  geom_point()+
  facet_wrap(~species, scales = "free_y")

#just mean counts vs lesions
df2 <- dtall2 %>% group_by(species, leafID) %>% summarise(lesion=mean(perc_lesion), countm=mean((count))) #log counts or not?
ggplot(df2, aes(lesion, countm))+
  geom_point()+
  facet_wrap(~species, scales = "free_y")
  #geom_smooth(method = 'lm') #probably need to remove ACMA's outlier because it has too much leverage


########################################
#try it with glmer first
m0 <- glmer(count ~ perc_lesion + (1|leafID) + (perc_lesion|species), data = dtall2, family = 'poisson')
summary(m0)

#thats cool it worked. Not sure how you check out model performance.
plot(m0)
plot(dtall2$count, predict(m0, type='response'))+abline(0,1)

#how do you do counterfactual plots again with lmer obj??

########################################
#statistical models!
#figure out priors. modeling counts ~ species + lesion 

#prior for intercept. mass should be under 20 or so, which is the maxish spores produced regardless of lesion area.
a <- rnorm(5000, 3, 1.5)
a2 <- exp(a)
dens(a2, xlim=c(0,100), xlab="sp intercept")
mean(a2)

#prior for slope. probably shouldnt go down (or much at least)
b <- rnorm(5000, 1, .5)
sp <- seq(0,1, length.out = 100)
lam <- sapply(1:length(sp), function(x) exp(a + b*sp[x]))
plot(NULL, xlim=c(0,1), ylim=c(0,200), xlab="lesion", ylab="#spores")
lapply(1:200, function(x) lines(sp, lam[x,], col=alpha(1, .3)))

#prior for rho, the correlation between int and slopes. should be centered around zero.
R <- rlkjcorr(5000, 2, 2) #array
dens(R[,1,2], xlab='correlation')
###################################
###################################
#model it?!? having issues. maybe make simpler by removing 3 subsamples (leafID) to troubleshoot and then add it back in?
dtall3 <- dtall2 %>% filter(!duplicated(leafID))
dat <- list(
  count=dtall3$count,
  species=as.integer(factor(dtall3$species)),
  #ID=dtall3$leafID,
  lesion=dtall3$perc_lesion)

#no leafID
m1 <- ulam(
  alist(
    #likelihood and linear model
    count ~ dpois(lambda),
    log(lambda) <- a[species] + b[species]*lesion,
    #varying effects
    c(a, b)[species] ~ multi_normal(c(abar, bbar), rho, sigma),
    #fixed priors
    abar ~ dnorm(3, 1.5),
    bbar ~ dnorm(1, .5),
    rho ~ dlkjcorr(2),
    sigma ~ dexp(1)
    ), 
  data=dat , chains=4 , cores=4 )
traceplot(m1)
precis(m1, depth = 3) #sampled well!

########################
#put leafID back in?
dat2 <- list(
  count=dtall2$count,
  species=as.integer(factor(dtall2$species)),
  ID=as.integer(factor(dtall2$leafID)),
  lesion=dtall2$perc_lesion)
#no leafID
m2 <- ulam(
  alist(
    #likelihood and linear model
    count ~ dpois(lambda),
    log(lambda) <- a[species] + g[ID] + b[species]*lesion,
    #varying effects
    g[ID] ~ dnorm(gbar, sigmag),#not multivariate cuz shouldnt vary with slope?
    c(a, b)[species] ~ multi_normal(c(abar, bbar), rho, sigma),
    #fixed priors
    gbar ~ dnorm(0, 2.2),
    sigmag ~ dexp(1),
    abar ~ dnorm(3, 1.5),
    bbar ~ dnorm(1, .5),
    rho ~ dlkjcorr(2),
    sigma ~ dexp(1)
  ), 
  data=dat2 , chains=4 , cores=4 )

#uhhh, did not sample well at all. model must be wrong.
traceplot(m2) 
precis(m2, pars = 'a', 'b') 
