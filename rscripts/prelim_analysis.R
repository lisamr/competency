#first stab at analyzing the (preliminary) data. need to run a power analysis.

setwd("~/Desktop/Competency project/")

library(dplyr)
library(ggplot2)
library(reshape2)
library(lme4)
library(coefplot)
library(emmeans)

wide <- read.csv("data/MASTERMERGED_wide.csv")
head(wide)

#Q1: does sporangia production vary by species? what is the rank order?
#Q2: does chlamydospore production vary by species? what is the rank order?
#Q3: is there a relationship between spore production (infectivity) and lesion area (susceptibility)? (nix this question. i don't think leaf lesion area is predictive of suspectilbity. looking ath the avg area vs species shows that typically unaffected hosts have lots of lesions--acma, umca, arme--while quag and quch are hardly affected. susceptility is probably better measured with lesions in the wood, not leaf.)

#QUESTION1
#sporangia count distribution is overdispersed. maybe zero-inflated poisson or negative binomial? 
#RV: countS, PV: species + ind OR species + (1|ind). with 3 individuals wouldnt make sense to make that a random effect, but with more it would.

#clean up data. remove CEOL (only 1 ind)
wideT <- filter(wide, trt=="T", !is.na(countS), species!="CEOL") %>% mutate(ind2=interaction(species, ind)) %>% droplevels()
wideT %>% group_by(species, ind) %>% summarise(min(countS), max(countS)) 

ggplot(wideT, aes(countS))+
  geom_histogram()
ggplot(wideT, aes(species, countS, group=interaction(species, ind))) +
  geom_boxplot() 
ggplot(wideT, aes(species, countS)) +
  geom_boxplot() 

#lets look at the relationship between variance and mean of the groups. that will inform how overdispersed the data is. 

mu.v <- wideT %>% group_by(species, ind) %>% summarise(mu = mean(countS), v = var(countS))
ggplot(mu.v, aes(mu, v))+
  geom_point()+
  geom_smooth(method="lm",formula=y~x-1)+##(quasi-Poisson/NB1) fit
  geom_smooth(colour="red")+## smooth (loess)
  geom_smooth(method="lm",formula=y~I(x^2)-1,colour="purple")+## semi-quadratic (NB2/LNP)
  geom_abline(a=0,b=1,lty=2)## Poisson (v=m)
#Poisson is obvious ill-suited. Either neg. binomial might work. leverage is high on that last point, which I bet is ACMA3. filter that out and redo.
mu.v2 <- wideT %>% filter(!(species=="ACMA"&ind==3)) %>%  group_by(species, ind) %>% summarise(mu = mean(countS), v = var(countS))
ggplot(mu.v2, aes(mu, v))+
  geom_point()+
  geom_smooth(method="lm",formula=y~x-1)+##(quasi-Poisson/NB1) fit
  geom_smooth(colour="red")+## smooth (loess)
  geom_smooth(method="lm",formula=y~I(x^2)-1,colour="purple")+## semi-quadratic (NB2/LNP)
  geom_abline(a=0,b=1,lty=2)## Poisson (v=m)

#start with just poisson as reference.
f0<-formula(countS ~ -1 + species+(1|species:ind)) 
m1.1 <- glmer(f0, data = wideT, family="poisson")
summary(m1.1)
#these two should also be the same. neg. binom.
m2 <- glmer.nb(f0, data = wideT)
m2.1 <- update(m2, control=glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))#update so model converges

#zeroinflated neg.binom. #wont run with random effects
m2.2 <- glmmADMB::glmmadmb(f0, data = wideT, zeroInflation=T, family="nbinom1")
coefplot(m2.1) #well that's uninformative.
m2.3 <- glmmadmb(f0, data = wideT, zeroInflation=T, family="nbinom1")#wont run.
multiplot( m1.1, m2, m2.1)

#graphically, they all look very similar, but AIC favors the NB model
plot(density(wideT$countS), lwd=1.5)
lines(density(fitted(m1.1)), col="red", lwd=1.5)
lines(density(fitted(m2.1)), col="blue", lwd=1.5)
AIC(m1.1, m2.1)

lsmeans(m2, pairwise~species) #pairwise contrasts

#how do the models look when I treat individuals as a fixed effect? Its identical to blue, the poisson.
f1<-formula(countS ~ species+ind2) 
m3.1 <- glm(f1, data = wideT, family = poisson)
m3.11 <- glm(f1, data = wideT, family = quasipoisson)
m3.2 <- glm.nb(f1, data = wideT)
plot(density(wideT$countS), lwd=1.5)
lines(density(fitted(m3.11)), col="red")
lines(density(fitted(m3.2)), col="pink")
multiplot(m1.1, m3.1) #poisson random vs fixed. RE more conservative.
multiplot(m2.1, m3.2) #NB random vs fixed. RE more conservative.
AIC(m3.1, m3.11, m3.2) #i dont think aic can compare these models properly
multiplot(m3.1, m3.11, m3.2)
coefplot(m3.11)
lsmeans(m3.11, pairwise~species) 
summary(m3.11)

#try a hurdle/zinfl model
test <- pscl::hurdle(countS ~ species+ind2|1, data = wideT, dist = "negbin")
############################################################
############################################################

#look at chlamydospores?
ggplot(wide, aes(countC))+
  geom_histogram()+
  facet_wrap(~trt)
ggplot(filter(wide, trt=="T"), aes(species, countC, group=interaction(species, ind))) +
  geom_boxplot() 
ggplot(filter(wide, trt=="T"), aes(species, countC)) +
  geom_boxplot() 

#lets look at the relationship between variance and mean of the groups. that will inform how overdispersed the data is. quasi-poisson/NB1 looks pretty good.
wideTC <- filter(wide, trt=="T", !is.na(countC)) %>% mutate(ind2=interaction(species, ind), countC=round(countC))
mu.v <- wideT %>% group_by(species, ind) %>% summarise(mu = mean(countC), v = var(countC))
ggplot(mu.v, aes(mu, v))+
  geom_point()+
  geom_smooth(method="lm",formula=y~x-1)+##(quasi-Poisson/NB1) fit
  geom_smooth(colour="red")+## smooth (loess)
  geom_smooth(method="lm",formula=y~I(x^2)-1,colour="purple")+## semi-quadratic (NB2/LNP)
  geom_abline(a=0,b=1,lty=2)## Poisson (v=m)

#models. ME models first. ZINB models won't run :(
mR.NBa <- glmer.nb(countC ~ -1+ species + (1|ind2), data = wideTC)#didnt actually converge
mR.NB <- glmmadmb(countC ~ -1+ species + (1|ind2), data = wideTC, zeroInflation=F, family="nbinom")
mR.NB1 <- glmmadmb(countC ~ -1+ species + (1|ind2), data = wideTC, zeroInflation=F, family="nbinom1")
mR.ZINB <- glmmadmb(countC ~ -1+ species + (1|ind2), data = wideTC, zeroInflation=T, family="nbinom")
mR.ZINB1 <- glmmadmb(countC ~ -1+ species + (1|ind2), data = wideTC, zeroInflation=T, family="nbinom1")
plot(density(wideTC$countC), lwd=1.5)
lines(density(fitted(mR.NBa)), col="green")
lines(density(fitted(mR.NB)), col="blue")
lines(density(fitted(mR.NB1)), col="purple")
summary(mR.NB)

AIC(mR.NBa, mR.NB, mR.NB1)
anova(mR.NB, mR.NB1) #NB1 is supposedly better
multiplot(mR.NB, mR.NB1)
coefplot(mR.NB1) #why are the error bars on CEOL and VAOV so large when those have the smallest sd (0)!!!
wideTC %>% group_by(species) %>% summarise(nind = length(unique(ind2)), samples = n(), mean(countC), sd(countC))

#ind as fixed effect. both models are identical.
f1<-formula(countC ~ -1+species+ind2) 
mF.QP <- glm(f1, data = wideTC, family = quasipoisson)
mF.NB <- glm.nb(f1, data = wideTC)
mm <- glmmadmb(f1, data = wideTC,family="nbinom")
multiplot(mF.QP, mF.NB)#i feel like these are wrong!
plot(density(wideTC$countC), lwd=1.5)
lines(density(fitted(mR.NB)), col="blue")
lines(density(fitted(mR.NB1)), col="purple")
lines(density(fitted(mF.QP)), col="red")
lines(density(fitted(mF.NB)), col="pink")
lsmeans(mF.QP, pairwise~species)
summary(mF.QP) #dispersion factor 94!
summary(mF.NB) #dispersion factor 94!

#try hurdle and zibn with pscl. I don't know how.
test <- zeroinfl(f1, data = wideTC, dist = "negbin")#fails because some coef. are not identifiable. omit them.
## set up version of data with non-identified regressors omitted
nas <- names(which(is.na(coef(mF.NB))))
nas <- substr(nas, 5, 10)
wideTC2 <-wideTC %>% filter(!ind2 %in% nas)
## re-fit glm.nb(). still doesn't work
mF.NB2 <- glm.nb(f1, data = wideTC2)
## fit zeroinfl(), still doesn't work
fm2a <- zeroinfl(f1, data = wideTC3, dist = "negbin")


plot(fitted(mR.NBa) ~ wideTC$countC,  col="green")
points(fitted(mR.NB) ~ wideTC$countC,  col="blue")
points(fitted(mR.NB1) ~ wideTC$countC,  col="purple")
abline(0,1)
abline(lm(fitted(mR.NBa) ~ wideTC$countC), col="green")
abline(lm(fitted(mR.NB) ~ wideTC$countC), col="blue")
abline(lm(fitted(mR.NB1) ~ wideTC$countC), col="purple")
points(fitted(mF.QP) ~ wideTC$countC,  col="red")
points(fitted(mF.NB) ~ wideTC$countC,  col="pink")
abline(0,1)
abline(lm(fitted(mF.QP) ~ wideTC$countC), col="red")
abline(lm(fitted(mF.NB) ~ wideTC$countC), col="pink")

#look at the variance within the different groups: 1. individuals nested within 2. species
ggplot(wideTC, aes(species, countC, group=interaction(species, ind))) +
  geom_boxplot() 
sumi <- wideTC %>% group_by(species, ind2) %>% summarise(m=mean(countC), v=var(countC), cv=sd(countC)/mean(countC))
sums <- wideTC %>% group_by(species) %>% summarise(m=mean(countC), v=var(countC), cv=sd(countC)/mean(countC))
ggplot(sumi, aes(species, cv, group=ind2))+
  geom_point(position=position_dodge(.8))
ggplot(sums, aes(species, v))+
  geom_point(position=position_dodge(.8))
##################################################
##################################################
#group individual level data and run glm
groupeddata <- wideT %>% group_by(species, ind2) %>% summarise(S = mean(countS, na.rm=T), C = mean(countC, na.rm=T))
mgrouped <- glm(S ~ -1+ species, data=groupeddata, family = quasipoisson)
summary(mgrouped)
coefplot(mgrouped)
rpois(10, 10)

plot(density(groupeddata$S), lwd=1.5)
lines(density(fitted(mgrouped)), col="blue")

mgroupedC <- glm(C ~ -1+ species, data=groupeddata, family = quasipoisson)
summary(mgroupedC)
coefplot(mgroupedC)
plot(density(na.omit(groupeddata$C) ), lwd=1.5)
lines(density(fitted(mgroupedC)), col="blue")
plot(fitted(mgroupedC) ~ na.omit(groupeddata$C))+abline(0,1)

#in the end, both sporangia and chlamydospore counts cannot be compared very well because theres so much intraspecific variation. treating ind as fixed effects is just as good (performs better because less conservative) but i would probably want to treat them as random effects with more of them.
#the distributions are overdispersed count data, which means negative binomial or quassipoisson, or ZI-NB, which I havent managed to do. 
#I need to do a power analysis. hopefully can fudge around with the data and choose how variable the data is. I think the plants will be less variable when they're all from Big Sur.
#power analysis: use a quassipoisson or NB distribution with individuals as fixed effects (I know that can be debated) for both sporangia and chlamydospore counts.
