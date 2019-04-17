#compare sporangia and chlamydospore production
library(dplyr)
library(ggplot2)
library(reshape2)
library(lme4)
library(coefplot)
library(emmeans)

wide <- read.csv("data/MASTERMERGED_wide.csv")
head(wide)

#Q1: what is the relationship between sporangia and chlamydospore production and does it vary by species?

#clean up data. remove CEOL (only 1 ind) 
wideT <- filter(wide, trt=="T", !is.na(countS), !is.na(countC), species!="CEOL") %>% mutate(ind2=interaction(species, ind), countC=round(countC)) %>% droplevels()
wideT %>% group_by(species, ind) %>% summarise(length(unique(leaf)))#all counts are paired, in that they came from the same leaf.
wideT %>% group_by(species, ind) %>% summarise(min(countS), max(countS)) 

ggplot(wideT, aes(countC, countS, color=as.factor(ind), group=ind))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  facet_wrap(~species, scales = "free", ncol=4)+
  scale_color_manual(values=viridis::viridis(3, begin = .1, end = .95))+
  labs(x="# chlamydospores", y="# sporangia", color="individual")

ggplot(wideT, aes(countC, countS, color=species, group=interaction(species, ind)))+
  geom_point()+
  geom_smooth(method='lm', se=F)
ggplot(wideT, aes(countC, countS, color=species, group=interaction(species)))+
  geom_point()+
  geom_smooth(method='lm', se=F)
ggplot(wideT, aes(log(countC+1), log(countS+1), color=species, group=interaction(species)))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  facet_wrap(~species)

#run model...
#fullest model wont converge.
m1 <- glmer.nb(countS ~ countC*species + (countC|species:ind), data = wideT)
#null
m2 <- glmer.nb(countS ~ 1 + (1|species:ind), data = wideT)
m2.1 <- glmer.nb(countS ~ 1 + (1|species) + (1|species:ind), data = wideT)
m2.2 <- glmer.nb(countS ~ countC + (countC|species) + (1|species:ind), data = wideT) #if youre not interested in seeing how species covary with C
m3 <- glmer.nb(countS ~ countC:species + (1|species:ind), data = wideT)
m4 <- glmer.nb(countS ~ countC*species + (1|species:ind), data = wideT) #doesn't run
m5 <- glmer.nb(countS ~ countC + (1|species) + (1|species:ind), data = wideT) #not converging 
m6 <- glmer.nb(countS ~ -1 + countC + species + (1|species:ind), data = wideT) #not converging, singular fit 
summary(m6)
AIC( m5, m3, m2.1, m2.2, m2)

#remove HEAR and RHCA since they don't make chlaydospores. maybe causing problems? Nope doesnt help. data too variable?
wideT2 <- wideT %>% filter(!species %in% c("HEAR", "RHCA"))
m1 <- glmer.nb(countS ~ countC + (countC|species/ind2), data = wideT2)#no converge
m1 <- glmer.nb(countS ~ countC + (1|species/ind2), data = wideT2)#no converge
m1 <- glmer.nb(countS ~ countC + species + (1|ind2), data = wideT2)#no converge

#treat ind as fixed effects
#fullest
m7 <- MASS::glm.nb(countS ~ countC*species*ind2, data = wideT)
m8 <- MASS::glm.nb(countS ~ countC*species + countC*ind2, data = wideT)
m8.1 <- MASS::glm.nb(countS ~ -1+countC*species + ind2, data = wideT)
m8.2 <- MASS::glm.nb(countS ~ species + countC*ind2, data = wideT)
m8.3 <- MASS::glm.nb(countS ~ -1+countC/species + ind2, data = wideT)
m9 <- MASS::glm.nb(countS ~ -1+ countC + species + ind2, data = wideT)
m10 <- MASS::glm.nb(countS ~ 1, data = wideT)
AIC(m10, m9, m8.3, m8.2, m8.1, m8, m7) #all models besides null are pretty much equal. will look at m9:simplest

summary(m9)
AIC(m4, m5) #no diff with interaction term. choose simpler.
anova(m4, m5)

#model evaluation
model=m9
hist(resid(model))
qqnorm(resid(model)) + qqline(resid(model)) 
plot(resid(model) ~ predict(model))
plot(log(model$y) ~ predict(model)) + abline(0,1)

wideT$pred <- exp(predict(model))
ggplot(wideT, aes(countC, pred, color=ind, group=interaction(ind, species)))+
  geom_point(alpha=.5)+
  geom_smooth(method='lm', se=F, aes(alpha=.2)) +
  geom_point(aes(countC, countS)) +
  geom_smooth(aes(countC, countS), method='lm', se=F) +
  facet_wrap(~species, scales = 'free')

#try a model with aggregated data?
agg <- wideT %>% group_by(species, ind) %>% summarise(mS = mean(countS), mC = mean(countC))
ggplot(agg, aes(mC, mS))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  facet_wrap(~species, scales="free")

m11 <- lm(mS ~ -1+mC*species, agg)
summary(m11)
