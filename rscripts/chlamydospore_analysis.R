#chlamydospore analysis.

setwd("~/Desktop/Competency project/")

library(dplyr)
library(ggplot2)
library(reshape2)
library(lme4)
library(coefplot)
library(emmeans)

wide <- read.csv("data/MASTERMERGED_wide.csv")
head(wide)

#Q1: does chlamydospore production vary by species? 

#QUESTION1
#chlamydospore count distribution is overdispersed. maybe zero-inflated poisson or negative binomial? 
#RV: countS, PV: species + ind OR species + (1|ind). with 3 individuals wouldnt make sense to make that a random effect, but with more it would.

#clean up data. remove CEOL (only 1 ind)
wideT <- filter(wide, trt=="T", !is.na(countC), !species=="CEOL") %>% mutate(ind2=interaction(species, ind), countC=round(countC)) %>% 
  filter(!(ind2=="UMCA.1"& leaf==4))%>% #drop fell off 
  droplevels()
wideT %>% group_by(species) %>% summarise(length(ind))
sumdat <- wideT %>% group_by(species, ind) %>% summarise(min=min(countC), max=max(countC), mu=mean(countC), sd=sd(countC), v=var(countC)) 

ggplot(wideT, aes(countC))+
  geom_histogram()
ggplot(wideT, aes(species, countC, fill=species, group=interaction(species, ind))) +
  geom_boxplot() 
ggplot(wideT, aes(species, countC)) +
  geom_boxplot() 

#figure out what kind of model. NB1 looks good.
ggplot(sumdat, aes(mu, v))+
  geom_point()+
  geom_smooth(method="lm",formula=y~x-1)+##(quasi-Poisson/NB1) fit
  geom_smooth(colour="red")+## smooth (loess)
  geom_smooth(method="lm",formula=y~I(x^2)-1,colour="purple")+## semi-quadratic (NB2/LNP)
  geom_abline(a=0,b=1,lty=2)## Poisson (v=m)

#start with just poisson as reference.
f0<-formula(countC ~ -1 + species+(1|species:ind)) 
m1 <- glmer(f0, data = wideT, family="poisson") #can't estimate VAOV
summary(m1)

#these two should also be the same. neg. binom.
m2 <- glmer.nb(f0, data = wideT)#26: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  ... :Hessian is numerically singular: parameters are not uniquely determined
#m2.1 <- update(m2, control=glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))#update so model converges
summary(m2)
multiplot( m1, m2)

#the above models having problems converging. what does removing VAOV, RHCA and HEAR do? Just model the ones that make chlamydospores.
wideTfilt <- wideT %>% filter(!species %in% c("VAOV"))
wideTfilt2 <- wideT %>% filter(!species %in% c("VAOV", "RHCA", "HEAR"))
m1.1 <- update(m1, data=wideTfilt)
summary(m1.1)
m2.1 <- update(m2, data=wideTfilt)
summary(m2.1)
m2.2 <- update(m2, data=wideTfilt2)
summary(m2.2)

multiplot( m1.1, m2.1, m2.2)
AIC(m1.1, m2.1)

#model eval
model=m2.2
plot(resid(model)~fitted(model))
hist(resid(model)) #model grossly underestimates some of the values
plot(log(model@resp$y) ~ predict(model))+abline(0,1)
plot(model@resp$y ~ exp(predict(model)))+abline(0,1)
plot(density( log(model@resp$y)), lwd=2)
lines(density( predict(model)))

#plot conf int
m2.1ci <- confint.merMod(m2.1, method = "boot")#singular fit and ... 49: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  ... :Model failed to converge with max|grad| = 1.16616 (tol = 0.001, component 1)50: In (function (fn, par, lower = rep.int(-Inf, n), upper = rep.int(Inf,  ... :failure to converge in 10000 evaluations
m2.2ci <- confint.merMod(m2.2, method = "boot")#model convergence failure
modelci=m2.2ci
mod=m2.2

modelcidf <- data.frame(coef=row.names(modelci), lower= modelci[,1], upper= modelci[,2])
modelcidf <- cbind(modelcidf, rbind(NA, coefficients(summary(mod))[,1:2]))
row.names(modelcidf) <- 1:nrow(modelcidf)
names(modelcidf) <- c("coef", "lower", "upper", "est", "SE")
modelcidf <- modelcidf %>% mutate_at(.vars =vars("lower":"SE"), .funs = exp) %>% 
  mutate(species=c(NA, substr(coef, 8, 11)[-1]))
#write.csv(modelcidf, "~/Desktop/Competency project/output/chlamydopars2.csv", row.names = F)

ggplot(modelcidf[-1,], aes(est, species))+
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.2)+
  #geom_errorbarh(aes(xmin=est-SE, xmax=est+SE), height=.1)+
  labs(x="# chlamydospores per leaf disc", y="species")


