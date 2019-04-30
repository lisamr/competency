#sporangia analysis.
library(dplyr)
library(ggplot2)
library(reshape2)
library(lme4)
library(coefplot)
library(emmeans)

wide <- read.csv("data/MASTERMERGED_wide.csv")
head(wide)

#Q1: does sporangia production vary by species? what is the rank order?

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
m2 <- glmer.nb(f0, data = wideT, control=glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))
#m2.1 <- update(m2, control=glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))#update so model converges
summary(m2)
multiplot( m1.1, m2)
plot(coefficients(summary(m2))[,3], 1:12)

#get confidence intervals. lots of warnings: singular fits, 3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :unable to evaluate scaled gradient 4: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
m2ci <- confint.merMod(m2, method = "boot")
m2cidf <- data.frame(coef=row.names(m2ci), lower= m2ci[,1], upper= m2ci[,2])
m2cidf <- cbind(m2cidf, rbind(NA, coefficients(summary(m2))[,1:2]))
row.names(m2cidf) <- 1:nrow(m2cidf)
names(m2cidf) <- c("coef", "lower", "upper", "est", "SE")
m2cidf <- m2cidf %>% mutate_at(.vars =vars("lower":"SE"), .funs = exp) %>% 
  mutate(species=c(NA, substr(coef, 8, 11)[-1]))
#write.csv(m2cidf, "~/Desktop/Competency project/output/sporangiapars.csv", row.names = F)

ggplot(m2cidf[-1,], aes(est, species))+
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.2)+
  #geom_errorbarh(aes(xmin=est-SE, xmax=est+SE), height=.1)+
  labs(x="# sporangia per leaf disc", y="species")

lsmeans(m2, "species") %>% pairs

hist(resid(m2)) #model overly underestimates counts (fat left tail)
qqnorm(resid(m2))+qqline(resid(m2))
plot(resid(m2)~predict(m2))#looks ok

#graphically, they all look very similar, but AIC favors the NB model
plot(density(wideT$countS), lwd=1.5)
lines(density(fitted(m1.1)), col="red", lwd=1.5)
lines(density(fitted(m2)), col="blue", lwd=1.5)
AIC(m1.1, m2)

lsmeans(m2, pairwise~species) #pairwise contrasts

#retry analysis removing ACMA3
wideT2 <- wideT %>% filter(!(species=="ACMA" & ind==3))
m3 <- update(m2, data=wideT2)
multiplot( m1.1, m2, m3)

#look at model fits
model=m2
plot(resid(model) ~ fitted(model))
plot(resid(model) ~ predict(model))
hist(resid(model))
qqnorm(resid(model))+qqline(resid(model))
plot(predict(model), log(model@resp$y)) +abline(0,1)

########################
########################
#simulate additional data to see how it would change results

#simulate data based on within- ind and species variation
#current model assumes constant variance across all individuals. variance increases with more counts.

#current number of treatment leaf discs
nrow(wideT) #12*3*6=216
#9 species (top + oaks) * 24. 8*3, 6*4, 4*6


library(simr)
simdata <- function(n.ind, n.leaf){
  #extend model to have new levels
  sim1 <- extend(m2, along = "ind", n=n.ind) 
  sim1 <- extend(sim1, along = "leaf", n=n.leaf) 
  #new data
  d <- cbind(getData(sim1), y=doSim(sim1)) 
  #fit the new data
  m <- doFit(d$y, sim1)
  return(m)
}
sim0 <- simdata(n.ind = 3, n.leaf = 6)
sim1 <- simdata(n.ind = 8, n.leaf = 3)
sim2 <- simdata(n.ind = 6, n.leaf = 4)
sim3 <- simdata(n.ind = 4, n.leaf = 6)
sim4 <- simdata(n.ind = 10, n.leaf = 6)
sim5 <- simdata(n.ind = 20, n.leaf = 6)

#need to do the simulation many times.
library(arm)
cbind(fixef(sim0), se.fixef(sim0))

multiplot(m2, sim0)
multiplot(m2, sim1, sim2, sim3, sim5)

#average SE of species
fse <-function(model) summary(model)$coefficients[,2] %>% mean
sapply(c(m2, sim1, sim2, sim3, sim4, sim5), fse) 

model=sim1
ggplot(getData(model), aes(species, countS, group=interaction(species, ind))) + geom_boxplot()
ggplot(getData(model), aes(species, countS, group=species)) +
  geom_boxplot()
coefplot(model)
lsmeans(model, "species") %>% pairs
plot(resid(model)~ predict(model))+abline(0,0)
qqnorm(resid(model))+qqline(resid(model))
plot(log(getData(model)$countS) ~ predict(model)) + abline(0,1)


#do it by hand
#simulate new data:
#y = species_effect + ind_effect + resids
#assume individuals are sampled from a normal distribution centered around each species' mean
#assume leaves are sampled from a normal distribution centered around each individual's mean

#first get parameters from the data
#species effect, individual effect, residuals
#resids = obs_i - sp_i - ind_i
mu.sp <- wideT %>% group_by(species) %>% 
  summarise(mu.sp=mean(countS))
mu.ind <- wideT %>% group_by(species, ind) %>% 
  summarise(mu.ind=mean(countS))
tmp <- wideT %>% select(species, ind, countS) %>% left_join(mu.sp) %>% left_join(mu.ind) %>% 
  mutate(ind.eff = mu.ind - mu.sp,
         resid = countS - mu.sp - ind.eff)
ggplot(tmp, aes(species, resid, group=interaction(species, ind)))+
  geom_boxplot() #since the residuals vary by species, adding residuals drawn from the same distribution doesn't seem right.

#resid variance increase with the response variable
ggplot(tmp, aes(countS, resid))+
  geom_point()


#establish simulated levels for within species and individual replication

sim.fun <- function(ninds, nleaves){
  ind.eff <- tmp %>% group_by(species) %>% summarise(m.ind = mean(ind.eff), sd.ind= sd(ind.eff)) #need sd of individual effect
  resids.eff <-  tmp %>% group_by(species) %>%  summarise(sd.resid= sd(resid)) #need sd of residuals, grouped by species
  
  #for each individual, simulate new individual effects from the mean and sd of the data's ind effects
  sim.ind <- function(x) rnorm(ninds, 0, ind.eff$sd.ind[x])
  sim.ind.eff <- sapply(1:nrow(ind.eff), sim.ind) 
  dfind <- data.frame(species=rep(unique(tmp$species), each=ninds ) , ind=1:ninds, ind.eff = as.vector(sim.ind.eff))
  #newdesign %>% left_join(mu.sp) %>% left_join(dfind) %>% head
  
  #for each leaf, simulate residuals based off the sd of the residuals, grouped by species and add to the dataframe
  sim.resid <- function(x) rnorm(ninds*nleaves, 0, resids.eff$sd.resid[x])
  sapply(1:nrow(ind.eff), sim.resid) 
  dfresid <- data.frame(species=rep(unique(tmp$species), each=ninds*nleaves ) , ind=1:ninds, leaf=rep(1:nleaves, each=ninds),  resid = as.vector(sapply(1:nrow(ind.eff), sim.resid) ))
  
  #merge simulated data with the simulated design
  newdesign <- expand.grid(species=unique(tmp$species), ind=1:ninds, leaf=1:nleaves)
  newdesign2 <- newdesign %>% left_join(mu.sp, by = "species") %>% 
    left_join(dfind, by = c("species", "ind")) %>% left_join(dfresid, by = c("species", "ind", "leaf")) %>% mutate(countS.est = mu.sp + ind.eff + resid)
  #newdesign2 %>% group_by(species, ind) %>% summarise(n())
  
  #hacky way of dealing with negative numbers: turn them into zeros? and round counts.
  newdesign3 <- newdesign2 %>% mutate(countS.est= ifelse(countS.est >= 0, countS.est, 0)) %>% mutate(countS.est = round(countS.est))
  
  return(newdesign3)
}

dat <- sim.fun(16, 2)
dat %>% arrange(species) %>% head
ggplot(dat, aes(species, countS.est, group=interaction(species, ind))) +geom_boxplot()
m2sim <- glmer.nb(countS.est ~ -1 + species + (1 | species:ind), data=dat, control=glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))
m2sim1 <- glmer(countS.est ~ -1 + species + (1 | species:ind), data=dat, family = poisson)
summary(m2sim1)
deviance(m2sim1)/df.residual(m2sim1) #not equal to 1, so overdispersed.


dat <- sim.fun(32, 1)
ggplot(dat, aes(species, countS.est)) +geom_boxplot()
ggplot(dat, aes(species, countS.est, group=interaction(species, ind))) +geom_point(alpha=.5)

m2sim <- MASS::glm.nb(countS.est ~ -1 + species, data=dat)
summary(m2sim)
m2sim2 <- glm(countS.est ~ -1 + species, data=dat, family=quasipoisson())
summary(m2sim2)
plot(m2sim2)
coefplot( m2sim2)

#after 1000 simulations, get mean and sd of countS for each species
#potential combinations within a tractable range: about 32 samples per species
#1:32, 2:16, 3:11, 4:8, 5:6 
L <- list(NULL)
for(i in 1:100){
  dat <- sim.fun(3, 50)
  datmean <- dat %>% group_by(species) %>% summarise(mu = mean(countS.est))
  L[[i]] <- datmean$mu
}
#rows are unique(dat$species) and in that order
test <- bind_cols(L) 
test2 <- data.frame(species = unique(dat$species), mean=apply(test, 1, mean), sd=apply(test, 1, sd))

ggplot(test2, aes(species, mean)) +
  geom_point()+
  geom_linerange(aes(ymin=mean-sd, ymax=mean+sd))+
  ylim(0,2000)

#real data
wideT %>% group_by(species) %>% summarise(mean=mean(countS), sd(countS), se=sd(countS)/sqrt(6)) %>% 
  ggplot(., aes(species, mean)) +
  geom_point()+
  geom_linerange(aes(ymin=mean-se, ymax=mean+se))+
  ylim(0,2000)

ggplot(wideT, aes(species, countS, group=interaction(species, ind)))+
  geom_boxplot()




###
#quick arithmetic to see how many samples we processed and what my options are.
wide %>% filter(trt=="T") %>% nrow #246 leaf samples (x 2 for controls too)
#I'll be inoculating 8 species, maybe 9
246/8 #30.75 samples per species. can round up to 32. 
#8 ind, 4 leaves; 7 ind, 4 or 5 leaves; 6 ind, 5 leaves; 5 ind, 6 lvs;

#what do the results look like with ind as fixed effect? Wacky estimates!
m3 <- glm.nb(countS ~ -1 + species+ species:ind, data = wideT)
summary(m3)
m3ci <- confint(m3, method="boot")

m3cidf <- data.frame(coef=row.names(m3ci), lower= m3ci[,1], upper= m3ci[,2])
m3cidf <- cbind(m3cidf, coefficients(summary(m3))[,1:2])
row.names(m3cidf) <- 1:nrow(m3cidf)
names(m3cidf) <- c("coef", "lower", "upper", "est", "SE")
m3cidf <- m3cidf %>% mutate_at(.vars =vars("lower":"SE"), .funs = exp) %>% 
  mutate(species=substr(coef, 8, 11))
#write.csv(m2cidf, "~/Desktop/Competency project/output/sporangiapars.csv", row.names = F)

#UHHHH, wow, those results are wrong!
ggplot(m3cidf[1:12,], aes(est, species))+
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=.2)+
  #geom_errorbarh(aes(xmin=est-SE, xmax=est+SE), height=.1)+
  labs(x="# sporangia per leaf disc", y="species")

####################################################
####################################################
#I'm thinking of removing some species for the redo experiment. ARME, HEAR, and VAOV seem like they didn't sporulate much. 
wideT %>% group_by(species) %>% summarise(
  muS = mean(countS), maxS = max(countS), 
  muC = mean(countC), maxC = max(countC))
ggplot(wideT, aes(species, countS, group=interaction(species, ind))) +
  geom_boxplot() 
ggplot(wideT, aes(species, countS+1, group=interaction(species, ind))) +
  geom_boxplot() +
  scale_y_continuous(trans='log10')
ggplot(wideT, aes(species, countC, group=interaction(species, ind))) +
  geom_boxplot() 
ggplot(wideT, aes(species, countC+1, group=interaction(species, ind))) +
  geom_boxplot() +
  scale_y_continuous(trans='log10')
