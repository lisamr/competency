#sporangia analysis.

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
#assume leaves are normally distributed around the individual mean
#assume individuals are normally distributed around the species mean
#generate a massive dataset with 1000 simulated ind and leaves
#sample from that dataset with schemes of interest. run model and evaluate output

#figure out what the distribution of the count data look like
ggplot(wideT, aes(countS, group=interaction(species, ind), color=ind))+
  geom_density()+
  facet_wrap(~species)
ggplot(wideT, aes(countS))+
  geom_density()
#maybe gamma?
fitdistr(wideT$countS[wideT$countS>0], "Gamma")
par.ind <- wideT %>% mutate(countS=countS+1) %>% group_by(species, ind) %>% summarise(shape=fitdistr(countS, "Gamma")$estimate[1], rate=fitdistr(countS, "Gamma")$estimate[2])

#for now, just assume individuals and leaves come from normal distribtuion

#get parameters from the data
#assume leaves are sampled from a normal distribution centered around each individual's mean
par.ind <- wideT %>% group_by(species, ind) %>% 
  summarise(mu=mean(countS), s=sd(countS), cv=s/mu)
#assume individuals are sampled from a normal distribution centered around each species' mean
par.sp <- wideT %>% group_by(species) %>% 
  summarise(mu=mean(countS), s=sd(countS), cv=s/mu)
#assume cv of variation is normally distributed. 
par.cv <- par.ind %>% group_by(species) %>% summarise(mu=mean(cv), s=sd(cv))

#generate additional data
set.seed(1)
nlvs=100
nind=100
l <- list(NULL)
for(s in 1:nrow(par.sp)){ #for each species
  
  #calculate mean of each individual
  mu_j <- rnorm(nind, par.sp$mu[s], par.sp$s[s]) 
  mat <- matrix(NA, nrow=nlvs, ncol=length(mu_j))
  
  for(i in 1:length(mu_j)){ #for each individual
    #generate more leaves 
    cv_j <- rnorm(1, par.cv$mu[s], par.cv$s[s])#calculate sd from cv
    y_ji <- rnorm(nlvs, mu_j[i], sd = abs(cv_j*mu_j[i]))
    mat[,i] <- y_ji
  }
  #turn the matrix into a dataframe
  m <- melt(mat) 
  names(m) <-  c("leaf", "ind", "countS")
  m$species <- par.sp$species[s]
  head(m)
  
  #add to a list
  l[[s]] <- m
}

master <- bind_rows(l)
master2 <- filter(master, countS>=0) %>%  #take only positive counts
  mutate(countS=round(countS) )#round all the counts to integers. a case for poisson or nbinom
ggplot(master2, aes(species, countS))+
  geom_boxplot()
ggplot(master2, aes( countS))+
  geom_density()+
  facet_wrap(~species)

#sample from the larger master dataset
#the number of individuals and leaves you're interested in sampling
indsamp=3
lvssamp=6
inds <- sample(1:nind, indsamp, replace = F)
master3 <- master2 %>% 
  filter(ind %in% inds) %>% #filter just the selected ind
  group_by(species, ind) %>%sample_n(lvssamp) #sample individuals

test <- function(indsamp, lvssamp){
  inds <- sample(1:nind, indsamp, replace = F)
  master3 <- master2 %>% 
    filter(ind %in% inds) %>% #filter just the selected ind
    group_by(species, ind) %>%sample_n(lvssamp) #sample individuals
  msim <- update(m2, data=master3)
  fixef(msim)
}

replicate(10, sample(1:nind, indsamp, replace = F))
msim2 <- update(m2, data=master3)
msim3 <- update(m2, data=master3)
multiplot(m2, msim, msim2, msim3)
multiplot(m2)


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


