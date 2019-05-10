#power analysis for telling 2 species apart. 
#assume 2 species are slightly different. how many samples do you need to tell them apart? Use the existing data to estimate alpha and sigma.

library(dplyr)
library(MASS)
library(lme4)
library(ggplot2)

#for model evaluation
eval <- function(model){
  par(mfrow=c(2,2))
  plot(resid(model) ~ fitted(model))
  hist(resid(model))
  qqnorm(resid(model))+qqline(resid(model))
  par(mfrow=c(1,1))
}

#load data
wide <- read.csv("data/MASTERMERGED_wide.csv")
wide$ind2 <- interaction(wide$species, wide$ind)
#clean up data. remove CEOL (only 1 ind)
wideT <- filter(wide, trt=="T", !is.na(countS), species!="CEOL") %>% mutate(ind2=interaction(species, ind)) %>% droplevels()
head(wideT)

#don't do this. uses pseudoreplicated data when I'm not simulating how within individual variation affects results
wideT %>% group_by(species) %>% summarise(mu=mean(countS), sigma=sd(countS))

#first collapse within individual variation and then calculate mean and sd for each species off of 3 data points. 
sumtab <- wideT %>% group_by(species, ind) %>% summarise(mu=mean(countS)) %>% group_by(species) %>% summarise(mean=mean(mu), sd=sd(mu))

#simulate additional data to see how many replicates you need to tell speciesA and speciesB apart.
ggplot(wideT %>% filter(!species=="ACMA"), aes(species, countS))+
  geom_boxplot()

#todi and umca seem like good species to compare
N=20
A="UMCA"
B="RHCA"
spA <- rnorm(N, sumtab$mean[sumtab$species==A], sumtab$sd[sumtab$species==A]) %>% round %>% ifelse(. < 0, 0, .)
spB <- rnorm(N, sumtab$mean[sumtab$species==B], sumtab$sd[sumtab$species==B]) %>% round %>% ifelse(. < 0, 0, .)
dat <- data.frame(species= c(rep(A, length(spA)), rep(B, length(spB))), spores=c(spA, spB))
dat$obs <- 1:nrow(dat)
ggplot(dat, aes(species, spores))+
  geom_boxplot()+
  geom_point()
m0 <- lm(spores ~ species, data=dat)
m1 <- glm(spores ~ species, data=dat, family = poisson)#too overdispersed (6352/18 >>>> 1)
m2 <- glm(spores ~ species, data=dat, family = quasipoisson)
m3 <- glm.nb(spores ~ species, data=dat)
summary(m2)
coefplot::coefplot(m2)
eval(m0)
eval(m1) #poisson
eval(m2) #quasipoisson: same as poisson
eval(m3) #nbinom: same as poisson
#seems like m2 is a safe bet to get a pvalue from.
coef(summary(m2))[2,4] #p-value for species effect

#now loop this. can change species easily.
N=seq(1, 50, by=5)
A="LIDE"
B="QUCH"
nsims <- 500
pvec <- rep(NA, nsims)

m <- matrix(NA, ncol = nsims, nrow=length(N))#matrix to be filled with simulated pvalues

for (j in 1:length(N)) {#for each # of replicates
  n=N[j]

  for(i in 1:nsims){ #simulate a run nsims times
    #simulate data
    spA <- rnorm(n, sumtab$mean[sumtab$species==A], sumtab$sd[sumtab$species==A]) %>% round %>% ifelse(. < 0, 0, .)
    spB <- rnorm(n, sumtab$mean[sumtab$species==B], sumtab$sd[sumtab$species==B]) %>% round %>% ifelse(. < 0, 0, .)
    dat <- data.frame(species= c(rep(A, length(spA)), rep(B, length(spB))), spores=c(spA, spB))
    
    #run model
    m2 <- glm(spores ~ species, data=dat, family = quasipoisson)
    
    #put p-value into a vector
    pvec[i] <- coef(summary(m2))[2,4] #p-value for species effect
  }
  #simulated pvalues. rows=reps, cols=different simulation
  m[j,] <- pvec
}

#power=proportion of times you get a pvalue < .05
pow <- apply(m, 1, function(x)sum(x<=.05)/length(x))
df <- data.frame(N, pow)
par(mfrow=c(1,1))
plot(df, type='b', ylim=c(0,1))+abline(.80, 0, lty=2)

