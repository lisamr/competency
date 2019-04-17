#relationship between sporangia production and ubiquity of species, controlling for replication (random effects) of species and individual
#ideal model: sporangia ~ ubiquity (at the species level) + 1|species + 1|individual
#sporangia ~ ubiquity + 1|species/ind ?

#sporangia analysis.
library(dplyr)
library(ggplot2)
library(reshape2)
library(lme4)
library(coefplot)
library(emmeans)

wide <- read.csv("data/MASTERMERGED_wide.csv")
head(wide)
wide$ind2 <- interaction(wide$species, wide$ind)

#clean up data. remove CEOL (only 1 ind)
wideT <- filter(wide, trt=="T", !is.na(countS), species!="CEOL") %>% mutate(ind2=interaction(species, ind)) %>% droplevels()

#first try with mixed evergreen forests
species <- read.csv("output/mixedeverspecies.csv")
head(species)
species <- species %>% mutate(rank = rank(-species$n.plots, ties.method = "min")) 
head(wideT)

#merge the two datasets
dfmerged <- wideT %>% left_join(species, by=c("species"="Species")) %>% 
  filter(!is.na(rank))
head(dfmerged)

ggplot(dfmerged, aes(n.plots, countS, group=interaction(species, ind))) +
  geom_boxplot()
ggplot(dfmerged, aes(rank, countS, group=interaction(species, ind))) +
  geom_boxplot()+
  geom_text(data=filter(dfmerged, !duplicated(species)), aes( label = species),  vjust=-10)

#testing relationship between species rank and spore production
#m1 and m2 are the same. was testing the syntax.
m0 <- glmer.nb(countS ~ 1 + (1|species) + (1|ind2), data=dfmerged) #null
m1 <- glmer.nb(countS ~ rank + (1|species/ind2), data=dfmerged)
m2 <- glmer.nb(countS ~ rank + (1|species) + (1|ind2), data=dfmerged)
summary(m0)
summary(m1)
summary(m2)
anova(m0, m2) #no difference between null and m2. ie. no relationship between rank and sporangia.

model=m0
plot(resid(model)~ predict(model))
plot(model@resp$y %>% log ~ predict(model)) + abline(0,1)
plot(density(model@resp$y))
lines(density(exp(predict(model))))

######
#repeat for redwood forests
species <- read.csv("output/redwoodspecies.csv")
head(species)
species <- species %>% mutate(rank = rank(-species$n.plots, ties.method = "min")) 
head(wideT)

#merge the two datasets
dfmerged <- wideT %>% left_join(species, by=c("species"="Species")) %>% 
  filter(!is.na(rank))
head(dfmerged)

ggplot(dfmerged, aes(n.plots, countS, group=interaction(species, ind))) +
  geom_boxplot()
ggplot(dfmerged, aes(rank, countS )) +
  geom_boxplot(aes(group=interaction(species, ind)))+
  geom_text(data=filter(dfmerged, !duplicated(species)), aes( label = species),  vjust=-10) 

#testing relationship between species rank and spore production
#m1 and m2 are the same. was testing the syntax.
m0 <- glmer.nb(countS ~ 1 + (1|species/ind2), data=dfmerged) #null
m1 <- glmer.nb(countS ~ rank + (1|species/ind2), data=dfmerged)
anova(m0,m1)
#no difference: rank does not predict countS
model=m0
plot(resid(model)~ predict(model))
plot(model@resp$y %>% log ~ predict(model)) + abline(0,1)
plot(density(model@resp$y), lwd=2)
lines(density(exp(predict(model))), col="blue")

############################################################
############################################################
#see if my proposed model tells me what I want to learn regarding how ubiquity affects sporangia counts. will simulate data.

library(dplyr)
library(ggplot2)
library(lme4)

#make data
spp <- LETTERS[1:10]
nplots <- (rlnorm(10, sdlog = .5)*20) %>% round %>% sort
spore.means <- runif(10, 0, 1000) %>% sort

N=32
spores <- sapply(1:length(spore.means), function(x) rnbinom(N, mu = spore.means[x], size = 10) ) %>% as.vector()

df <- data.frame(spp, nplots, spore.means)
df2 <- df[rep(seq_len(nrow(df)), each=N),]
df2$spores <- spores
row.names(df2) <- 1:nrow(df2)
head(df2)

#plot data
#ggplot(df2, aes(spp, spores))+ geom_boxplot()
ggplot(df2, aes(log(nplots), spores, group=spp))+
  geom_boxplot()

#model 
m2 <- glmer.nb(spores ~ (nplots) + (1|spp), data=df2)
summary(m2)
arm::display(m2)

m1 <- glmer(spores ~ nplots + (1|spp), family = 'poisson', data=df2)
summary(m1)

plot(m1)
eval <- function(mod){
  plot(mod)
  hist(resid(mod))
  qqnorm(resid(mod))
  qqline(resid(mod))
  plot(fitted(mod), df2$spores) + abline(0,1)
}
eval(m2)
eval(m1)
anova(m1, m2)


