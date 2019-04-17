#how do you analyze data that varies in effect direction at different scales? Ex/across all species, the relationship between sporangia and chlamydospres is negative, but within species, it is positive?
#create fake data to explore it

library(dplyr)
library(ggplot2)
library(lme4)
library(emmeans)

#relationship within each of the species is positive, but the relationship among species is negative
N=9
meanC <- sample(100:1000, N)
f <- function(x, n){
  spC <- rnorm(n, meanC[x], 50) %>% round #generate n samples with chlamydospore counts normally distr. around the mean
  mu <- rnorm(1, 20*mean(spC), 1000) #the species intercepts are positively relatived to the mean chlamydospore counts
  spS <- rnorm(n,mu + (-10)*spC, 1000) %>% round #generate n samples with sporangia counts normally distr. around the mean (int + x*chlam)
  species <- rep(LETTERS[x], length(spS))
  data.frame(species, spS, spC)
}
df <- lapply(1:N,function(x) f(x,32)) %>% bind_rows()
df1 <- df %>% group_by(species) %>%  summarise(S=mean(spS), C=mean(spC))
ggplot(df, aes(spC, spS, group=species, color=species))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  geom_point(data=df1, aes(C, S), color="black")
ggplot(df, aes(spC, spS, group=species, color=species))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  facet_wrap(~species, scales = "free")

#how would you analyze this? 
head(df)
m1 <- lm(spS ~ spC*species, data=df)
summary(m1) #this model says that at the species level, there's neg relationship between C and S (-5) and the species estimates refer to the differences in intercepts (mean values of S)
emtrends(m1, pairwise ~ species, var = "spC")
#how do you measure the relationship between S and C across all species? The output still says it's a neg. relationship after treating species as a random effect. Do you need to manually aggregate the data?
m2 <- lmer(spS ~ spC + (1|species), data=df)
summary(m2)
cim2 <- confint.merMod(m2, method = 'boot')

#data grouped by species. this would work if each species has the same number of individuals, but what about an unbalanced design? surely you would want to weight the ones with more samples.
m3 <- lm(S ~ C, df1)
summary(m3)
ggplot(df1, aes(C, S))+
  geom_point()+
  geom_smooth(method='lm')



