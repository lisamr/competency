#starting new script to analyze chlamydos. using a multilevel model with a observation level random effect to account for the dispersion and a poisson likelihood. Im hoping i can ignore the zero inflation but i dont think i can.
setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(scales)
library(rethinking)
library(tidybayes)
library(brms)
library(ggridges)
library(forcats)

library(dplyr)
library(tidyr)
library(cowplot)
library(ggstance)

reset <- function(x) par(mfrow=c(1,1))#for turning pars back 
options(scipen = 999) #turns off scientific notation
theme_set(theme_bw())

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf chl. treatment assay
df1 <- df %>% filter(leafID<=530, trt=='T', spore_assay=="C", !is.na(count1)) %>% select(-count2, -count3) %>% rename(count=count1) %>% mutate(perc_lesion=perc_lesion/100)

#plot data (estimated counts due to offsets)
df1 %>% mutate(count_est=count/prop) %>% 
  ggplot( aes(perc_lesion, (count_est)))+
  geom_point(alpha=.4) +
  facet_wrap(~species, scales = 'free')

########################################################
#model count~lesion with offsets and RE

#remove HEAR cuz all the zeros are causing problems
df2 <- filter(df1, !species=="HEAR")
df2$species <- droplevels(df2$species)

#model formula
f1 <- bf(count ~ perc_lesion + (1|leafID) + (perc_lesion|species) + offset(log(prop)), family = poisson())

#check out default priors
get_prior(f1, df2)

#set your own
prior1 <- c(
  set_prior("normal(0,1)", class = "b", coef='perc_lesion'),
  set_prior("normal(0,1)", class = "Intercept"),
  set_prior("exponential(1)", class = "sd") )
prior1.1 <- c(
  set_prior("normal(0,1)", class = "Intercept"),
  set_prior("exponential(1)", class = "sd") )

#model it
m1 <- brm(
  count ~ perc_lesion + (1|leafID) + (perc_lesion|species) + offset(log(prop)),
  data = df2, family = poisson(), prior = prior1,
  chains = 3, cores = 4, sample_prior = T,
  control = list(adapt_delta = .95, max_treedepth=15))
#sampled much better after restricting the priors a bit

#try without the leafID random effect. hopefully we can "assign" more of the response to the lesion size rather than having a random effect for overdispersion soak it up.
m2 <- brm(
  count ~ perc_lesion + (perc_lesion|species) + offset(log(prop)),
  data = df2, family = poisson(), prior = prior1,
  chains = 3, cores = 4, sample_prior = T,
  control = list(adapt_delta = .95, max_treedepth=15))

#should go with m1
m1 <- add_loo(m1)
m2 <- add_loo(m2)
loo_compare(m1, m2)

summary(m1)
plot(m2, pars = "^b_")#looks good
pp_check(m1) #one with the the leafID RE fit better

#see model fit. seems pretty good.
#simulate prediction data using same dataset to see fit
msim <- df2 %>% 
  select(species, leafID, perc_lesion, count, prop) %>% 
  add_predicted_draws(m1) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(50) 
msim %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.5, color='grey')+
  stat_density(data=df2, aes(x=count), geom="line", color='slateblue')+
  facet_wrap(~species, scales = 'free')

#coefficient plot
parnames(m1)
m1 %>% 
  gather_draws(b_Intercept, b_perc_lesion, sd_leafID__Intercept, sd_species__Intercept, sd_species__perc_lesion, cor_species__Intercept__perc_lesion, r_species[species,]) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')

#check out species-specific slopes and ints
m1.1 %>% 
  spread_draws(r_species[species,term]) %>%
  ggplot(aes(y = fct_rev(species), x = r_species)) +
  geom_halfeyeh() +
  facet_grid(~term)

#check out model predictions
xx <- seq(0,1,length.out = 100)
newd <- data.frame(perc_lesion=xx, leafID=1, prop=1, species=rep(unique(df2$species), each=length(xx)))
pp1 <- fitted(m1, newdata = newd, re_formula = ~(1|species), summary = F, scale='response')
pp2 <- fitted(m2, newdata = newd, re_formula = ~(1|species), summary = F, scale='response')

plotmodfit <- function(i, pp, data=df2, ...){
  sp <- as.vector(unique(newd$species))
  use <- as.integer(factor(newd$species))==i
  #estimate observed counts with offset
  data$counte <- data$count/data$prop
  plot(data$perc_lesion[data$species==sp[i]], data$counte[data$species==sp[i]], ...,
       xlab='lesion', ylab='# chamydospores', main=sp[i], 
       pch=16, col=alpha('slateblue', .5), xlim=c(0,1))
  apply(pp[,use], 2, median) %>% lines(xx, .)
  apply(pp[,use], 2, PI, .9) %>% shade(., xx)
}
#these are pretty shitty predictions
par(mfrow=c(2,5))
sapply(1:length(unique(newd$species)), function(x) plotmodfit(x, pp1))
sapply(1:length(unique(newd$species)), function(x) plotmodfit(x, pp2))
reset()
plotmodfit(1, pp1, ylim=c(0,4000))
plotmodfit(2, pp1)
plotmodfit(3, pp1)
plotmodfit(4, pp1)
plotmodfit(5, pp2, ylim=c(0,10))
