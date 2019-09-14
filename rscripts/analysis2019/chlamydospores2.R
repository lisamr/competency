#starting new script to analyze chlamydos. using a multilevel model with a observation level random effect to account for the dispersion and a poisson likelihood. Im hoping i can ignore the zero inflation but i dont think i can.
setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)
library(scales)
library(rethinking)
library(tidybayes)
library(modelr)
library(brms)
library(ggridges)

#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
options(scipen = 999) #turns off scientific notation
theme_set(theme_bw())

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf chl. treatment assay
df1 <- df %>% filter(leafID<=530, trt=='T', spore_assay=="C", !is.na(count1)) %>% select(-count2, -count3) %>% rename(count=count1)

#plot data (estimated counts due to offsets)
df1 %>% mutate(count_est=count/prop) %>% 
  ggplot( aes(species, log(count_est)))+
  geom_boxplot()+
  geom_point(alpha=.4)

#model it
#remove hear because it seems to be causing problems and we know it's zero.
df2 <- filter(df1, !species=="HEAR")
df2$species <- droplevels(df2$species)
m1 <- brm(
  count ~ -1 + species + (1|leafID) + offset(log(prop)),
  data = df2, family = poisson(),
  chains = 3, cores = 4, iter = 4000,
  control = list(adapt_delta = .95, max_treedepth=15))
#no warnings :)

#check out coefs
summary(m1)
fixef(m1)
launch_shinystan(m1)
#population level effects only. traceplots look ok
plot(m1, pars = "^b_")#looks good
postm1 <- posterior_samples(m1, pars = c("^b_", "sd"))
head(postm1)
pp_check(m1)

#see model fit. seems pretty good.
#simulate prediction data using same dataset to see fit
m1sim <- df2 %>% 
  select(species, leafID, count, prop) %>% 
  add_predicted_draws(m1) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(50) 
head(m1sim)
#plot against your own data
#pdf('plots/chlamydospores/modelfit_brms.pdf')
m1sim %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.5, color='grey')+
  stat_density(data=df2, aes(x=count), geom="line", color='slateblue')+
  facet_wrap(~species, scales = 'free')
#dev.off()  

#coef plot
#pdf('plots/chlamydos/coefplot_brms.pdf', 6, 4)
parnames(m1)
m1 %>% 
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesLIDE, b_speciesTODI, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')
#dev.off()

#summarize coefs
m1 %>% 
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesLIDE, b_speciesTODI, sd_leafID__Intercept) %>%
  median_qi(.width=c(.9, .95))

#get backtransformed mean counts
#observed estimated means
obsm <- df2 %>% group_by(species) %>% 
  summarise(counte=median(count/prop))
obs <- df2 %>% mutate(counte=count/prop)
head(obs)
ggplot()+
  geom_density_ridges(obs, aes(counte, species), jittered_points = TRUE,position = position_points_jitter(width = 0.05, height = 0), fill=NA, color=NA, point_color='slateblue', point_shape = '|', point_size = 3, point_alpha = .3, alpha = 0.7) #+ scale_x_log10()

#predicted means
pr <- df2 %>%
  data_grid(species) %>% mutate(prop=1)%>%
  add_fitted_draws(m1, re_formula = ~0, scale = 'response') 
pr %>% median_qi(.width=c(.9, .95))
#plot1, just predicted wiht median obs
#pdf('plots/chlamydos/predicted_brms.pdf', 6, 4)
ggplot(pr, aes(y = fct_rev(species), x = .value)) +
  geom_density_ridges(lwd=0, color=NA)+
  stat_pointintervalh(.width = c(.9), size=.2)+
  geom_point(data=obsm, aes(counte, species), fill='slateblue', color='slateblue', shape=24) +
  labs(x='# chlamydospores', y='Species')+
  scale_x_log10()
#dev.off()

#plot it with raw observed data
#pdf('plots/chlamydos/predictedwdata_brms.pdf', 6, 4)
ggplot(pr, aes(y = fct_rev(species), x = .value)) +
  geom_density_ridges(col=NA)+
  stat_pointintervalh(.width = c(.9), size=.2)+
  geom_point(data=obsm, aes(counte, species), fill='slateblue', color='slateblue', shape=24) +
  labs(x='# chlamydospores', y='Species')+
  #observed data
  geom_density_ridges(
    data=obs, aes(counte+.1, species),
    jittered_points = TRUE,position = position_points_jitter(width = 0.05, height = 0), 
    fill='slateblue', alpha=.2, color=NA, 
    point_color='slateblue', point_shape = 'o', point_size = 3, point_alpha = .3)+
  scale_x_log10()
#dev.off()
#Should I be concerned that the model fit doesn't really line up with the observed means? when i do it with medians, it makes more sense and lines up better. 

#contrasts between species
#extract posterior of following coefficients
ex <- m1 %>% 
  spread_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesLIDE, b_speciesTODI) %>% 
  select(-c(.chain, .iteration, .draw)) %>% as.data.frame()
head(ex)

#can't figure out how to do it the tidy way, so doing it kinda clunky

pairs <- combn(1:ncol(ex),2)
f <- function(x){
  spdiff <-  (ex[,pairs[1,x]]-ex[,pairs[2,x]])
  PI(spdiff, .9)
}
diffs <- sapply(1:ncol(pairs), f)
#put contrasts into a readable dataframe
spp <- colnames(ex) %>% str_replace('b_species', '')
diffsdf <- rbind(pairs, diffs) %>% round(2) %>%  t %>% as.data.frame
diffsdf <- diffsdf %>% 
  mutate(sp1=spp[diffsdf[,1]],
         sp2=spp[diffsdf[,2]],
         sig=sign(diffsdf[,3])==sign(diffsdf[,4]))
diffsdf
