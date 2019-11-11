setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(tidyverse) #includes dplyr, ggplot, forcats
library(reshape2)
library(scales)
library(ggridges)
library(rethinking)
library(tidybayes)
library(brms)

#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
theme_set(theme_bw()) #set ggplot theme

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just conifer assay, trt, sporangia
df1 <- df %>% 
  filter(leafID>=531, trt=="T", spore_assay=="S") 

#viz sporangia data
#start with averages
df1<- df1 %>% rowwise() %>% mutate(countm=mean(c(count1, count2, count3), na.rm = T) ) #rowwise "groups" each row so each mean calucation is unique by row

#plot it. 
#raw counts
df1 %>% 
  ggplot( aes(species, countm))+
  geom_boxplot()+
  geom_point(alpha=.4)+
  scale_y_log10()
df1 %>% 
  ggplot( aes(species, log(countm/leaf_area_cm2)))+
  geom_boxplot()+
  geom_point(alpha=.4)
#approximate spores/cm2
df1 %>% 
  ggplot( aes(species, (countm/leaf_area_cm2*102.5)))+
  geom_boxplot()+
  geom_point(alpha=.4)+
  scale_y_log10()

df1[is.na(df1$countm),] #the samples without counts were too low in vol or too goopy


################################################
################################################
#model sporangia counts. include leaf area as an offset or standardize the posteriors after modeling? I think it should be an offset.

#make data tall in regards to counts
dtall <- melt(df1, id.vars = c("species", "leafID", 'leaf_area_cm2'), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count), !is.na(leaf_area_cm2)) %>% 
  arrange(leafID) 
dtall$species <- droplevels(dtall$species)
head(dtall); tail(dtall)

#figure out priors
#species int. 0-1000 sp/cm2? based on larch study. larch was above (~2000) but other species were way below. divide that by 102.5: 0-10 mostly
#same prior as broadleaf model?
dens(exp(rnorm(1000, 2, 1.5)), xlim=c(0,100)) 

#model formula
f1 <- bf(count ~ -1 + species + (1|leafID) + offset(log(leaf_area_cm2)), family = poisson())

#check out default priors
get_prior(f1, dtall)

#set your own
prior1 <- c(
  set_prior("normal(2, 1.5)", class = "b"),
  set_prior("exponential(1)", class = "sd") )

#model finally
m1 <- brm(
  formula = f1, data = dtall, family = poisson(),
  chains = 4, cores = 4, prior = prior1,
  control = list(adapt_delta = .95, max_treedepth=15))
#no warnings :)

#check out coefs
summary(m1)
fixef(m1)
#population level effects only. traceplots look good
plot(m1, pars = "^b_")
post <- posterior_samples(m1, pars = c("^b_", "sd"))
head(post)
pp_check(m1)

#see model fit. seems good.
#simulate prediction data using same dataset to see fit
m1sim <- add_predicted_draws(dtall, m1) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(100) 
#plot against your own data
m1sim %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.5, color='grey')+
  stat_density(data=dtall, aes(x=count), geom="line", color='slateblue')+
  facet_wrap(~species, scales = 'free')

#coef plot
m1 %>% 
  gather_draws(b_speciesLIDED, b_speciesPIPO, b_speciesPSME, b_speciesSESE, b_speciesUMCAD, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')

#get predicted fit based on just the species coef.
postp <- dtall %>%
  modelr::data_grid(species, leaf_area_cm2=1) %>%
  add_fitted_draws(m1, re_formula = ~0, scale = 'response') %>%
  mutate(.valueSTD = .value*102.5) #standardized sp/cm2

#get the summarized values
postp %>% 
  median_hdci(x=.valueSTD, .width = c(.9, .95))

#plot them
ggplot(postp, aes(y = fct_rev(species), x = .valueSTD)) +
  geom_density_ridges(lwd=.1)+
  stat_pointintervalh(point_interval = median_hdcih, .width = c(.9), size=.2)+
  labs(x='# sporangia', y='Species')


