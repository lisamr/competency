#sporangia analysis for both assays--broadleaf and conifers. pull the posteriors from the two models, standardize the units and plot the results together.

setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(tidyverse) #includes dplyr, ggplot, forcats
library(reshape2)
library(scales)
library(ggridges)
library(rethinking)
library(tidybayes)
library(brms)
library(ggthemes)

source("rscripts/analysis2019/sporangia_lesion_plotfuncs.R")#load plotting functions for these models
reset <- function(x) par(mfrow=c(1,1))#turns pars back 
theme_set(theme_light(base_size = 12)) #set ggplot theme

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')

#check out lesion area and contrast with the controls
df %>% filter(spore_assay=="S") %>% 
  group_by(species, trt) %>% 
  summarise(lesion=mean(perc_lesion, na.rm = T)) %>% 
  spread(trt, lesion)
df %>% filter(spore_assay=="S") %>% 
  rowwise() %>% 
  mutate(countm=mean(c(count1, count2, count3), na.rm = T)) %>% ungroup() %>% 
  ggplot(., aes(species, perc_lesion))+
  geom_point(aes(color=log(countm), alpha=.3))+
  scale_color_viridis_c()+
  facet_grid(~trt)
   
#analyze just trt, sporangia
df1 <- df %>% 
  filter(trt=="T", spore_assay=="S") 

#viz sporangia data
#start with averages
df1<- df1 %>% rowwise() %>% mutate(countm=mean(c(count1, count2, count3), na.rm = T) ) #rowwise "groups" each row so each mean calucation is unique by row

#plot it. 
#raw counts
df1 %>% 
  ggplot( aes(lesion_area_cm2, countm))+
  geom_point(alpha=.4)+
  facet_wrap(~species, scales = 'free')
df1 %>% 
  ggplot( aes(species, (countm/leaf_area_cm2)))+
  geom_boxplot()+
  geom_point(alpha=.4)+
  scale_y_log10()
df1 %>% 
  ggplot( aes(lesion_area_cm2, countm))+
  geom_point(alpha=.4)

#get summary statistics for lesion size by species
df1 %>% group_by(species) %>% 
  filter(!is.na(lesion_area_cm2)) %>% 
  summarise(mean(lesion_area_cm2), sd(lesion_area_cm2))
df1 %>% group_by(species) %>% 
  filter(!is.na(perc_lesion)) %>% 
  summarise(mean(perc_lesion), sd(perc_lesion)) %>% 
  mutate_if(is.numeric, signif, 3) #%>% 
  #write_csv('output/lesion/summary_lesionsize.csv')
ggplot(df1, aes(perc_lesion, species)) +
  geom_halfeyeh(.width = .9, size = .1)
df1 %>% filter(perc_lesion<1) %>% 
  select(species, ind, perc_lesion, count1, count2, count3) %>% 
  print(n = Inf)
#not sure the best way to model:
#1. broadleaf and conifers in separate models
#2. model together with a fixed effect for assay
#3. don't bother with the fixed effect

#trying option 1 for now.
#make data tall in regards to counts
dtall2 <- melt(df1, id.vars = c("species", "leafID", "lesion_area_cm2"), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count), !is.na(lesion_area_cm2)) %>%
  rename(lesion=lesion_area_cm2) %>% 
  arrange(leafID) 
#separate by assay
dtallb <- dtall2 %>% filter(leafID<=530) %>% droplevels()
dtallc <- dtall2 %>% filter(leafID>530) %>% droplevels()

#separate models first
#set prior first
f1 <- bf(count ~ lesion + (1|leafID) + (lesion|species), family = poisson())
get_prior(f1, dtallb)
prior1 <- c(
  set_prior("normal(2, 2)", class="Intercept"),
  set_prior("normal(0, 20)", class="b"),
  set_prior("exponential(1)", class = "sd") )

#models
m1b <- brm(
  f1, data = dtallb, family = poisson(), prior = prior1,
  chains = 4, cores = 4, iter=3000,
  control = list(adapt_delta = .99, max_treedepth=15))
#no warnings :)
m1c <- brm(
  f1, data = dtallc, family = poisson(), prior = prior1,
  chains = 4, cores = 4, iter=3000,
  control = list(adapt_delta = .99, max_treedepth=15))
#no warnings :)

#traceplots
#broadleaf model
plot(m1b, pars = "^b_")#looks good
#sink('output/lesion/summary_broad.txt')
summary(m1b) 
#sink()
#check posterior with plots. functions from other script. source code called at beginning.
psims(dtallb, m1b)#good
#pdf('plots/lesions/ranef_coefplot_broad.pdf', 8, 6)
pcoef(m1b)
#dev.off()
#pcoef(m1b)[[1]] %>%  mutate_if(is.numeric, signif, 3) %>%write_csv("output/lesion/randcoef_sporangia_broad.csv")
#pdf('plots/lesions/predict_broad.pdf', 8, 6)
ppredict(dtallb, m1b)
#dev.off()

#conifer model
plot(m1c, pars = "^b_")
#sink('output/lesion/summary_conifer.txt')
summary(m1c) #sd_lesion kinda high
#sink()
psims(dtallc, m1c)
#pdf('plots/lesions/ranef_coefplot_conifer.pdf', 8, 6)
pcoef(m1c)
#dev.off()
#pcoef(m1c)[[1]] %>% mutate_if(is.numeric, signif, 3)%>% write_csv("output/lesion/randcoef_sporangia_conifer.csv")
#pdf('plots/lesions/predict_conifer.pdf', 10, 3)
ppredict(dtallc, m1c, prow = 1, pcol = 5)
#dev.off()

#get values of other coefs 
parnames(m1b)
m1b %>% 
  gather_draws(b_Intercept, b_lesion, sd_leafID__Intercept,sd_species__Intercept,sd_species__lesion, cor_species__Intercept__lesion) %>% 
  median_hdci(.width = .9) %>% 
  mutate_if(is.numeric, signif, 3)#%>% 
  #write_csv('output/lesion/fixedcoef_sporangia_broad.csv')
m1c %>% 
  gather_draws(b_Intercept, b_lesion, sd_leafID__Intercept,sd_species__Intercept,sd_species__lesion, cor_species__Intercept__lesion) %>% 
  median_hdci(.width = .9)%>% 
  mutate_if(is.numeric, signif, 3)#%>% 
  #write_csv('output/lesion/fixedcoef_sporangia_conifer.csv')

#should you multiply the slope coefs by 102.5 so you can say something about how each cm2 lesion affects the counts on the leaf (not the sample)?


######
#other models below where I try to model both assays together. Could compare the coefs, but don't think its necessary. 

#throw everything in there toghther.
m2 <- brm(
  f1, data = dtall2, family = poisson(),
  chains = 3, cores = 4, iter=2000,
  control = list(adapt_delta = .95, max_treedepth=15))
#16 dt, .95, no custom priors

#check it out.
plot(m2, pars = "^b_")
summary(m2) 
psims(dtall2, m2)
pcoef(m2)
ppredict(dtall2, m2, prow = 3, pcol = 5)

#add in a fixed effect for assay
#add in assay
dtall3 <- melt(df1, id.vars = c("species", "leafID", "lesion_area_cm2"), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count), !is.na(lesion_area_cm2)) %>%
  rename(lesion=lesion_area_cm2) %>% 
  arrange(leafID) %>% 
  mutate(detached=if_else(leafID >530, 1, 0))
#add in fixed effect
f2 <- bf(count ~ lesion + detached + (1|leafID) + (lesion|species), family = poisson())
#model
m3 <- brm(
  f2, data = dtall3, family = poisson(),
  chains = 3, cores = 4, iter=2000,
  control = list(adapt_delta = .95, max_treedepth=15))
#16 dt, .95, no custom priors

#check it out.
plot(m3, pars = "^b_")
summary(m3) #sd_lesion crazy high
psims(dtall3, m3)
pcoef(m3)
ppredict(dtall3, m3, prow = 3, pcol = 5)

#compare m2 and m3. m1 can't be compared since its different data.
m2 <- add_criterion(m2, criterion = c('loo', 'waic'))
m3 <- add_criterion(m3, criterion = c('loo', 'waic'))
loo_compare(m2, m3, criterion = 'waic') #m3 not any better so the one without the fixed effect would be fine. Still think seperate models completely is appropriate.

