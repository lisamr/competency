#sporangia analysis for both assays--broadleaf and conifers. pull the posteriors from the two models, standardize the units and plot the results together.
setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(tidybayes)
library(modelr)
library(brms)
library(ggridges)
library(ggstance)
library(forcats)

#set ggplot theme
theme_set(theme_bw(base_size = 9) + theme(panel.grid=element_blank()))
#functions
inv_logit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x/(1-x))
median_hdi_sd <- function(data, value=.value, width=.9){
  value <- enquo(value)
  data %>% 
    summarise(estimate=median(!!value), 
              lower=hdci(!!value, width)[1], 
              upper=hdci(!!value, width)[2],
              sd=sd(!!value),
              width=.9)
}
#plotting functions
psims <- function(data, model) {
  #sample from posterior
  sims<- add_predicted_draws(data, model) %>% 
    group_by(species, .draw) %>% 
    sample_draws(30) 
  #plot draws against data
  sims %>% ggplot() +
    geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.2, color='grey')+
    stat_density(data=data, aes(x=count), geom="line", color='steelblue') +
    facet_wrap(~species, scales = 'free')
}
#coef plot
pcoef <- function(model){
  #coef plot
  coefs <- model %>% 
    spread_draws(b_lesion, r_species[species,term]) %>%
    mutate(slope=r_species + b_lesion) %>% 
    filter(term=="lesion")
  coeftable <- coefs %>% median_hdci(slope, .width = .9) 
  #plot
  PLOT <- coefs %>% 
    ggplot(aes(y = fct_rev(species), x = slope)) +
    geom_halfeyeh(.width = .9, size=.5, point_interval = median_hdcih) +
    geom_point(data = coeftable, aes(slope, species), color='steelblue', size=2) +
    geom_vline(xintercept=0, lty=2, color='grey50') +
    labs(y='coefficient') 
  return(list(coeftable, PLOT))
}
#posterior predictions
fpred <- function(data, model){
  #get post pred
  draws <- data %>%
    group_by(species) %>%
    data_grid(
      lesion = seq(0,1.23,length.out = 100),
      area_samp=1) %>%
    mutate(broad=ifelse(species %in% c('PSME', 'LIDED', 'SESE', 'PIPO', 'UMCAD'), 0, 1)) %>% 
    add_fitted_draws(model, n = 100, re_formula = ~(lesion|species)) 
  #plot
  ggplot(draws, aes(x = lesion,y = .value)) +
    geom_line(aes(group = paste(species, .draw)), alpha = .1)+
    facet_wrap(~species, scales = 'free_y')+
    geom_point(data = data, aes(lesion, countm), 
               color='steelblue', alpha=.4) +
    labs(x=expression(paste("Lesion area (", cm^{2}, ")")),
         y="Number of sporangia")
}
#############################################
#############################################
#source("rscripts/analysis2019/lesion_plotfuncs.R")#load plotting functions for these models

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just trt, sporangia
df1 <- df %>% 
  filter(trt=="T", spore_assay=="S") 
#start with averages
df1<- df1 %>% rowwise() %>% mutate(countm=mean(c(count1, count2, count3), na.rm = T) ) #rowwise "groups" each row so each mean calucation is unique by row

#make data tall in regards to counts
dtall <- reshape2::melt(df1, id.vars = c("species", "leafID", "lesion_area_cm2", 'countm'), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>%
  filter(!is.na(count), !is.na(lesion_area_cm2)) %>%
  rename(lesion=lesion_area_cm2) %>% 
  arrange(leafID) %>% 
  mutate(broad = ifelse(leafID>530, 0, 1),
         counte = countm*102.5)

#read final model
m1 <- readRDS('output/models/spor_lesion_m1.rds')

#viz sporangia data
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

#going with option 2
f1 <- bf(count ~ lesion + broad + (1|leafID) + (lesion|species), family = poisson())
get_prior(f1, dtall)
prior1 <- c(
  set_prior("normal(1.5, 1.5)", class="Intercept"),
  set_prior("normal(0, 20)", class = "b", coef = 'lesion'),
  set_prior("normal(0, .5)", class = "b", coef = 'broad'),
  set_prior("exponential(1)", class = "sd") )

#models
m1 <- brm(
  f1, data = dtall, family = poisson(), prior = prior1,
  chains = 4, cores = 4, iter=3000,
  control = list(adapt_delta = .99, max_treedepth=15))

#check out model
summary(m1)
psims(dtall, m1)+scale_x_continuous(limits=c(0,500))
pp_check(m1, type = "intervals_grouped", group = "species")#pretty bad fit for some species
pcoef(m1)#looks similar to the zinb

#prediction plot. have to manually change the y scale limits
draws <- dtall %>%
  group_by(species) %>%
  data_grid(
    lesion = seq(0,1.3,length.out = 100)) %>%
  mutate(broad=ifelse(species %in% c('PSME', 'LIDED', 'SESE', 'PIPO', 'UMCAD'), 0, 1)) %>% 
  add_fitted_draws(m1, n = 100, re_formula = ~(lesion|species)) %>% 
  mutate(.value2=.value*102.5)
#plot
ylims <- c(15, 2.5, 15, 5, 100, 100, 2.5, 15, 15, 30, 15, 15, 15, 100, 30)*102.5
lookup <- data.frame(species=unique(draws$species), ylims)

#those that need asterisks: acma, todi, lided, umcad, quch

#get shades and zoom in without losing data
#http://www.zachburchill.ml/ggplot_facets/
source("rscripts/analysis2019/zoom_facets.R")
plot_shade <- draws %>% 
  inner_join(lookup, by = 'species') %>%
  select(-.chain, -.iteration) %>% 
  ggplot(., aes(x = lesion,y = .value2)) +
  stat_lineribbon(aes(y = .value2), .width = .9, fill='grey', lwd=.5, color='gray20', alpha=.5) +
  #add in asterisks for the 'significant' slopes
  facet_wrap(~species, scales='free_y', nrow = 3, ncol = 5, labeller = labeller(
    species=c(ACMA="ACMA*", 
              ARME="ARME",
              CEOL='CEOL',
              HEAR='HEAR',
              LIDE='LIDE',
              LIDED='LIDED*',
              PIPO='PIPO',
              PSME='PSME',
              QUAG='QUAG',
              QUCH='QUCH*',
              QUPA='QUPA',
              SESE='SESE',
              TODI='TODI*',
              UMCA='UMCA',
              UMCAD='UMCAD*'))) +
  geom_point(data =filter(dtall, !(species=="ACMA"&countm>100)), aes(lesion, counte), color='steelblue', alpha=.2) +
  labs(x=expression(paste("Lesion area (", cm^{2}, ")")),
       y="Number of sporangia") +
  scale_x_continuous(breaks=c(0, .5, 1)) + 
  #define x and y limits for each panel
  coord_panel_ranges((
    panel_ranges=lapply(1:length(ylims), function(i)list(x=c(-.05, 1.3), y=c(-10, ylims[i])))
    ))
plot_shade



################################################
#saving model, outputs, and figures
saveRDS(m1, file = 'output/models/spor_lesion_m1.rds')
#saving outputs
sink('output/lesion/spor_summary_lesion.txt')
summary(m1, prob = .9)
sink()
pdf('plots/lesions/spor_multiple_m1.pdf')
psims(dtall, m1)
pp_check(m1, type = "intervals_grouped", group = "species")#huge improvement
PREDplot
pcoef(m1)
dev.off()
ggsave('plots/lesions/spor_predictions.jpg', PREDplot, width = 7, height = 5.5, dpi = 600, units = 'in')
ggsave('plots/lesions/spor_predictions_tall.jpg', PREDplot, width =3.5, height = 5, dpi = 600, units = 'in')
ggsave('plots/lesions/spor_predictions_wide.jpg', plot_shade, width =7, height = 4.5, dpi = 600, units = 'in')
slopes <- pcoef(m1)
slopes[[1]] %>% mutate_if(is.numeric, signif, 3) %>% 
  write_csv('output/lesion/spor_randcoef.csv') 

#model coefficents into dataframe
vars <- get_variables(m1)[1:14]
m1 %>% 
  gather_draws(!!!syms(vars)) %>% 
  median_hdi_sd() %>% 
  mutate_if(is.numeric, signif, 3) %>% 
  write_csv('output/lesion/spor_coef_estimates.csv')


################################################
#############################################
#older stuff

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

