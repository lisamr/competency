#chlamydospores estimations models for both conifers and broadleaves

#starting new script to analyze chlamydos. using a multilevel model with a observation level random effect to account for the dispersion and a poisson likelihood. I might be running into overfitting issues since the fit is good, but the predictions suck.
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

#load plotting functions for these models
#model fit
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
  lesionmin <- min(data$lesion)
  draws <- data %>%
    group_by(species) %>%
    data_grid(
      lesion = seq(lesionmin,3,length.out = 100),
      area_samp=1,
      leafID=686) %>%
    mutate(broad=ifelse(species %in% c('PIPO', 'PSME', 'SESE', 'LIDED', 'UMCAD'), 0, 1)) %>% 
    add_fitted_draws(model, n = 100, re_formula = ~(lesion|species)) 
  #plot
  ggplot(draws, aes(x = lesion,y = .value)) +
    geom_line(aes(group = paste(species, .draw)), alpha = .1)+
    facet_wrap(~species)+
    geom_point(data = data, aes(lesion, counte+.01), 
               color='steelblue', alpha=.4)+
    scale_y_continuous(limits = c(0,5000))
}

#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
options(scipen = 999) #turns off scientific notation
theme_set(theme_light(base_size = 12)) #set ggplot theme


#for easy reporting predicted values with ci and sd
median_hdi_sd <- function(data, value=.value, width=.9){
  data %>% 
    summarise(estimate=median(.value), 
              lower=hdci(.value, .width = width)[1], 
              upper=hdci(.value, .width = width)[2],
              sd=sd(.value),
              width=width)
}

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf chl. treatment assay
df1 <- df %>% 
  filter(trt=='T', spore_assay=="C", !is.na(count1), species !="SESE") %>% 
  dplyr::select(-count2, -count3) %>% 
  rename(count=count1) %>% 
  mutate(area_samp=prop*leaf_area_cm2) %>% 
  droplevels() 

#plot data (estimated counts due to offsets)
df1 %>% mutate(count_est=count/area_samp) %>% 
  ggplot( aes(lesion_area_cm2, log(count_est)))+
  geom_point(alpha=.4)+
  facet_wrap(~species, scales = 'free')


################################################
#model it!
#run seperate models for each assay.

#remove hear and pipo because they seems to be causing problems and we know it's zero.
df2 <- filter(df1, !species %in% c("HEAR", "PIPO")) %>% 
  select(species, ind, leafID, count, area_samp, prop, lesion_area_cm2, leaf_area_cm2) %>% 
  rename(lesion=lesion_area_cm2) %>% 
  mutate(broad=ifelse(leafID>530, 0, 1),
         counte=count/area_samp,
         species=droplevels(species))

#all together model
f2 <- bf(count ~ lesion + broad + (lesion|species) + (1|leafID) + offset(log(area_samp)), family = poisson())
get_prior(f2, df2)
prior2 <- c(
  set_prior("normal(2, 3.5)", class = "Intercept"), 
  set_prior("normal(0, 2)", class = "b", coef = 'lesion'),
  set_prior("normal(0, .5)", class = "b", coef = 'broad'), #assay
  set_prior("exponential(1)", class = "sd")) #all the variations
m4 <- brm(formula = f2, prior = prior2,
           data = df2, family = poisson(),
           chains = 4, cores = 4, iter = 3000,
           control = list(adapt_delta = .99, max_treedepth=15))
#quick checks on model
mod=m4
dat=df2
summary(mod)

plot(mod, pars = "^b_")#looks good
#visualization
pdf('plots/lesions/multiple_chlreg.pdf')
psims(dat, mod)#model fit
pcoef(mod)#species coefs
fpred(df2, m4)#predictions
dev.off()

################################################
#centering/scaling and using fixed effects slopes

#what does centering lesion by group(species) do? I know you typically scale by taking the mean of all lesions, but i want to do it with the mean of the lesion_sp[i]
#df3 <- df2 %>% group_by(species) %>% mutate(lesion=(lesion-mean(lesion)))#centered
df3 <- df2 %>% group_by(species) %>% mutate(lesion=(lesion-mean(lesion))/sd(lesion))#scaled
#centered/scaled
m5 <- brm(formula = f2, prior = prior2,
          data = df3, family = poisson(),
          chains = 4, cores = 4, iter = 3000,
          control = list(adapt_delta = .99, max_treedepth=15))
#check it out. centering made no difference. try scaling?
summary(m5)
psims(df3, m5)#model fit
pcoef(m5)#species coefs
fpred(df3, m5)

#species slopes are fixed effects
f4 <- bf(count ~ -1 + lesion*species + broad + (1|leafID) + offset(log(area_samp)), family = poisson())
get_prior(f4, df3)
exp(rnorm(10000, 2, 3.5)) %>% rethinking::dens(xlim=c(0,20000))
prior3 <- c(
  set_prior("normal(2, 3.5)", class = "b"), #species int 
  set_prior("normal(0, 2)",class="b",coef='lesion:speciesARME'),
  set_prior("normal(0, 2)",class="b",coef='lesion:speciesCEOL'),
  set_prior("normal(0, 2)",class="b",coef='lesion:speciesLIDE'),
  set_prior("normal(0, 2)",class="b",coef='lesion:speciesLIDED'),
  set_prior("normal(0, 2)",class="b",coef='lesion:speciesPSME'),
  set_prior("normal(0, 2)",class="b",coef='lesion:speciesTODI'),
  set_prior("normal(0, 2)", class = "b", coef = 'lesion'),
  set_prior("normal(0, .5)", class = "b", coef = 'broad'), #assay
  set_prior("exponential(1)", class = "sd")) #all the variations
m6 <- brm(formula = f4, prior = prior3,
          data = df2, family = poisson(),
          chains = 4, cores = 4, iter = 3000,
          control = list(adapt_delta = .99, max_treedepth=15))

#check out model
summary(m6)
#see model fit
psims(df2, m6) #looks good
#slope coef plot
post <- m6 %>% 
  spread_draws(b_lesion, `b_lesion:speciesARME`, `b_lesion:speciesCEOL`, `b_lesion:speciesLIDE`, `b_lesion:speciesLIDED`, `b_lesion:speciesPSME`, `b_lesion:speciesTODI`) %>% 
  transmute(.chain, .iteration, .draw,
        beta_ACMA = b_lesion,
         beta_ARME = b_lesion+`b_lesion:speciesARME`,
         beta_CEOL = b_lesion+`b_lesion:speciesCEOL`,
         beta_LIDE = b_lesion+`b_lesion:speciesLIDE`,
         beta_LIDED = b_lesion+`b_lesion:speciesLIDED`,
         beta_PSME = b_lesion+`b_lesion:speciesPSME`,
         beta_TODI = b_lesion+`b_lesion:speciesTODI`) %>% 
  tidyr::gather("slope", "value", -.chain, -.iteration, -.draw) %>% group_by(slope)
post %>% median_hdci(value, .width = .9) 
post %>% 
  ggplot(aes(y = fct_rev(slope), x = value)) +
  geom_halfeyeh(.width=.9, size=.5,point_interval=median_hdcih) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='species specific slopes', title = 'species = fixed effects') 
#see predictions 
fpred(df2, m6)+scale_y_log10()

#try out the fixed effect model with only a subset of the data. include only species where you know the slopes should be around zero. #LIDE, PSME
#yeah it looks more in line with what I'm expecting. Not sure why the interaction slope parameters (lesion:species) aren't stronger when more species are included.
df4 <- df2 %>% filter(species %in% c("LIDE", "LIDED", "PSME")) %>% droplevels()
prior4 <- c(
  set_prior("normal(2, 3.5)", class = "b"), #species int 
  set_prior("normal(0, 2)",class="b",coef='lesion:speciesLIDED'),
  set_prior("normal(0, 2)",class="b",coef='lesion:speciesPSME'),
  set_prior("normal(0, 2)", class = "b", coef = 'lesion'),
  set_prior("normal(0, .5)", class = "b", coef = 'broad'), #assay
  set_prior("exponential(1)", class = "sd")) #all the variations
m7 <- brm(formula = f4, prior = prior4,
          data = df4, family = poisson(),
          chains = 4, cores = 4, iter = 3000,
          control = list(adapt_delta = .99, max_treedepth=15))
#species slopes 
post <- m7 %>% 
  spread_draws(b_lesion, `b_lesion:speciesLIDED`, `b_lesion:speciesPSME`) %>% 
  transmute(.chain, .iteration, .draw,
            beta_LIDE = b_lesion,
            beta_LIDED = b_lesion+`b_lesion:speciesLIDED`,
            beta_PSME = b_lesion+`b_lesion:speciesPSME`) %>% 
  tidyr::gather("slope", "value", -.chain, -.iteration, -.draw) %>% group_by(slope)
post %>% median_hdci(value, .width = .9) 
#see predictions 
fpred(df4, m7)+scale_y_continuous(limits=c(0,350))+scale_x_continuous(limits = c(0,1))
fpred(df4, m7)+scale_y_log10(limits=c(0.001,350))+scale_x_continuous(limits = c(0,1))

#ok, going to remove the population level slope variable for lesion. can I do that?
f5 <- bf(count ~ broad + (lesion|species) + (1|leafID) + offset(log(area_samp)), family = poisson())
get_prior(f5, df2)
prior5 <- c(
  set_prior("normal(2, 3.5)", class = "b"), #species int 
  set_prior("normal(0, .5)", class = "b", coef = 'broad'), #assay
  set_prior("exponential(1)", class = "sd")) #all the variations
m8 <- brm(formula = f5, prior = prior5,
          data = df2, family = poisson(),
          chains = 4, cores = 4, iter = 3000,
          control = list(adapt_delta = .99, max_treedepth=15))
#check out model
summary(m8)
fpred(df2, m8)
fpred(df2, m8)+scale_y_log10()
psims(df2, m8)
summary(m8)
#CEOL and TODI only significantly positive slopes in this model
coefs <- m8 %>% 
  spread_draws(b_broad, r_species[species,term]) %>%
  mutate(slope=r_species + b_broad) %>% 
  filter(term=="lesion")
coeftable <- coefs %>% median_hdci(slope, .width = .9) 
#plot
coefs %>% 
  ggplot(aes(y = fct_rev(species), x = slope)) +
  geom_halfeyeh(.width = .9, size=.5, point_interval = median_hdcih) +
  geom_point(data = coeftable, aes(slope, species), color='steelblue', size=2) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient') 

#save all your models. can recover with readRDS()
saveRDS(m4, file = 'output/models/ch_lesion_m4.rds')
saveRDS(m6, file = 'output/models/ch_lesion_m6.rds')
saveRDS(m7, file = 'output/models/ch_lesion_m7.rds')
saveRDS(m8, file = 'output/models/ch_lesion_m8.rds')

####################################
#seperate models
#subset data for seperate models. 
df3b <- filter(df2, leafID<=530) %>% droplevels()
df3c <- filter(df2, leafID>530) %>% droplevels()

#formulas and priors. priors same as chlamydo sporulation model
f1 <- bf(count ~ lesion + (lesion|species) + (1|leafID) + offset(log(area_samp)), family = poisson())
prior1 <- c(
  set_prior("normal(2, 3.5)", class = "b"), #species int 
  set_prior("exponential(1)", class = "sd")) #all the variations
#models
m3b <- brm(formula = f1, prior = prior1,
           data = df3b, family = poisson(),
           chains = 4, cores = 4, iter = 3000,
           control = list(adapt_delta = .99, max_treedepth=15))
m3c <- brm(formula = f1, prior = prior1,
           data = df3c, family = poisson(),
           chains = 4, cores = 4, iter = 3000,
           control = list(adapt_delta = .99, max_treedepth=15))
