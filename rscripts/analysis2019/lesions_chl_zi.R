#uggh why am i still working on this. lets try using a zeroinflated model to deal with the chlamydospores. I read in Harrison (2015, PeerJ) that observation level random effects can be counterproductive for overdispersion caused by zero-inflation. they suggest using a model specific for zi responses.

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
    mutate(broad=ifelse(species %in% c('PSME', 'LIDED'), 0, 1)) %>% 
    add_fitted_draws(model, n = 100, re_formula = ~(lesion|species)) 
  #plot
  ggplot(draws, aes(x = lesion,y = .value)) +
    geom_line(aes(group = paste(species, .draw)), alpha = .1)+
    facet_wrap(~species)+
    geom_point(data = data, aes(lesion, counte+.01), 
               color='steelblue', alpha=.4) +
    scale_y_continuous(limits = c(0,5000))+
    labs(x=expression(paste("Lesion area (", cm^{2}, ")")),
         y="Number of chlamydospores")
}

############################
############################

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf chl. treatment assay
df1 <- df %>% 
  filter(trt=='T', spore_assay=="C", !is.na(count1), !is.na(lesion_area_cm2), !species %in%c("SESE", "HEAR", "PIPO")) %>% 
  dplyr::select(-count2, -count3) %>% 
  rename(count=count1, lesion=lesion_area_cm2) %>% 
  mutate(area_samp=prop*leaf_area_cm2, 
         broad=ifelse(leafID>530,0,1),
         counte=count/area_samp) %>% 
  select(species, leafID, count, lesion, area_samp, broad, counte) %>% 
  droplevels() 

#plot data (estimated counts due to offsets)
df1 %>% 
  ggplot( aes(lesion, log(counte)))+
  geom_point(alpha=.4)+
  facet_wrap(~species, scales = 'free')

#MODEL 0 (no zero-inflation)
f0 <- bf(count ~ broad + lesion + (lesion|species) + offset(log(area_samp)), family = negbinomial())
get_prior(f0, df1)
prior0 <- c(
  set_prior("normal(2, 3.5)", class = "Intercept"), 
  set_prior("normal(0, 2)", class = "b", coef = 'lesion'),
  set_prior("normal(0, .5)", class = "b", coef = 'broad'),
  set_prior("exponential(1)", class = "sd"),
  set_prior("exponential(1)", class = "shape"))#dispersion
mnb0 <- brm(f0,prior0, data = df1, 
             family = negbinomial(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))
psims(df1, mnb0)+scale_x_continuous(limits=c(0,500))
pp_check(mnb0, type = "intervals_grouped", group = "species")#pretty bad fit for some species
pcoef(mnb0)#looks similar to the zinb

#MODEL 1
#priors
f1 <- bf(count ~ broad + lesion + (lesion|species) + offset(log(area_samp)), zi~1, family = zero_inflated_poisson())
get_prior(f1, df1)
prior1 <- c(
  set_prior("normal(2, 3.5)", class = "Intercept"), 
  set_prior("normal(0, 2)", class = "b", coef = 'lesion'),
  set_prior("normal(0, .5)", class = "b", coef = 'broad'),
  set_prior("exponential(1)", class = "sd"),
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi')) #flat prior
mzin1 <- brm(f1,prior1, data = df1, 
             family = zero_inflated_poisson(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))

#viz prob zero. flat 
rlogis(10000, 0,1) %>% inv_logit() %>%  density() %>% plot

#model fit looks like it has problems with adding too many zeros to some of the species. directly model zi~species?
psims(df1, mzin1)
pp_check(mzin1, type = "intervals_grouped", group = "species")#pretty bad fit for some species
fpred(df1, mzin1)+ scale_y_log10()#at least predictions better even if fit is awful
pcoef(mzin1)#whoa that's crazy


#MODEL 2
f2 <- bf(count ~ broad + lesion + (lesion|species) + offset(log(area_samp)), zi ~ species, family = zero_inflated_poisson())
get_prior(f2, df1)
prior2 <- c(
  set_prior("normal(2, 3.5)",class ="Intercept"), #species int
  set_prior("normal(0, 2)", class = "b", coef ='lesion'),
  set_prior("normal(0, .5)",class="b",coef='broad'), #assay
  set_prior("exponential(1)", class = "sd"),#all variation
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi'), #flat int prior
  set_prior("normal(0, .5)",class="b",dpar='zi'))#effect of species on zi probability
mzin2 <- brm(f2,prior2, data = df1, 
             family = zero_inflated_poisson(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))
summary(mzin2)
psims(df1, mzin2)#a little better. maybe try negbinom?
pp_check(mzin2, type = "intervals_grouped", group = "species")#also pretty bad
fpred(df1, mzin2)+scale_y_log10()
pcoef(mzin2)#whoa that's crazy

#MODEL 3 ZINB1 zi~sp
f3 <- bf(count ~ broad + lesion + (lesion|species) + offset(log(area_samp)), zi ~ species, family = zero_inflated_negbinomial())
get_prior(f3, df1)
prior2 <- c(
  set_prior("normal(2, 3.5)",class ="Intercept"),
  set_prior("normal(0, 2)", class = "b", coef ='lesion'),
  set_prior("normal(0, .5)",class="b",coef='broad'), #assay
  set_prior("exponential(1)", class = "sd"),#all variation
  set_prior("exponential(1)", class = "shape"),#dispersion
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi'),
  set_prior("normal(0, .5)",class="b",dpar='zi'))#effect of species on zi probability
mzinb3 <- brm(f3,prior2, data = df1, 
              family = zero_inflated_negbinomial(),
              chains = 4, cores = 4, iter = 2000, 
              control = list(adapt_delta = .95, max_treedepth=15))
summary(mzinb3, prob = .9)
psims(df1, mzinb3)
pp_check(mzinb3, type = "intervals_grouped", group = "species")#huge improvement
fpred(df1, mzinb3) + scale_y_log10()
pcoef(mzinb3)
slopes <- pcoef(mzinb3)

#MODEL 4 ZINB2 zi~1
f4 <- bf(count ~broad + lesion + (lesion|species) + offset(log(area_samp)), zi ~ 1, family = zero_inflated_negbinomial())
get_prior(f4, df1)
prior4 <- c(
  set_prior("normal(2, 3.5)",class ="Intercept"),
  set_prior("normal(0, 2)", class = "b", coef ='lesion'),
  set_prior("normal(0, .5)",class="b",coef='broad'), #assay
  set_prior("exponential(1)", class = "sd"),#all variation
  set_prior("exponential(1)", class = "shape"),#dispersion
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi')) #flat prior
mzinb4 <- brm(f4,prior4, data = df1, 
              family = zero_inflated_negbinomial(),
              chains = 4, cores = 4, iter = 2000, 
              control = list(adapt_delta = .95, max_treedepth=15))
summary(mzinb4, prob = .9)
psims(df1, mzinb4)
pp_check(mzinb4, type = "intervals_grouped", group = "species")
fpred(df1, mzinb4) + scale_y_log10()
pcoef(mzinb4)

#compare models with loo. zinb zi~species wins!
loo0 <- loo(mnb0)
loo1 <- loo(mzin1)
loo2 <- loo(mzin2)
loo3 <- loo(mzinb3)
loo4 <- loo(mzinb4)
loo_compare(loo0, loo1, loo2, loo3, loo4) #%>% write.csv(., 'output/lesion/loocomparison_chlamydo.csv', row.names = F)
#elpd_diff se_diff
#mzinb3     0.0       0.0
#mnb0     -19.9       7.9
#mzinb4   -21.5       5.5
#mzin2  -3609.3     972.0
#mzin1  -3626.7     965.7

#save models
saveRDS(mnb0, file = 'output/models/ch_lesion_nb0.rds')
saveRDS(mzin1, file = 'output/models/ch_lesion_zi1.rds')
saveRDS(mzin2, file = 'output/models/ch_lesion_zi2.rds')
saveRDS(mzinb3, file = 'output/models/ch_lesion_zinb3.rds')
saveRDS(mzinb4, file = 'output/models/ch_lesion_zinb4.rds')

mnb0 <- readRDS(file = 'output/models/ch_lesion_nb0.rds')
mzin1 <- readRDS(file = 'output/models/ch_lesion_zi1.rds')
mzin2 <- readRDS(file = 'output/models/ch_lesion_zi2.rds')
mzinb3 <- readRDS(file = 'output/models/ch_lesion_zinb3.rds')
mzinb4 <- readRDS(file = 'output/models/ch_lesion_zinb4.rds')

#saving outputs
sink('output/lesion/summary_lesionchl.txt')
summary(mzinb3, prob = .9)
sink()
pdf('plots/lesions/multiple_chl_mzinb3.pdf')
psims(df1, mzinb3)
pp_check(mzinb3, type = "intervals_grouped", group = "species")#huge improvement
fpred(df1, mzinb3) + scale_y_log10(limits=c(.001, 10000))
pcoef(mzinb3)
dev.off()
slopes <- pcoef(mzinb3)
slopes[[1]] %>% mutate_if(is.numeric, signif, 3) %>% 
  write_csv('output/lesion/randcoef_chlamydos.csv') 

#model coefficents into dataframe
vars <- get_variables(mzinb3)[1:14]
out <- mzinb3 %>% 
  gather_draws(!!!syms(vars)) %>% 
  median_hdi_sd() %>% 
  mutate_if(is.numeric, signif, 3) %>% 
  write_csv('output/lesion/chl_coef_estimates.csv')
