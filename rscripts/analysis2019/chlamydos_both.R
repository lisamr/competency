#chlamydospores estimations models for both conifers and broadleaves

#starting new script to analyze chlamydos. using a multilevel model with a observation level random effect to account for the dispersion and a poisson likelihood. Im hoping i can ignore the zero inflation but i dont think i can.
setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(tidybayes)
library(modelr)
library(brms)
library(ggridges)
library(ggstance)
library(eemeans)

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
  select(-count2, -count3) %>% 
  rename(count=count1) %>% 
  mutate(area_samp=prop*leaf_area_cm2) %>% 
  droplevels() 

#plot data (estimated counts due to offsets)
df1 %>% mutate(count_est=count/area_samp) %>% 
  ggplot( aes(species, log(count_est)))+
  geom_boxplot()+
  geom_point(alpha=.4)


################################################
#model it!
#still think I should run seperate models. But I'll see what happens when I throw all the species into one. Should I include a fixed effect for that?

#remove hear and pipo because they seems to be causing problems and we know it's zero.
df2 <- filter(df1, !species %in% c("HEAR", "PIPO")) %>% 
  select(species, ind, leafID, count, area_samp, prop, leaf_area_cm2)
#df2 <- filter(df1, !is.na(count|area_samp)) %>% select(species, ind, leafID, count, area_samp, prop, leaf_area_cm2)
df2$species <- droplevels(df2$species)
df2$broad <- ifelse(df2$leafID>530, 0, 1)
#subset data for seperate models. 
df3b <- filter(df2, leafID<=530) %>% droplevels()
df3c <- filter(df2, leafID>530) %>% droplevels()

#check out default priors
f1 <- bf(count ~ -1 + species + (1|leafID) + offset(log(area_samp)), family = poisson())
f2 <- bf(count ~ -1 + species + broad + (1|leafID) + offset(log(area_samp)), family = poisson())
get_prior(f2, df2)

#set your own
rnorm(1000, 2, 3.5) %>% dens
prior1 <- c(
  set_prior("normal(2, 3.5)", class = "b"), #species int 
  set_prior("exponential(1)", class = "sd")) #all the variations
prior2 <- c(
  set_prior("normal(2, 3.5)", class = "b"), #species int 
  set_prior("normal(0, .5)", class = "b", coef = 'broad'), #assay
  set_prior("exponential(1)", class = "sd")) #all the variations

#run model
#no warnings when i remove hear and pipo
#no assay coeffient
m1 <- brm(formula = f1, prior = prior1,
          data = df2, family = poisson(),
          chains = 4, cores = 4, iter = 4000,
          control = list(adapt_delta = .95, max_treedepth=15))
m2 <- brm(formula = f2, prior = prior2,
          data = df2, family = poisson(),
          chains = 4, cores = 4, iter = 4000,
          control = list(adapt_delta = .95, max_treedepth=15))
#seperate models
m3b <- brm(formula = f1, prior = prior1,
          data = df3b, family = poisson(),
          chains = 4, cores = 4, iter = 4000,
          control = list(adapt_delta = .95, max_treedepth=15))
m3c <- brm(formula = f1, prior = prior1,
           data = df3c, family = poisson(),
           chains = 4, cores = 4, iter = 4000,
           control = list(adapt_delta = .95, max_treedepth=15))

#compare models m1 and m2
#m1 generally has slightly smaller intervals
drawsm1 <- m1 %>% 
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesLIDE, b_speciesLIDED, b_speciesPSME, b_speciesTODI, sd_leafID__Intercept) %>% 
  median_qi(estimate = .value) %>%
  mutate(model="m1")
drawsm2 <- m2 %>% 
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesLIDE, b_speciesLIDED, b_speciesPSME, b_speciesTODI, sd_leafID__Intercept) %>% 
  median_qi(estimate = .value) %>%
  mutate(model="m2")
drawsm3b <- m3b %>% 
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesLIDE, b_speciesTODI, sd_leafID__Intercept) %>% 
  median_qi(estimate = .value) %>%
  mutate(model="m3b")
drawsm3c <- m3c %>% 
  gather_draws(b_speciesLIDED, b_speciesPSME, sd_leafID__Intercept) %>% 
  median_qi(estimate = .value) %>%
  mutate(model="m3c")

bind_rows(drawsm1, drawsm2, drawsm3b, drawsm3c) %>% 
  ggplot(aes(y = .variable, x = estimate, color=model)) +
  geom_pointintervalh(position = position_dodgev(.3)) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient') %>% 
  scale_color_viridis_d(end = .9)

#######################################################
#check out models individually
#check out coefs
mod=m3b
dat=df3b
mod=m3c
dat=df3c

summary(mod)
fixef(mod)
#population level effects only. traceplots look ok
plot(mod, pars = "^b_")#looks good
post <- posterior_samples(mod, pars = c("^b_", "sd"))
pp_check(mod)

#see model fit. seems pretty good.
#simulate prediction data using same dataset to see fit
msim <- dat %>% 
  select(species, leafID, count, broad, area_samp) %>% 
  add_predicted_draws(mod) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(50) 
msim %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.5, color='grey')+
  stat_density(data=dat, aes(x=count), geom="line", color='slateblue')+
  facet_wrap(~species, scales = 'free')

#coefplot
parnames(mod)
mod %>% 
  #gather_draws(b_speciesLIDED, b_speciesPSME,  sd_leafID__Intercept) %>%
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesLIDE, b_speciesTODI, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')

########################################
#need to plot estimates backtransformed to spores/cm2. create a plot with 2 panels like the one with sporangia.
#predicted means
prm3b <- df3b %>%
  data_grid(species) %>% mutate(area_samp=1)%>%
  add_fitted_draws(m3b, re_formula = ~0, scale = 'response') %>% 
  mutate(assay='leaf disc')
prm3c <- df3c %>%
  data_grid(species) %>% mutate(area_samp=1)%>%
  add_fitted_draws(m3c, re_formula = ~0, scale = 'response') %>% 
  mutate(assay='detached leaf')
pr <- bind_rows(prm3b, prm3c)
#pr <- pr %>% ungroup() %>% add_row(species="HEAR", .value=0.009, assay="leaf disc") %>% add_row(species="HEAR", .value=0.009, assay="leaf disc") %>% add_row(species="PIPO", .value=0.009, assay="detached leaf") %>%  add_row(species="PIPO", .value=0.009, assay="detached leaf")
prm3b %>% 
  median_hdi_sd(.value) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  bind_rows(
    prm3c %>% 
      median_hdi_sd(.value) %>% 
      mutate_if(is.numeric, round, 2)
  ) %>% 
  write_csv('output/chlamydo_both/estimates.csv')

#plot with observed data
obs <- df2 %>% mutate(counte=count/area_samp)
ggplot(pr, aes(y = fct_rev(species), x = .value)) +
  geom_density_ridges(lwd=0, color=NA)+
  stat_pointintervalh(.width = c(.9), size=.2)+
  labs(x=expression(paste("Mean chlamydospores/", cm^{2})), y='Species')+
  scale_x_log10()+
  geom_point(data=obs, aes(counte, species), fill='slateblue', color='slateblue', shape=24, alpha=.2) 

pdf('plots/chlamydos_both/predicted_2panels_HEARPIPO.pdf', 3.25, 2.75)
ggplot(pr, aes(y = fct_rev(species), x = .value)) +
  geom_density_ridges(lwd=0, color=NA, panel_scaling = F)+
  stat_pointintervalh(point_interval = median_hdi, .width = c(.9),shape=21, size=.1, fill='white')+
  labs(x=expression(paste("Mean chlamydospores/", cm^{2})), y='Species')+
  #scale_x_log10()+
  scale_x_log10(limits=c(.009, 5000), breaks=c(.1, 1, 10, 100, 1000))+
  facet_grid(rows = vars(fct_rev(assay)), scales='free_y', space = 'free_y') +
  theme(text = element_text(size=8), panel.grid = element_blank())
dev.off()

#do pairwise contrasts
m3b %>%
  emmeans( ~ species) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_hdci(.width=.9) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(signif=ifelse(.lower*.upper > 0, T, F)) %>% 
  bind_rows(
    m3c %>%
      emmeans( ~ species) %>%
      contrast(method = "pairwise") %>%
      gather_emmeans_draws() %>%
      median_hdci(.width=.9) %>% 
      mutate_if(is.numeric, round, 2) %>% 
      mutate(signif=ifelse(.lower*.upper > 0, T, F))
  ) %>% 
  write_csv('output/chlamydo_both/contrasts.csv')

#print summary of models
sink('output/chlamydo_both/summary_modelbroad.txt')
summary(m3b)
sink()
sink('output/chlamydo_both/summary_modelconifer.txt')
summary(m3c)
sink()


