#chlamydos ~ species
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
write.csv <- function(...)write.csv(..., row.names = F)
inv_logit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x/(1-x))
#for easy reporting predicted values with ci and sd
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
psim <- function(data, model){
  msim <- data %>% 
    select(species, leafID, count, broad, area_samp) %>% 
    add_predicted_draws(mod) %>% 
    select(-.chain, -.iteration) %>% 
    group_by(species, .draw) %>% 
    sample_draws(50) 
  msim %>% ggplot() +
    geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.5, color='grey')+
    stat_density(data=data, aes(x=count), geom="line", color='slateblue')+
    facet_wrap(~species, scales = 'free')
}
#predicted species values
fpred <- function(data, model, observed=F){
  predict <- data %>%
    data_grid(species) %>% mutate(area_samp=1)%>%
    add_fitted_draws(model, scale = 'response') %>%
    mutate(assay='leaf disc')
  coef <- predict %>% median_hdi_sd()
  Plot <- ggplot(predict, aes(y=fct_rev(species), x = .value)) +
    geom_density_ridges(lwd=0, color=NA)+
    stat_pointintervalh(.width = c(.9), size=.2)+
    labs(x=expression(paste("Mean chlamydospores/", cm^{2})), y='Species')+
    scale_x_log10()
  
  #add observed data
  if(observed){
    Plot <- Plot + geom_point(data=data, aes(counte, species), fill='slateblue', color='slateblue', shape=24, alpha=.2) 
  }else{
    Plot
  }
  
  return(list(coef, Plot))
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
  ggplot( aes(species, log(counte+.01)))+
  geom_point(alpha=.4)

#subset the data
df1b <- df1 %>% filter(broad==1) %>% droplevels()
df1c <- df1 %>% filter(broad==0) %>% droplevels()

#MODEL 00 (pois no ZI)
f00 <- bf(count ~ -1 + species + offset(log(area_samp)), family = poisson())
get_prior(f00, df1b)
prior00 <- c(
  set_prior("normal(2, 3.5)", class = "b")) #species intercepts
mp00 <- brm(f00,prior00, data = df1b, 
            family = negbinomial(),
            chains = 4, cores = 4, iter = 2000,
            control = list(adapt_delta = .95, max_treedepth=15))
summary(mp00, prob = .9)
pp_check(mp00, type = "intervals_grouped", group = "species")#wow poisson is bad! underestimates variance a ton!
psim(df1b, mp00)
fpred(df1b, mp00, T)

#MODEL 0 (nb no ZI)
f0 <- bf(count ~ -1 + species + offset(log(area_samp)), family = negbinomial())
get_prior(f0, df1b)
prior0 <- c(
  set_prior("normal(2, 3.5)", class = "b"), #species intercepts
  set_prior("exponential(1)", class = "shape"))#dispersion
mnb0 <- brm(f0,prior0, data = df1b, 
            family = negbinomial(),
            chains = 4, cores = 4, iter = 2000,
            control = list(adapt_delta = .95, max_treedepth=15))
summary(mnb0, prob = .9)
pp_check(mnb0, type = "intervals_grouped", group = "species")
psim(df1b, mnb0)
fpred(df1b, mnb0, T)

#MODEL 1
#priors
f1 <- bf(count ~ -1 + species + offset(log(area_samp)), zi~1, family = zero_inflated_poisson())
get_prior(f1, df1b)
#viz what logistic curve looks like after backtransforming through the logit link function
rlogis(1000, 0, .5) %>% inv_logit() %>%  density() %>% plot
prior1 <- c(
  set_prior("normal(2, 3.5)", class = "b"), 
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi')) #flat prior
mzip1 <- brm(f1,prior1, data = df1b, 
             family = zero_inflated_poisson(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))
summary(mzip1, prob = .9)
pp_check(mzip1, type = "intervals_grouped", group = "species")#pretty bad
psim(df1b, mzip1)
fpred(df1b, mzip1, T)

#MODEL 2
f2 <- bf(count ~ -1 + species + offset(log(area_samp)), zi~species, family = zero_inflated_poisson())
get_prior(f2, df1b)
prior2 <- c(
  set_prior("normal(2, 3.5)", class = "b"), 
  set_prior("normal(2, 3.5)", class = "b", dpar = 'zi'), 
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi')) #flat prior
mzip2 <- brm(f2,prior2, data = df1b, 
             family = zero_inflated_poisson(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))
summary(mzip2, prob = .9)
pp_check(mzip2, type = "intervals_grouped", group = "species")
psim(df1b, mzip2)
fpred(df1b, mzip2, T)

#MODEL 3
f3 <- bf(count ~ -1 + species + offset(log(area_samp)), zi~1, family = zero_inflated_negbinomial())
get_prior(f3, df1b)
prior3 <- c(
  set_prior("normal(2, 3.5)", class = "b"), 
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi'),
  set_prior("exponential(1)", class = "shape"))
mzinb3 <- brm(f3,prior3, data = df1b, 
             family = zero_inflated_negbinomial(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))
summary(mzinb3, prob = .9)
pp_check(mzinb3, type = "intervals_grouped", group = "species")
psim(df1b, mzinb3)
fpred(df1b, mzinb3, T)

#MODEL 4
f4 <- bf(count ~ -1 + species + offset(log(area_samp)), zi~-1 + species, family = zero_inflated_negbinomial())
get_prior(f4, df1b)
prior4 <- c(
  set_prior("normal(2, 3.5)", class = "b"), 
  set_prior("logistic(0,1)", class = "b", dpar = 'zi'), 
  set_prior("exponential(1)", class = "shape"))
mzinb4 <- brm(f4,prior4, data = df1b, 
             family = zero_inflated_negbinomial(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))
summary(mzinb4, prob = .9)
pp_check(mzinb4, type = "intervals_grouped", group = "species")
psim(df1b, mzinb4)
fpred(df1b, mzinb4, T)

#compare all the models with LOO
loo00 <- loo(mp00)#18 obs pareto_k>.7
loo0 <- loo(mnb0)#1 obs pareto_k>.7
loo1 <- loo(mzip1)#22 obs pareto_k>.7
loo2 <- loo(mzip2)#25 obs pareto_k>.7
loo3 <- loo(mzinb3)#1 obs pareto_k>.7
loo4 <- loo(mzinb4)#1 obs pareto_k>.7
loo_compare(loo00, loo0, loo1, loo2, loo3, loo4)

#       elpd_diff se_diff
#mzinb4     0.0       0.0
#mzinb3   -39.3       7.3
#mnb0     -43.4      10.1
#mzip2  -6293.5    1191.9
#mzip1  -6324.5    1188.4
#mp00   -7328.6    1259.5

#save all the models
saveRDS(mp00, 'output/models/chlamydo_mp00.RDS')
saveRDS(mnb0, 'output/models/chlamydo_mnb0.RDS')
saveRDS(mzip1, 'output/models/chlamydo_mzip1.RDS')
saveRDS(mzip2, 'output/models/chlamydo_mzip2.RDS')
saveRDS(mzinb3, 'output/models/chlamydo_mzinb3.RDS')
saveRDS(mzinb4, 'output/models/chlamydo_mzinb4.RDS')
#read models
mp00 <- readRDS('output/models/chlamydo_mp00.RDS')
mnb0 <- readRDS('output/models/chlamydo_mnb0.RDS')
mzip1 <- readRDS('output/models/chlamydo_mzip1.RDS')
mzip2 <- readRDS('output/models/chlamydo_mzip2.RDS')
mzinb3 <- readRDS('output/models/chlamydo_mzinb3.RDS')
mzinb4 <- readRDS('output/models/chlamydo_mzinb4.RDS')

#get desired output
sink('output/chlamydo_both/summary_speciesbroad.txt')
summary(mzinb4, prob = .9)
sink()
pdf('plots/chlamydos_both/multiple_chl_broad.pdf')
pp_check(mzinb4, type = "intervals_grouped", group = "species")
psim(df1b, mzinb4)
fpred(df1b, mzinb4, T)
fpred(df1b, mzinb4, F)
dev.off()
coefs <-fpred(df1b, mzinb4)
coefs[[1]] %>% mutate_if(is.numeric, signif, 3) %>% 
  write.csv(., row.names = F, file = 'output/chlamydo_both/chlamydosbroad_predicted.csv') 

############################################################
#CONIFER MODELS! REPEAT PROCESS WITH OTHER DATASET.
############################################################
#MODEL 00 (pois no ZI)
f00 <- bf(count ~ -1 + species + offset(log(area_samp)), family = poisson())
get_prior(f00, df1c)
prior00 <- c(
  set_prior("normal(2, 3.5)", class = "b")) #species intercepts
mp00 <- brm(f00,prior00, data = df1c, 
            family = negbinomial(),
            chains = 4, cores = 4, iter = 2000,
            control = list(adapt_delta = .95, max_treedepth=15))

#MODEL 0 (nb no ZI)
f0 <- bf(count ~ -1 + species + offset(log(area_samp)), family = negbinomial())
get_prior(f0, df1c)
prior0 <- c(
  set_prior("normal(2, 3.5)", class = "b"), #species intercepts
  set_prior("exponential(1)", class = "shape"))#dispersion
mnb0 <- brm(f0,prior0, data = df1c, 
            family = negbinomial(),
            chains = 4, cores = 4, iter = 2000,
            control = list(adapt_delta = .95, max_treedepth=15))

#MODEL 1
#priors
f1 <- bf(count ~ -1 + species + offset(log(area_samp)), zi~1, family = zero_inflated_poisson())
get_prior(f1, df1c)
#viz what logistic curve looks like after backtransforming through the logit link function
rlogis(1000, 0, .5) %>% inv_logit() %>%  density() %>% plot
prior1 <- c(
  set_prior("normal(2, 3.5)", class = "b"), 
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi')) #flat prior
mzip1 <- brm(f1,prior1, data = df1c, 
             family = zero_inflated_poisson(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))

#MODEL 2
f2 <- bf(count ~ -1 + species + offset(log(area_samp)), zi~species, family = zero_inflated_poisson())
get_prior(f2, df1c)
prior2 <- c(
  set_prior("normal(2, 3.5)", class = "b"), 
  set_prior("normal(2, 3.5)", class = "b", dpar = 'zi'), 
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi')) #flat prior
mzip2 <- brm(f2,prior2, data = df1c, 
             family = zero_inflated_poisson(),
             chains = 4, cores = 4, iter = 2000,
             control = list(adapt_delta = .95, max_treedepth=15))

#MODEL 3
f3 <- bf(count ~ -1 + species + offset(log(area_samp)), zi~1, family = zero_inflated_negbinomial())
get_prior(f3, df1c)
prior3 <- c(
  set_prior("normal(2, 3.5)", class = "b"), 
  set_prior("logistic(0,1)", class = "Intercept", dpar='zi'),
  set_prior("exponential(1)", class = "shape"))
mzinb3 <- brm(f3,prior3, data = df1c, 
              family = zero_inflated_negbinomial(),
              chains = 4, cores = 4, iter = 2000,
              control = list(adapt_delta = .95, max_treedepth=15))

#MODEL 4
f4 <- bf(count ~ -1 + species + offset(log(area_samp)), zi~-1 + species, family = zero_inflated_negbinomial())
get_prior(f4, df1c)
prior4 <- c(
  set_prior("normal(2, 3.5)", class = "b"), 
  set_prior("logistic(0,1)", class = "b", dpar = 'zi'), 
  set_prior("exponential(1)", class = "shape"))
mzinb4 <- brm(f4,prior4, data = df1c, 
              family = zero_inflated_negbinomial(),
              chains = 4, cores = 4, iter = 2000,
              control = list(adapt_delta = .95, max_treedepth=15))

#compare all the models with LOO
loo00 <- loo(mp00)#18 obs pareto_k>.7
loo0 <- loo(mnb0)#1 obs pareto_k>.7
loo1 <- loo(mzip1)#22 obs pareto_k>.7
loo2 <- loo(mzip2)#25 obs pareto_k>.7
loo3 <- loo(mzinb3)#1 obs pareto_k>.7
loo4 <- loo(mzinb4)#1 obs pareto_k>.7
loo_compare(loo00, loo0, loo1, loo2, loo3, loo4)
loo_compare(loo00, loo0, loo1, loo2, loo3)

#differences between always between model_i and first model
#elpd_diff se_diff
#mzinb3    0.0       0.0 
#mzinb4   -0.1       0.8 
#mnb0     -2.7       1.8 
#mzip1  -235.7     126.7 
#mzip2  -241.8     130.7 
#mp00   -701.1     251.6 

#save models
saveRDS(mp00, 'output/models/chlamydo_mp00_conif.RDS')
saveRDS(mnb0, 'output/models/chlamydo_mnb0_conif.RDS')
saveRDS(mzip1, 'output/models/chlamydo_mzip1_conif.RDS')
saveRDS(mzip2, 'output/models/chlamydo_mzip2_conif.RDS')
saveRDS(mzinb3, 'output/models/chlamydo_mzinb3_conif.RDS')
saveRDS(mzinb4, 'output/models/chlamydo_mzinb4_conif.RDS')


#average the models? or just go with the best and simplest, mnb0
mod=mzinb3
pdf('plots/chlamydos_both/multiple_chl_conifer.pdf')
summary(mod, prob = .9)
pp_check(mod,type = "intervals_grouped", group = "species")
psim(df1c, mod)
fpred(df1c, mod, T)
dev.off()

##################################################
#final touches. merging the two models together and getting outputs and figures with both assays. 

#merge the sporulation predictions for both assays into 1 figure
modB <- readRDS('output/models/chlamydo_mzinb4.RDS')
modC <- readRDS('output/models/chlamydo_mzinb3_conif.RDS')
sink('output/chlamydo_both/summaries.txt')
summary(modB, prob =.9)
summary(modC, prob =.9)
sink()

#model coefficents into dataframe
vars <- get_variables(modB)[1:11]
outB <- modB %>% 
  gather_draws(!!!syms(vars)) %>% 
  median_hdi_sd() %>% 
  mutate_if(is.numeric, signif, 3) %>% 
  mutate(assay='broad')
vars <- get_variables(modC)[1:4]
outC <- modC %>% 
  gather_draws(!!!syms(vars)) %>% 
  median_hdi_sd() %>% 
  mutate_if(is.numeric, signif, 3) %>% 
  mutate(assay='conif')
bind_rows(outB, outC) %>% write.csv('output/chlamydo_both/coef_estimates.csv', row.names = F)

#predicted values into dataframe
predB <-fpred(df1b, modB)
predC <-fpred(df1c, modC)
pred <- bind_rows(predB[[1]], predC[[1]]) %>% 
  mutate_if(is.numeric, signif, 3) 
readr::write_csv(pred, 'output/chlamydo_both/chlamydosbroad_predicted_both.csv') 

#do pairwise contrasts
library(emmeans)
modB %>%
  emmeans( ~ species) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_hdci(.width=.9) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(signif=ifelse(.lower*.upper > 0, T, F)) %>% 
  bind_rows(
    modC %>%
      emmeans( ~ species) %>%
      contrast(method = "pairwise") %>%
      gather_emmeans_draws() %>%
      median_hdci(.width=.9) %>% 
      mutate_if(is.numeric, round, 2) %>% 
      mutate(signif=ifelse(.lower*.upper > 0, T, F))
  ) %>% 
  readr::write_csv('output/chlamydo_both/contrasts.csv') 

#make predictions
prB <- df1b %>%
  data_grid(species) %>% mutate(area_samp=1)%>%
  add_fitted_draws(modB, scale = 'response') %>% 
  mutate(assay='Leaf disc')
prC <- df1c %>%
  data_grid(species) %>% mutate(area_samp=1)%>%
  add_fitted_draws(modC, scale = 'response') %>%
  mutate(assay='Leaf dip')
pr <- bind_rows(prB, prC)

#change the name of LIDED
pr$species2 <- recode_factor(pr$species, LIDED="LIDE-D")

#plot 2 panel predictions
options(scipen=10000)#no scientific notation!
PLOT <- ggplot(pr, aes(y = fct_rev(species2), x = .value)) +
  geom_density_ridges(lwd=0, color=NA, panel_scaling = F, scale=1.35)+
  stat_pointintervalh(point_interval = median_hdi, .width = c(.9),shape=16, size=.1)+
  labs(x=expression(paste("Mean chlamydospores/", cm^{2})), y='Species')+
  #scale_x_log10()+
  scale_x_log10(limits=c(.1, 5000), breaks=c(.1, 1, 10, 100, 1000))+
  facet_grid(rows = vars(fct_rev(assay)), scales='free_y', space = 'free_y') 
PLOT
ggsave('plots/chlamydos_both/predicted_2panels.jpg', 
       PLOT, 
       width = 3.5, height = 2., units = 'in',
       dpi = 600)


