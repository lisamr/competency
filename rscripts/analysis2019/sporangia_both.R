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

#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
theme_set(theme_light(base_size = 12)) #set ggplot theme

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just trt, sporangia
df1 <- df %>% 
  filter(trt=="T", spore_assay=="S") 

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
  ggplot( aes(species, (countm/leaf_area_cm2)))+
  geom_boxplot()+
  geom_point(alpha=.4)+
  scale_y_log10()

df1[is.na(df1$countm),] %>% View#the samples without counts were too low in vol or too goopy


################################################
################################################
#model sporangia counts. include leaf area as an offset or standardize the posteriors after modeling? I think it should be an offset.

#make data tall in regards to counts
dtall <- melt(df1, id.vars = c("species", "leafID", 'leaf_area_cm2'), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count), !is.na(leaf_area_cm2)) %>% 
  arrange(leafID) 
#split data by assay--broadleaf and conifer. modeling them seperately.
dtallb <- dtall %>% filter(leafID<=530) %>% droplevels()
dtallc <- dtall %>% filter(leafID>530) %>% droplevels()

#figure out priors
#species int. 0-1000 sp/cm2? based on larch study. larch was above (~2000) but other species were way below. divide that by 102.5: 0-10 mostly
#same prior as broadleaf model?
dens(exp(rnorm(1000, 2, 1.5)), xlim=c(0,100)) 
mean(exp(rnorm(1000, 2, 1.5)))
#model formula. offset in there for conifer model and doesn't change anything if its left in for the broadleaf model.
f1 <- bf(count ~ -1 + species + (1|leafID) + offset(log(leaf_area_cm2)), family = poisson())

#check out default priors
get_prior(f1, dtall)

#set your own
prior1 <- c(
  set_prior("normal(2, 1.5)", class = "b"),
  set_prior("exponential(1)", class = "sd") )

#model finally
#broadleaf model. had to increase the n.iterations because Bulk Effective Samples Size was too low.
mb <- brm(
  formula = f1, data = dtallb, family = poisson(),
  chains = 4, cores = 4, prior = prior1, iter=4000,
  control = list(adapt_delta = .95, max_treedepth=15))
#conifer model
mc <- brm(
  formula = f1, data = dtallc, family = poisson(),
  chains = 4, cores = 4, prior = prior1,
  control = list(adapt_delta = .95, max_treedepth=15))
#no warnings :)

#check out coefs and save summaries
#sink('output/sporangia_both/summary_modelbroad.txt')
summary(mb)
#sink()
#sink('output/sporangia_both/summary_modelconifer.txt')
summary(mc)
#sink()

#population level effects only. traceplots look good in both models.
plot(mb, pars = "^b_")
plot(mc, pars = "^b_")
#see model fit. seems good.
pp_check(mb)
pp_check(mc)

#simulate prediction data using same dataset to see fit
sim_mb <- add_predicted_draws(dtallb, mb) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(50) 
sim_mc <- add_predicted_draws(dtallc, mc) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(50) 
#bind predictions together before plotting fit
sims <- bind_rows(sim_mb, sim_mc)
#set levels order of species
sims$species <- as.factor(sims$species)
neworder <- c(1:5, 9:11, 13, 14, 6, 15, 7, 8, 12)
levels(sims$species) <- levels(sims$species)[neworder]
levels(dtall$species) <- levels(dtall$species)[neworder]

#plot against your own data
#pdf('plots/sporangia_both/modelfit_brms.pdf', 10, 8)
sims %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.2, color='grey')+
  stat_density(data=dtall, aes(x=count), geom="line", color='steelblue') +
  facet_wrap(~species, scales = 'free')
#dev.off()  

#coef plot. keep models seperate.
coefs_mb <- mb %>% 
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  mutate(par2 = if_else(
    par=="sd_leafID__Intercept", "sd_leafID_broad", par))
coefs_mc <- mc %>% 
  gather_draws(b_speciesLIDED, b_speciesPIPO, b_speciesPSME, b_speciesSESE, b_speciesUMCAD, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  mutate(par2 = if_else(
    par=="sd_leafID__Intercept", "sd_leafID_conifer", par))

#plots
coefs_mb %>% 
  ggplot(aes(y = fct_rev(par2), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')
coefs_mc %>% 
  ggplot(aes(y = fct_rev(par2), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')

#get predicted fit based on just the species coef.
post_mb <- dtallb %>%
  modelr::data_grid(species, leaf_area_cm2=1) %>%
  add_fitted_draws(mb, re_formula = ~0, scale = 'response') %>%
  mutate(.valueSTD = .value*102.5, assay='Leaf disc') #standardized sp/cm2
post_mc <- dtallc %>%
  modelr::data_grid(species, leaf_area_cm2=1) %>%
  add_fitted_draws(mc, re_formula = ~0, scale = 'response') %>%
  mutate(.valueSTD = .value*102.5, assay='Detached leaf') #standardized sp/cm2

#plot them together
head(post_mb)
post <- bind_rows(post_mb, post_mc)
#fix species factor order
post$species <- factor(post$species, levels = unique(post$species)[c(1:11, 15, 12:14)])
#plot it!
#pdf('plots/sporangia_both/meansporangia_brms.pdf', 10, 8)
ggplot(post, aes(y = fct_rev(species), x = .valueSTD)) +
  geom_density_ridges(lwd=.1, col=alpha('black', .5))+
  stat_pointintervalh(point_interval = median_hdcih, .width = c(.9), size=.2)+
  labs( x=expression(paste("Mean sporangia/", cm^{2})), 
    y='Species') +
  facet_grid(rows = vars(fct_rev(assay)), scales='free_y', space = 'free_y') + 
  scale_x_continuous(limits=c(-25, 1250), breaks = seq(0,1250, 250))
#dev.off()

#get the summarized values
postmb <- post_mb %>% 
  median_hdci(x=.valueSTD, .width = c(.9, .95)) %>% 
  mutate_if(is.numeric, signif, 4)
postmc <- post_mc %>% 
  median_hdci(x=.valueSTD, .width = c(.9, .95))%>% 
  mutate_if(is.numeric, signif, 4)
#write_csv(postmb, 'output/sporangia_both/meanspeciesint_broadleaf.csv')
#write_csv(postmc, 'output/sporangia_both/meanspeciesint_conifer.csv')

#calculate contrasts. keep the assays seperate.
#can't figure out how to do it the tidy way, so doing it kinda clunky.

#extract posterior of following coefficients
exb <- mb %>% 
  spread_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA) %>% 
  select(-c(.chain, .iteration, .draw)) %>% as.data.frame()
exc <- mc %>% 
  spread_draws(b_speciesLIDED, b_speciesPIPO, b_speciesPSME, b_speciesSESE, b_speciesUMCAD) %>% 
  select(-c(.chain, .iteration, .draw)) %>% as.data.frame()

#get pairwise comparisons
pairs_b <- combn(1:10,2)
pairs_c <- combn(1:5,2)
f <- function(ex, pairs, x){
  spdiff <-  (ex[,pairs[1,x]]-ex[,pairs[2,x]])
  hdci(spdiff, .9)
}
diffs_b <- sapply(1:ncol(pairs_b),function(x) f(exb, pairs_b, x))
diffs_c <- sapply(1:ncol(pairs_c),function(x) f(exc, pairs_c, x))

#put contrasts into a readable dataframe
sppb <- colnames(exb) %>% str_replace('b_species', '')
sppc <- colnames(exc) %>% str_replace('b_species', '')
diffsdf_b <- rbind(pairs_b, diffs_b) %>% t %>% as.data.frame
diffsdf_c <- rbind(pairs_c, diffs_c) %>% t %>% as.data.frame
diffsdf_b <- diffsdf_b %>% 
  mutate(sp1=sppb[diffsdf_b[,1]],
         sp2=sppb[diffsdf_b[,2]],
         sig=sign(diffsdf_b[,3])==sign(diffsdf_b[,4]))
diffsdf_c <- diffsdf_c %>% 
  mutate(sp1=sppc[diffsdf_c[,1]],
         sp2=sppc[diffsdf_c[,2]],
         sig=sign(diffsdf_c[,3])==sign(diffsdf_c[,4]))
#write_csv(diffsdf_b, 'output/sporangia_both/contrasts_broad.csv') 
#write_csv(diffsdf_c, 'output/sporangia_both/contrasts_conifer.csv') 
