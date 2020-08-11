#redo spores vs. species with a nested structure. Run two modeles, considering species as a fixed effect or random effect.

rm(list = ls())#clear environment
library(tidyverse)
library(scales)
library(rethinking)
library(forcats)
library(brms)
library(tidybayes)
library(bayesplot)
library(loo)

#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
zscore <- function(x) (x-mean(x))/sd(x)
theme_set(theme_bw())#set ggplot theme


#load data
spores <- read.csv('data2019/master_tall.csv')

#clean data----

#clean up spore dataset
spores2 <- spores %>% 
  filter(spore_assay=="S", trt=="T") %>% 
  select(Species = species, ind, count1:count3, leaf_area_cm2) %>% 
  pivot_longer(cols = count1:count3, names_prefix = 'count',  names_to = 'rep', values_to = 'count') %>% 
  droplevels() %>% 
  mutate(rep = as.integer(rep),
         leafID = interaction(Species, ind)) %>% 
  filter(!is.na(count), !is.na(leaf_area_cm2)) 
print(spores2, width = Inf)


#how many are zeros?
dens(spores2$count)
sum(spores2$count==0)/length(spores2$count) #44% are zeros.



#model it----

#ignore forest types and look at both n.plots and rank
#make most counts under 100
N=10000
a0 <- rnorm(N, 2.2, 1)

xrange <- seq(-2, 2, length.out = 100)
dens(exp(a0), xlim=c(0,300), show.HPDI = .9) + abline(v=c(50,100))
mean(exp(a0))




#using a nested model----


#do it in Stan...

#prep datasets----

disc_spp <- c('ACMA', 'ARME', 'CEOL', 'HEAR', 'LIDE', 'QUAG', 'QUCH', 'QUPA', 'TODI', 'UMCA', 'CONTROL')
dip_spp <- c('PIPO', 'SESE', 'PSME', 'UMCAD', 'LIDED')

#leaf disc assay
wholetab_disc <- spores2 %>% 
  filter(Species %in% disc_spp) %>% 
  droplevels() %>% 
  mutate(SpID = as.integer(as.factor(Species)),
         leafID = as.integer(as.factor(as.numeric(leafID)))) 
Xi <- wholetab_disc %>% select(leaf_area_cm2, count, leafID)
Xj <- wholetab_disc %>% group_by(leafID) %>% summarise(SpID = unique(SpID)) 
df_list_disc <- list(N = nrow(Xi),
                J = nrow(Xj),
                S = max(Xj$SpID),
                count = Xi$count,
                leafID = Xi$leafID,
                leaf_area_cm2 = Xi$leaf_area_cm2,
                spID = Xj$SpID
)
str(df_list_disc)


#leaf dip assay
wholetab_dip <- spores2 %>% 
  filter(Species %in% dip_spp) %>% 
  droplevels() %>% 
  mutate(SpID = as.integer(as.factor(Species)),
         leafID = as.integer(as.factor(as.numeric(leafID)))) 
Xi <- wholetab_dip %>% select(leaf_area_cm2, count, leafID)
Xj <- wholetab_dip %>% group_by(leafID) %>% summarise(SpID = unique(SpID)) 
df_list_dip <- list(N = nrow(Xi),
                     J = nrow(Xj),
                     S = max(Xj$SpID),
                     count = Xi$count,
                     leafID = Xi$leafID,
                     leaf_area_cm2 = Xi$leaf_area_cm2,
                     spID = Xj$SpID
)
str(df_list_dip)




#model it!----

#species as random effect or fixed effect
mod_spRE <- stan_model(file = 'rscripts/analysis2019/sporangia_nested_spRE.stan')
mod_spFE <- stan_model(file = 'rscripts/analysis2019/sporangia_nested_spFE.stan')

#leaf disc assay
fit_spRE <- sampling(mod_spRE, data = df_list_disc, iter = 2000, chains = 4, cores = 4)
fit_spFE <- sampling(mod_spFE, data = df_list_disc, iter = 2000, chains = 4, cores = 4)
precis(fit_spRE, depth = 2, pars = c('a0', 'sdID', 'sdSp', 'zSp')) 
precis(fit_spFE, depth = 2, pars = c('a0', 'sdID', 'aSp')) 
wholetab_disc %>% distinct(Species, SpID) %>% arrange(SpID)


#leaf dip assay
fit_spRE_dip <- sampling(mod_spRE, data = df_list_dip, iter = 2000, chains = 4, cores = 4)
precis(fit_spRE_dip, depth = 2, pars = c('a0', 'sdID', 'sdSp', 'zSp')) 
wholetab_dip %>% distinct(Species, SpID) %>% arrange(SpID)



#quickly evaluate models
post <- extract.samples(fit_spFE)
ppc_dens_overlay(df_list_disc$count, post$y_rep[1:50,])
loofit <- loo(fit_spFE)
loofit
post <- extract.samples(fit_spRE_dip)
ppc_dens_overlay(df_list_dip$count, post$y_rep[1:50,])


#test difference between individual-level and species-level variation. not significant statistically speaking.
postRE <- extract.samples(fit_spRE)
diffSD <- postRE$sdID - postRE$sdSp
dens(diffSD, show.HPDI = .9)
median(diffSD); HPDI(diffSD, .9)



#test difference of individual-level between assays. is the dip variation greater? Yes, signficantly so.
postREdip <- extract.samples(fit_spRE_dip)
par(mfrow=c(2,1))
postREdip$sdID %>% dens(xlim=c(0,2))
postRE$sdID %>% dens(xlim=c(0,2))
par(mfrow=c(1,1))
diffsdID <- postREdip$sdID - postRE$sdID 
dens(diffsdID, show.HPDI = .9)
median(diffsdID); HPDI(diffsdID, .9)
