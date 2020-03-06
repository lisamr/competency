#Presentation for NUDE

#sporangia analysis of recent inoculation tests. I (Lisa) inoculated common CA coast species with P. ramorum, the pathogen that causes sudden oak death, so I can know which species are best at transmitting the pathogen. This data will help me learn how community compostion affects disease risk using a forest disease system.

#SETUP----
setwd("~/Box/Competency project/competency.git/NUDE/")#set your workspace to whereever the files are.
rm(list = ls())#clear environment
library(tidyverse) #includes dplyr, ggplot, forcats
library(ggridges)
library(tidybayes)
library(brms)

#set theme
theme_set(theme_bw(base_size = 9) + theme(panel.grid=element_blank())) 
#IGNORE CHUNK. I just filtered out the data I'll be presenting on.
#df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just trt, sporangia, leaf disc assay. Can give the workshop this data.
#df1 <- df %>% filter(trt=="T", spore_assay=="S", leafID<=530|species=="CONTROL", is.na(omit)) %>% mutate(species = droplevels(species))
#write.csv(df1, 'NUDE/master_tall_NUDE.csv', row.names = F)

#READ DATA----
#read in master file
df1 <- read.csv('master_tall_NUDE.csv') #fixed typo in directory.

#check out data
head(df1)
summary(df1)

#make data tall in regards to counts
dtall <- df1 %>% 
  reshape2::melt(df1, id.vars = c("species", "leafID", 'leaf_area_cm2'), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count), !is.na(leaf_area_cm2)) %>% 
  arrange(leafID) 

#viz sporangia data. just look at means
dtall %>% 
  group_by(species, leafID) %>% 
  summarise(spores = mean(count)) %>% 
  ggplot(., aes(species, spores)) +
  geom_boxplot() +
  geom_point() +
  scale_y_log10()

#QUESTION AND MODEL----
################################################
#Question: how does sporulation vary among host species?
################################################

#model sporangia counts as a function of species. Since I have 3 counts per sample, I need leafID as a varying intercept (random effect).

#spores ~ Poisson(lambda) #likelihood of counts
#log(lambda) = alpha_species[i] + alpha_leafID[i] #linear regression 
# alpha_species ~ Normal(1.5, 1.5) #see priors below
# alpha_leafID ~ Normal(0, sigma) #random effect
# sigma ~ exp(1) #variation around the mean due to sampling

#figure out priors
#species int. 0-1000 sp/cm2? based on larch study. larch was above (~2000) but other species were way below. divide that by 102.5: 0-10 mostly

#same prior as broadleaf model?
#put the distribution on the outcome scale to get an intuitive sense of what you're telling the model. plot it!
rnorm(10000, 7, 1) %>% exp %>%  density() %>% plot(xlim=c(0,10000))

#model formula in brms
f1 <- bf(count ~ -1 + species + (1|leafID), family = poisson())

#check out default priors
get_prior(f1, dtall)

#set your own. I think I may have done this right, possibly it's wrong? anyone out there to pipe in?
prior1 <- c(
  set_prior("normal(1.5, 1.5)", class = "b"),
  set_prior("exponential(1)", class = "sd") )

#model finally
#broadleaf model. had to increase the n.iterations because Bulk Effective Samples Size was too low.
#commenting out since it takes a while to run.
m1 <- brm(formula = f1, data = dtall, family = poisson(), chains = 4, cores = 4, prior = prior1, iter=4000, control = list(adapt_delta = .98, max_treedepth=15))

#save models
saveRDS(m1, 'sporangia_NUDE.RDS')
m1 <- readRDS(file = 'sporangia_NUDE.RDS')

#INITIAL LOOKS AT MODEL----
#check out coefs and save summaries
summary(m1, prob=.9)

#population level effects only. traceplots look good in both models. good diagnostic tool to make sure you're model is running correctly.
plot(m1, pars = "^b_")

#see model fit. seems good.
pp_check(m1) #all data together. x-axis is #spores
pp_check(m1,type = "intervals_grouped", group = "species")#grouped by species. looks great.

#simulate prediction data using same dataset to see fit. takes a little time.
sim_m1 <- add_predicted_draws(dtall, m1) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(50) 

#plot against your own data. Another way to check out model fit.
sim_m1 %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.2, color='grey')+
  stat_density(data=dtall, aes(x=count), geom="line", color='steelblue') +
  facet_wrap(~species, scales = 'free')

#CHECKING OUT POSTERIOR----
#coef plot to see what the model thinks about your species.
#draw parameter values from the posterior.
coefs_m1 <- m1 %>% 
  gather_draws(b_speciesCONTROL, b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value)

#print a summary table of your coefs or just plot it
coefs_m1 %>% median_hdci(.width = c(.5, .9)) %>% 
  ggplot(., aes(value, par)) +
  geom_pointintervalh(size_range = c(.1,.5))

#I want to know what these numbers actually mean on the outcome scale--how many spores does it predict on average for each species.

#PREDICTIONS----
#get predicted fit based on just the species coef. Ignoring random effects cuz I don't care about predicting sampling error.
post_pred <- dtall %>%
  modelr::data_grid(species) %>% #say what data to predict with. here it's just the species names.
  add_fitted_draws(m1, re_formula = ~0, scale = 'response') 
print(post_pred, width=Inf) #check it out. similar to the draws, but posterior predictions out outcome scale rather than model parameters.

#get the summarized values
median_hdi(post_pred, .value) %>% print(width=Inf)

#plot predictions showing the distributions.
#halfeye function. geom_density_ridges plots fake data. use halfeye from tidybayes instead.
post_pred %>% 
  ggplot(., aes(y = fct_rev(species), x = .value)) +
  geom_halfeyeh(.width = .9, size=.1, fatten_point = .1, relative_scale = 2, point_interval = median_hdi) +
  labs( x=expression(paste("Mean sporangia/", cm^{2})), 
        y='Species') +
  scale_x_continuous(limits=c(0, 15))

#CONTRASTS----
#calculate contrasts.
#never figured out how to do it the tidy way, so doing it kinda clunky

#extract posterior of following coefficients
get_variables(m1)
draws <- m1 %>% 
  spread_draws(b_speciesCONTROL, b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA) %>% 
  select(-c(.chain, .iteration, .draw)) %>% as.data.frame()

#get pairwise comparisons
pairs <- combn(1:ncol(draws),2) #get pairwise combination for all 11 species

#calculates pairwise differences between each posterior combo. outputs HDCI of difference.
f <- function(draws, pairs, x){ 
  spdiff <-  (draws[,pairs[1,x]]-draws[,pairs[2,x]])
  hdci(spdiff, .9)
}
constrasts <- sapply(1:ncol(pairs),function(x) f(draws, pairs, x))

#put contrasts into a readable dataframe
spp <- colnames(draws) %>% str_replace('b_species', '')
constrasts_df <- rbind(pairs, constrasts) %>% t %>% as.data.frame
constrasts_df <- constrasts_df %>% 
  mutate(sp1=spp[constrasts_df[,1]],
         sp2=spp[constrasts_df[,2]],
         sig=sign(constrasts_df[,3])==sign(constrasts_df[,4]))
head(constrasts_df)
