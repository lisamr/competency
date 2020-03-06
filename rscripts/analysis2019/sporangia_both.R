#sporangia analysis for both assays--broadleaf and conifers. pull the posteriors from the two models, standardize the units and plot the results together.

setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(tidyverse) #includes dplyr, ggplot, forcats
library(reshape2)
library(scales)
library(ggridges)
library(tidybayes)
library(brms)
library(ggthemes)

#functions
reset <- function(x) par(mfrow=c(1,1))#turning pars back 
theme_set(theme_bw(base_size = 9) + theme(panel.grid=element_blank())) #set ggplot theme
dens <- function(x, ...)plot(density(x), ...)
median_hdi_sd <- function(data, value=.value, width=.9){
    value <- enquo(value)
    data %>% 
      summarise(estimate=median(!!value), 
                lower=hdci(!!value, width)[1], 
                upper=hdci(!!value, width)[2],
                sd=sd(!!value),
                width=.9)
}

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
dtall <- df1 %>% 
  melt(df1, id.vars = c("species", "leafID", 'leaf_area_cm2'), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count), !is.na(leaf_area_cm2)) %>% 
  arrange(leafID) 
#split data by assay--broadleaf and conifer. modeling them seperately.
dtallb <- dtall %>% filter(leafID<=530 | species=="CONTROL") %>% droplevels()
dtallc <- dtall %>% filter(leafID>530, species!="CONTROL") %>% droplevels()

#figure out priors
#species int. 0-1000 sp/cm2? based on larch study. larch was above (~2000) but other species were way below. divide that by 102.5: 0-10 mostly
#same prior as broadleaf model?
rnorm(10000, 1.5, 1.5) %>% exp %>% dens(xlim=c(0,40))
rnorm(10000, 1.5, 1.5) %>% dens()
exp(3)*102.5
20*102.5
#model formula. offset in there for conifer model and doesn't change anything if its left in for the broadleaf model.
f1 <- bf(count ~ -1 + species + (1|leafID) + offset(log(leaf_area_cm2)), family = poisson())

#check out default priors
get_prior(f1, dtall)

#set your own
prior1 <- c(
  set_prior("normal(1.5, 1.5)", class = "b"),
  set_prior("exponential(1)", class = "sd") )

#model finally
#broadleaf model. had to increase the n.iterations because Bulk Effective Samples Size was too low.
mb <- brm(
  formula = f1, data = dtallb, family = poisson(),
  chains = 4, cores = 4, prior = prior1, iter=4000,
  control = list(adapt_delta = .98, max_treedepth=15))
#conifer model
mc <- brm(
  formula = f1, data = dtallc, family = poisson(),
  chains = 4, cores = 4, prior = prior1,
  control = list(adapt_delta = .95, max_treedepth=15))
#no warnings :)

#save models
saveRDS(mb, 'output/models/sporangia_broad.RDS')
saveRDS(mc, 'output/models/sporangia_conif.RDS')
mb <- readRDS('output/models/sporangia_broad.RDS')
mc <- readRDS('output/models/sporangia_conif.RDS')

#check out coefs and save summaries
sink('output/sporangia_both/summaries.txt')
summary(mb, prob =.9)
summary(mc, prob =.9)
sink()

#population level effects only. traceplots look good in both models.
plot(mb, pars = "^b_")
plot(mc, pars = "^b_")
#see model fit. seems good.
pp_check(mb)
pp_check(mc)
pp_check(mb,type = "intervals_grouped", group = "species")
pp_check(mc,type = "intervals_grouped", group = "species")

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

#plot against your own data
pdf('plots/sporangia_both/modelfit_brms.pdf', 10, 8)
sims %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.2, color='grey')+
  stat_density(data=dtall, aes(x=count), geom="line", color='steelblue') +
  facet_wrap(~species, scales = 'free')
dev.off()  

#coef plot. keep models seperate.
coefs_mb <- mb %>% 
  gather_draws(b_speciesCONTROL, b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  mutate(assay='broad')
coefs_mc <- mc %>% 
  gather_draws(b_speciesLIDED, b_speciesPIPO, b_speciesPSME, b_speciesSESE, b_speciesUMCAD, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  mutate(assay='conif')
bind_rows(coefs_mb, coefs_mc) %>% 
  group_by(par, assay) %>% 
  median_hdi_sd(value = value, width = .9) %>% 
  arrange(assay) %>% 
  mutate_if(is.numeric, signif, 3) %>% 
  write_csv('output/sporangia_both/coefficients.csv')

#get predicted fit based on just the species coef.
post_mb <- dtallb %>%
  modelr::data_grid(species, leaf_area_cm2=1) %>%
  add_fitted_draws(mb, re_formula = ~0, scale = 'response') %>%
  mutate(.valueSTD = .value*102.5, assay='Leaf disc')
post_mc <- dtallc %>%
  modelr::data_grid(species, leaf_area_cm2=1) %>%
  add_fitted_draws(mc, re_formula = ~0, scale = 'response') %>%
  mutate(.valueSTD = .value*102.5, assay='Leaf dip') #standardized sp/cm2
post <- bind_rows(post_mb, post_mc)

#get the summarized values
median_hdi_sd(post, .valueSTD) %>%
  mutate_if(is.numeric, signif, 3) %>% 
  write_csv('output/sporangia_both/predictions.csv')

#plot them together
#fix species factor order
#post$species <- factor(post$species, levels = unique(post$species)[c(1:11, 15, 12:14)])
#plot it!
#change control and detached name for figure
post$species2 <- recode_factor(post$species, 
  CONTROL="Inoculum only",
  LIDED="LIDE-D",
  UMCAD="UMCA-D")
# Relevel control to the end
post$species2 <- fct_relevel(post$species2, 'ACMA', 'ARME', 'CEOL', 'HEAR', 'LIDE', 'QUAG', 'QUCH', 'QUPA', 'TODI', 'UMCA', 'LIDE-D', 'PIPO', 'PSME', 'SESE', 'UMCA-D', 'Inoculum only')

#halfeye function. geom_density_ridges plots fake data. use halfeye from tidybayes instead.
post_plot <- post %>% 
  #filter(species!="CONTROL") %>% 
  ggplot(., aes(y = fct_rev(species2), x = .valueSTD)) +
  geom_halfeyeh(.width = .9, size=.1, fatten_point = .1, relative_scale = 2, point_interval = median_hdi) +
  #geom_density_ridges(scale=1.35, lwd=0, panel_scaling = F, rel_min_height = 0.1)+
  #stat_pointintervalh(point_interval = median_hdi, .width =.9, shape=16, size=.1)+
  labs( x=expression(paste("Mean sporangia/", cm^{2})), 
        y='Species') +
  facet_grid(rows = vars(fct_rev(assay)), scales='free_y', space = 'free_y') + 
  scale_x_continuous(limits=c(0, 1250), breaks = seq(0,1250, 250))
post_plot 

ggsave('plots/sporangia_both/predicted_2panels.jpg', post_plot, width = 7, height = 4, units = 'in')

#want latin names for species for poster figure
latin_names <- data.frame(
  species = unique(post$species)[-which(unique(post$species)=="CONTROL")],
  latin = c('Acer macrophyllum', 'Arbutus menziesii', 'Ceanothus oliganthus', 'Heteromeles arbutifolia', 'Notholithocarpus densiflorus', 'Quercus agrifolia', 'Quercus chrysolepis', 'Quercus parvula', 'Toxicodendron diversilobum', 'Umbellularia californica', 'Notholithocarpus densiflorus', 'Pinus ponderosa', 'Pseudotsuga menziesii', 'Sequoia sempervirens', 'Umbellularia california'))
post2 <- left_join(post, latin_names, by='species') %>% 
  filter(species !="CONTROL")
#add in line break to latin name
levels(post2$latin) <- gsub(" ", "\n", levels(post2$latin))
#plot
post_plotlatin <- ggplot(post2, aes(y = fct_rev(latin), x = .valueSTD)) +
  geom_density_ridges(lwd=0, panel_scaling = F)+
  stat_pointintervalh(point_interval = median_hdi, .width = c(.9),shape=16, size=.1)+
  labs( x=expression(paste("Mean sporangia/", cm^{2})), 
        y='Species') +
  facet_grid(rows = vars(fct_rev(assay)), scales='free_y', space = 'free_y') + 
  scale_x_continuous(limits=c(-25, 1250), breaks = seq(0,1250, 250)) +
  theme(axis.text.y = element_text(face = "italic"))
ggsave('plots/sporangia_both/predicted_2panels_latin.jpg', post_plotlatin, width = 7, height = 5.5, dpi = 600, units = 'in')

#calculate contrasts. keep the assays seperate.
#can't figure out how to do it the tidy way, so doing it kinda clunky.

#extract posterior of following coefficients
get_variables(mb)
exb <- mb %>% 
  spread_draws(b_speciesCONTROL, b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA) %>% 
  select(-c(.chain, .iteration, .draw)) %>% as.data.frame()
exc <- mc %>% 
  spread_draws(b_speciesLIDED, b_speciesPIPO, b_speciesPSME, b_speciesSESE, b_speciesUMCAD) %>% 
  select(-c(.chain, .iteration, .draw)) %>% as.data.frame()

#get pairwise comparisons
pairs_b <- combn(1:11,2)
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
         sig=sign(diffsdf_b[,3])==sign(diffsdf_b[,4]),
         assay="broad") %>% 
  mutate_if(is.numeric, signif, 3)
diffsdf_c <- diffsdf_c %>% 
  mutate(sp1=sppc[diffsdf_c[,1]],
         sp2=sppc[diffsdf_c[,2]],
         sig=sign(diffsdf_c[,3])==sign(diffsdf_c[,4]),
         assay="conif")%>% 
  mutate_if(is.numeric, signif, 3)
bind_rows(diffsdf_b, diffsdf_c) %>% 
  write_csv('output/sporangia_both/contrasts.csv') 
