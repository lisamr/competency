#is there a relationship between ubiquity in the plot network and number of sporangia?
setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(lme4)
library(merTools)
library(reshape2)
library(scales)
library(ggridges)
library(rethinking)
library(forcats)
library(cowplot)
library(tidyr)
#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
theme_set(theme_classic())#set ggplot theme

#table to describe number of plots each species is present
ubiq <- data.frame(
  species=c("UMCA", "QUAG", "LIDE", 'ARME', 'QUPA', 'QUCH', 'TODI', 'HEAR', 'CEOL', 'PIPO', 'SESE', 'ACMA', 'PSME'),
  mixed = c(116, 101, 83, 73, 56, 52, 35, 29, 25, 16, 5, 11, 6),
  red = c(68, 12, 86, 18, 26, 2, 8, 4, 10, 0, 111, 14, 6))

#viz plot data
p1=ggplot(ubiq, aes(fct_reorder(species, -mixed), mixed)) +
  geom_col(fill='darkolivegreen')+
  labs(title='Mixed Evergreen Forests', x='species', y='# plots present')
p2=ggplot(ubiq, aes(fct_reorder(species, -red), red)) +
  geom_col(fill='firebrick4')+
  labs(title='Redwood Forests', x='species', y='# plots present')
plot_grid(p1, p2, nrow = 2)

#merge that data into the masterfile
#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf assay
df1 <- df %>% 
  filter(spore_assay=="S", trt=="T") %>% 
  left_join(ubiq, by='species') %>% 
  melt(id.vars = c("species", "leafID", "leaf_area_cm2", "mixed", 'red'),
       measure.vars = c("count1", "count2", "count3"),
       variable.name = "sample", value.name = "count") %>% 
  rename(area=leaf_area_cm2) %>% 
  drop_na() %>% #filters out rows with na
  arrange(leafID) 

#crude summary
df1 %>% group_by(species) %>% summarise(ubiq=mean(mixed), count=mean(count)) %>% 
  ggplot(., aes(ubiq, count)) +
  geom_point() +
  geom_smooth(method='lm')

#quick plot
ggplot(df1, aes(mixed, (count))) +
  geom_point()+
  geom_smooth(formula = y~x, method = 'glm')
df1 %>% group_by(species, leafID, mixed) %>% summarise(count=mean(count)) %>% 
  ggplot(., aes(mixed, count))+
  geom_point()+
  geom_smooth(formula = y~x, method = 'glm')

######################
#filter out na's
#model the sporangia against number of plots 
f1 <- bf(count ~ mixed + offset(log(area)) + (1|leafID) + (1|species), family = poisson())
f2 <- bf(count ~ red + offset(log(area)) + (1|leafID) + (1|species), family = poisson())
f3 <- bf(count ~ 1 + offset(log(area)) + (1|leafID) + (1|species), family = poisson())

get_prior(f1, df1)
prior1 <- c(
  set_prior("normal(0, 2)", class="Intercept"),
  set_prior("normal(0, .5)", class="b"),
  set_prior("exponential(1)", class = "sd") )
prior2 <- c(
  set_prior("normal(0, 2)", class="Intercept"),
  set_prior("exponential(1)", class = "sd") )

#mixed
m1m <- brm(formula = f1, prior = prior1,
  data = df1, family = poisson(), iter = 2000,
  chains = 4, cores = 4,
  control = list(adapt_delta = .95, max_treedepth=15))
#red
m1r <- brm(
  formula = f2, prior = prior1,
  data = df1, family = poisson(),
  chains = 4, cores = 4,
  control = list(adapt_delta = .95, max_treedepth=15))
#intercept only
m0 <- brm(
  formula = f3, prior = prior2,
  data = df1, family = poisson(),
  chains = 4, cores = 4,
  control = list(adapt_delta = .95, max_treedepth=15))

#compare models
m1m <- add_criterion(m1m, c('loo', 'waic'))
m1r <- add_criterion(m1r, c('loo', 'waic'))
m0 <- add_criterion(m0, c('loo', 'waic'))
loo_compare(m1m, m1r, m0)
LOO(m1m, m1r, m0)


#check out models
plot(m1m, pars = "^b_")
summary(m1m)
plot(m1r, pars = "^b_")
summary(m1r)
pp_check(m1m)
pp_check(m1r)

#coef values
#mixed
m1m %>% 
  gather_draws(b_Intercept, b_mixed, sd_leafID__Intercept, sd_species__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')
m1m %>% 
  gather_draws(b_Intercept, b_mixed, sd_leafID__Intercept, sd_species__Intercept) %>%
  median_hdci(.width=c(.9, .95))
#redwood
m1r %>% 
  gather_draws(b_Intercept, b_red, sd_leafID__Intercept, sd_species__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')
m1r %>% 
  gather_draws(b_Intercept, b_red, sd_leafID__Intercept, sd_species__Intercept) %>%
  median_hdci(.width=c(.9, .95))

#see model prediction
#get posterior predictions
xx <- seq(0,120,length.out = 121)
#mixed
newm <- data.frame(mixed=xx, area=1, leafID=1, species="UMCA")
pp1 <- fitted(m1m, newdata = newm, scale = 'response', re_formula = ~0, summary = F)
pp1m <- apply(pp1, 2, median)
pp1pi <- apply(pp1, 2, PI, .9)
#redwood
newr <- data.frame(red=xx, area=1, leafID=1, species="UMCA")
pp2 <- fitted(m1r, newdata = newr, scale = 'response', re_formula = ~0, summary = F)
pp2m <- apply(pp2, 2, median)
pp2pi <- apply(pp2, 2, PI, .9)

#plot it!
#mixed
pdf("plots/ubiquity/ubiquity_sporangia_mixed.pdf", 8, 6)
plot(df1$mixed, df1$count/df1$area, col=alpha('slateblue', .2), pch=16,ylim=c(-5,85), xlab='# plots (mixed)', ylab='#sporangia')
abline(a=0, b=0, lty=2)
lines(xx, pp1m, col=alpha('slateblue4', .9))
shade(pp1pi, xx, col=alpha('slateblue4', .3))
dev.off()
#redwood
pdf("plots/ubiquity/ubiquity_sporangia_red.pdf", 8, 6)
plot(df1$red, df1$count/df1$area, col=alpha('slateblue', .2), pch=16,ylim=c(-5,85), xlab='# plots (mixed)', ylab='#sporangia')
abline(a=0, b=0, lty=2)
lines(xx, pp2m, col=alpha('slateblue4', .9))
shade(pp2pi, xx, col=alpha('slateblue4', .3))
dev.off()
