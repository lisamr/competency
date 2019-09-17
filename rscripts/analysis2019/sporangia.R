setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(tidyverse) #do I not need all the other packages below?
library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)
library(scales)
library(ggridges)
library(rethinking)
library(forcats)
library(tidybayes)
library(modelr)

#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
theme_set(theme_bw()) #set ggplot theme

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf assay
df1 <- df %>% filter(leafID<=530) 

###########################################
#DATA VIZ

#plot sporangia counts on a map for select species
df1 %>% filter(spore_assay=="S", trt=="T", species=="TODI") %>%
ggplot(aes(easting, northing)) +
  geom_point(aes(color=count1))+
  scale_color_viridis_c()

#let's try to see how sporangia counts differ across species
#2 choices: average the sporangia counts (gamma distribution?) OR use 3 different counts (mixed model with poisson likelihood)

#start with averages
df1<- df1 %>% rowwise() %>% mutate(countm=mean(c(count1, count2, count3), na.rm = T) ) #rowwise "groups" each row so each mean calucation is unique by row

#plot it. well they aren't all that different!
df1 %>% filter(spore_assay=="S", !countm==0, trt=="T") %>% 
ggplot( aes(species, log(countm)))+
  geom_boxplot()+
  geom_point(alpha=.4)


########################################
#STATISTICAL MODELS?!

#I think it would be most right to treat these counts as COUNTS rather than averaging.
#zero-inflated? poisson mixed model. start with lme4 functions.
df2 <- filter(df1, spore_assay=="S", trt=="T")

#make data tall in regards to counts
dtall <- melt(df2, id.vars = c("species", "leafID", "leafID2"), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count)) %>% 
  arrange(leafID) 
dtall$species <- droplevels(dtall$species)

#I get convergence issues and I forget how to read these model outputs
#m1 <- glmer(count ~ species + (1|leafID), data = dtall, family = "poisson")
#m1.1 <- glmer(count ~ 1 + (1|leafID), data = dtall, family = "poisson")
#summary(m1.1)
#anova(m1, m1.1)
#m2 <- glmer.nb(count ~ species + (1|leafID), data = dtall)
#summary(m2)

########################################
#try with brms first before rethinking
m3 <- brm(
  count ~ -1 + species + (1|leafID),
  data = dtall, family = poisson(),
  chains = 4, cores = 4,
  control = list(adapt_delta = .99, max_treedepth=15))
#no warnings :)

#check out coefs
summary(m3)
fixef(m3)
launch_shinystan(m3)
#population level effects only. traceplots look ok
plot(m3, pars = "^b_")
postm3 <- posterior_samples(m3, pars = c("^b_", "sd"))
head(postm3)
pp_check(m3)

#check out posterior with tidybayes
get_variables(m3) #b_ are fixed, r_ are rand
get_variables(m3)[1:11]
#would be nice to know how to include a vector of params, but dont know how atm.

#see model fit
#simulate prediction data using same dataset to see fit
m3sim <- add_predicted_draws(dtall, m3) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(30) 
#plot against your own data
#pdf('plots/sporangia/modelfit_brms.pdf')
m3sim %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.5, color='grey')+
  stat_density(data=dtall, aes(x=count), geom="line", color='slateblue')+
  facet_wrap(~species, scales = 'free')
#dev.off()  

#coef plot
#pdf('plots/sporangia/coefplot_brms.pdf')
m3 %>% 
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA, sd_leafID__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')
#dev.off()

#summarize coefs
m3 %>% 
  gather_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA, sd_leafID__Intercept) %>%
  median_qi(.width=.9)

#contrasts between species
#extract posterior of following coefficients
ex <- m3 %>% 
  spread_draws(b_speciesACMA, b_speciesARME, b_speciesCEOL, b_speciesHEAR, b_speciesLIDE, b_speciesQUAG, b_speciesQUCH, b_speciesQUPA, b_speciesTODI, b_speciesUMCA) %>% 
  select(-c(.chain, .iteration, .draw)) %>% as.data.frame()
head(ex)

#can't figure out how to do it the tidy way, so doing it kinda clunky

pairs <- combn(1:10,2)
f <- function(x){
  spdiff <-  (ex[,pairs[1,x]]-ex[,pairs[2,x]])
  PI(spdiff, .9)
}
diffs <- sapply(1:ncol(pairs), f)
#put contrasts into a readable dataframe
spp <- colnames(ex) %>% str_replace('b_species', '')
diffsdf <- rbind(pairs, diffs) %>% t %>% as.data.frame
diffsdf <- diffsdf %>% 
  mutate(sp1=spp[diffsdf[,1]],
         sp2=spp[diffsdf[,2]],
         sig=sign(diffsdf[,3])==sign(diffsdf[,4]))
diffsdf

#get predicted fit based on just the species coef.
postp <- dtall %>%
  modelr::data_grid(species) %>%
  add_fitted_draws(m3, re_formula=~0, scale='response') %>%
  #standardized sp/cm2
  mutate(.valueSTD = .value*102.5/1.227) 

#pdf('plots/sporangia/predictions_brms.pdf')
ggplot(postp, aes(y = fct_rev(species), x = .valueSTD)) +
  geom_density_ridges(lwd=.1)+
  stat_pointintervalh(point_interval = median_hdcih, .width = c(.9), size=.2)+
  labs(x='# sporangia', y='Species')
#dev.off()
#get the values
postp %>% 
  median_qi(x=.valueSTD, .width = c(.9, .95))

write.csv(results, 'output/sporangia/backtransformed_results_brms.csv', row.names = F)
########################################
#use rethinking package?
#goal: create a model with a poisson likelihood and leafID as random intercept
library(rethinking)

#analyze data with stan
#make sure index variables are consecutive integers
dat2 <- list(
  species=as.integer(factor(dtall$species)),
  leafID=as.integer(factor(dtall$leafID)),
  count=dtall$count
)

#figure out prior for the species intercept
curve( dlnorm( x, 3, 1.5) , from=0 , to=100 , n=200, xlab = "mean # spores (lambda)")
exp(rnorm(1000, 3, 1.5)) %>% PI(.5)
curve( dlnorm( x, 0, 5) , from=0 , to=50 , n=200, xlab = "deviation around the mean")
#b <- exp(rnorm(5000, 0, 1.5))
#mean(b)
#dens(b, xlim=c(0,100))

#model it
m4 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- a[species]+ b[leafID],
    #priors
    a[species] ~ dnorm(2, 1.5), 
    #adaptive priors
    b[leafID] ~ dnorm(0, sigma_b),
    sigma_b ~ dexp(1)
  ), data=dat2, chains=4, log_lik = T
)
traceplot(m4) #the alphas dont look amazing
par(mfrow=c(1,1))
stancode(m4)#check out stan code
ptablem4 <-  precis(m4, depth = 2, pars=c('a', 'sigma_b'), prob = .9) #%>% plot
parnames <- c('ACMA', 'ARME', 'CEOL', 'HEAR', 'LIDE', 'QUAG', 'QUCH', 'QUPA', 'TODI', 'UMCA', 'sigma_leafID')
plot(ptablem4, labels=parnames, pch=16)
#write.csv(round(ptablem4, 2), 'output/sporangia/precis_m4.csv')

#contrast it with intercept only model. isn't that different actually.
dat2.1 <- dat2
dat2.1$species <- NULL
m4.1 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- a + b[leafID],
    #priors
    b[leafID] ~ dnorm(0, sigma_b),
    a ~ dnorm(2, 1.5),
    sigma_b ~ dexp(1)
  ), data=dat2, chains=4, log_lik = T
)
precis(m4.1)
traceplot_ulam(m4.1, pars = 'a')#looks fine
dev.off()
#compare with WAIC
waic_comp <- compare(m4.1, m4)
waic_comp
dwaic <- rnorm(10000, 12.3, 9.43)
dens(dwaic, show.zero = T,show.HPDI = .95)
#compare with LOO
loo_comp <- compare(m4.1, m4, func = LOO)
loo_comp
dloo <- rnorm(10000, 21.8, 10.56)
dens(dloo, show.zero = T,show.HPDI = .95)

#try to reparameterize it so it's non-centered? hopefully it'll sample better.
#model it
m5 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- a[species] + z[leafID]*sigma_b, 
    #priors
    a[species] ~ dnorm(2, 1.5), #species intercepts
    z[leafID] ~ dnorm(0, .5), #variation in subsamples
    sigma_b ~ dexp(1)
  ), data=dat2, chains=4, iter = 2000, cores = 4
)
spp <- as.vector(unique(dtall$species))
precis(m4, depth = 2, pars=c('a', 'sigma_b'))
precis(m5, depth=2, pars=c('a', 'sigma_b', 'z')) #sampling definitely better, but still low for some

#compare both models. pretty comparable and m4 is more interpretable.
CT <- coeftab(m4, m5) 
CTc <- CT@coefs
CTse <- CT@se
colnames(CTse) <- colnames(CTc)
CTd <-
  left_join(
    as.data.frame(CTc) %>%
      mutate(param = row.names(CTc)) %>%
      tidyr::gather(model, estimate, 1:ncol(CTc)),
    as.data.frame(CTse) %>%
      mutate(param = row.names(CTse)) %>%
      tidyr::gather(model, se, 1:ncol(CTse))
  )
CTd %>% filter(grepl('a', param)) %>% 
  ggplot(., aes(param, estimate, group=model, color=model)) +
  #geom_col(position = "dodge") +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), position="dodge")

#both predict really well
postcheck(m5)
postcheck(m4)

#go with m4.
precis(m4, depth = 2, pars=c('a', 'sigma_b')) %>% plot(labels=c(spp, "sigma_ind"))

#check out how well the model fits the data. predict new data points from model and contrast with real data.
pred <- sim(m4) #951 samples (cols) simulated 1500 times (rows)
#all of the samples
dens(pred[1,], col=alpha('black', .1))
lapply(1:500, function(x)dens(pred[x,], add = T, col=alpha('black', .05)))
dens(dtall$count, add = T, col='blue')#observed data

#split among species
colnames(pred) <- dtall$species
fdens <- function(species, xmax=200, ymax=.2, nsim=100){
  use <- which(colnames(pred)==species)
  #plot predicted lines
  plot(NULL, xlim=c(0,xmax), ylim=c(0,ymax), xlab="#spores", ylab="density", main=species)
  lapply(1:nsim, function(x) dens(pred[x, use], add=T, col=alpha('black', .1)))
  #plot observed line
  dens(dtall$count[dtall$species==species], add = T, col='blue')
}

#viz model fits for each species
#wow, those are some amazing model fits!
#pdf('plots/sporangia/post_pred.pdf', width = 10, height = 5)
par(mfrow=c(2,5))
fdens("ACMA", xmax=250, ymax=.5)
fdens("LIDE", ymax=.15, xmax=50)
fdens("UMCA", ymax=.1, xmax=100)
fdens("TODI", ymax=1, xmax=20)
fdens("QUCH", ymax=1, xmax=20)
fdens("QUAG", ymax=1.5, xmax=20)
fdens("QUPA", ymax=1, xmax=20)
fdens("HEAR", ymax=2, xmax=10)
fdens("ARME", ymax=3, xmax=10)
fdens("CEOL", ymax=3, xmax=10)
#dev.off()
reset()

#understand what the model is saying...
#pdf('plots/sporangia/alpha_coef.pdf', width = 6, height = 4)
precis(m4, depth = 2, pars = c("a")) %>% plot(labels=c(spp), xlab="species coef (log-sporangia)")
#dev.off()
####################
#contrasts between species intercept coefficients
#ex$a %>% head #posterior of species effects
pairs <- combn(1:10,2)
ex <- extract.samples(m4)
f <- function(x){
  spdiff <-  ex$a[,pairs[1,x]]-ex$a[,pairs[2,x]]
  PI(spdiff, .95)
}
diffs <- sapply(1:ncol(pairs), f)
#put contrasts into a readable dataframe
diffsdf <- rbind(pairs, diffs) %>% t %>% as.data.frame
diffsdf <- diffsdf %>% 
  mutate(sp1=spp[diffsdf[,1]],
         sp2=spp[diffsdf[,2]],
         sig=sign(diffsdf[,3])==sign(diffsdf[,4]))
diffsdf
#write.csv(diffsdf, row.names = F, "~/Box/Competency project/competency.git/output/sporangiacontrasts.csv")
##############################

#compare to data. not sure why alpha is so different from observed data. maybe there's just a lot of intraindividual variation? 
head(df2)
summ <- df2 %>% group_by(species) %>% summarise(obsmean=mean(countm, na.rm=T))
am <- sapply(1:10, function(x) mean(exp(ex$a[,x]))) 
api <- sapply(1:10, function(x) PI(exp(ex$a[,x]))) %>% t
cbind(summ, am, api)
####################
#get back transformed mean and sd values of counts
#simulate data and get the means and sd???
#log(lambda) <- a[species] + b[leafID]
loglam <- with(ex, (a + rnorm(nrow(a), 0, sigma_b)))
praw <- exp(loglam)
head(praw)
plam <- apply(praw, 2, mean)
pPI <- apply(praw, 2, PI, .9) %>% t %>% as.data.frame
names(pPI) <- c("lo5", 'hi95')
psd <- apply(praw, 2, sd)
estcounts <- data.frame(summ, plam, pPI, psd)
estcounts <- data.frame(species=spp, apply(estcounts[,-1], 2, round, 2))
#write.csv(estcounts, 'output/sporangia_postmeans.csv', row.names = F)

#plot 1
#pdf('plots/sporangia/lam_ptrange.pdf', width=6, height=4)
ggplot(estcounts, aes(species, plam)) +
  geom_pointrange(aes(ymin=lo5, ymax=hi95), shape=22, fill='white') +
  ylim(c(0,75)) +
  labs(y='90th PI mean counts (lambda)')
#dev.off()

#plot 2
par(mfrow=c(2,5))
lapply(1:10, function(x) dens(praw[,x], main=spp[x], xlab="lambda"))
reset()

#plot 3
praw2 <- as.data.frame(praw)
colnames(praw2) <- spp
prawtall <- melt(praw2, variable.name = "species", value.name = "lambda")
prawtall$speciesrev <- factor(prawtall$species, levels = rev(levels(prawtall$species)))
#pdf('plots/sporangia/lam_ridges.pdf', width = 10, height = 5)
ggplot(prawtall, aes(lambda, speciesrev)) +
  geom_density_ridges(aes(fill=species), alpha=.9) +
  xlim(c(0,30))+
  scale_fill_viridis_d()+
  labs(y='species')
#dev.off()

##########################################
#RELATIONSHIP WITH LESION SIZE?

#visualize first
#make data tall in regards to counts
dtall2 <- melt(df2, id.vars = c("species", "leafID", "perc_lesion"), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count)) %>% 
  arrange(leafID) 
head(dtall2)

#all counts
ggplot(dtall2, aes(count, perc_lesion))+
  geom_point()+
  facet_wrap(~species, scales = "free_x")
#do it with the mean counts
df1 %>% filter(trt=="T") %>% 
ggplot(aes(countm, perc_lesion))+
  geom_point()+
  facet_wrap(~species, scales = "free_x") +
  labs(x="mean # sporangia")

##########################################
#RELATIONSHIP WITH chlamydospores? theres a paucity of species (6) with chlamydo counts.
df1 %>% filter(spore_assay=="C", !is.na(count1)) %>% count(species, trt)

#do it with the mean counts. need to make dataframe wide
df3 <- df1 %>% rowwise() %>% 
  mutate(count_est=case_when(spore_assay=="C"~count1/prop,
                   spore_assay=="S"~countm)) %>% 
  select(species, ind, trt, spore_assay, count_est) %>% 
  dcast(species+ind+trt~spore_assay, value.var = "count_est")
head(df3)

df3 %>% filter(trt=="T", !is.na(C)) %>% 
ggplot(aes(S, C))+
  geom_point()+
  facet_wrap(~species, scales = "free")

