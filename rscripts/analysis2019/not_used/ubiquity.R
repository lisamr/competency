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
#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))
theme_set(theme_bw())#set ggplot theme

#table to describe number of plots each species is present
ubiq <- data.frame(
  species=c("UMCA", "QUAG", "LIDE", 'ARME', 'QUPA', 'QUCH', 'TODI', 'HEAR', 'CEOL', 'PIPO', 'SESE', 'ACMA'),
  mixed = c(116, 101, 83, 73, 56, 52, 35, 29, 25, 16, 5, 11),
  red = c(68, 12, 86, 18, 26, 2, 8, 4, 10, 0, 111, 14))

#merge that data into the masterfile
#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf assay
df1 <- df %>% 
  filter(leafID<=530, spore_assay=="S", trt=="T") %>% 
  left_join(ubiq, by='species') %>% 
  melt(id.vars = c("species", "leafID", "mixed", 'red'),
       measure.vars = c("count1", "count2", "count3"),
       variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count)) %>% 
  arrange(leafID) 
#add plot rank too
sprank <- rank(unique(df1$mixed))
df1$mixedr <- sprank[as.integer(factor(df1$species))]

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

################################################
#count ~ ubiq + (1|sp) + (1|ID)
m0 <- glmer(count ~ mixed + (1|species) + (1|leafID), data=df1, family='poisson')
summary(m0)
#how do you see the distribution of all these params?
ranef(m0)$species 
ranef(m0)$leafID %>% dens
fixef(m0)
str(m0)
exp(fixef(m0))#slope on response scale

#check out fit. it's good
plot(predict(m0, type='response'), df1$count)+abline(0,1)

#plot model fit against data
xseq <- seq(0,120, length.out = 1000)
pr <- predictInterval(m0, data.frame(mixed=xseq, species="UMCA", leafID=1), which='fixed', type = 'linear.prediction')
pr$mixed <- xseq 
pr2 <- predictInterval(m0, data.frame(mixed=xseq, species="UMCA", leafID=1), which='fixed', type = 'probability')
pr2$mixed <- xseq 

#plot1
ggplot(df1, aes(mixed, log(count))) +
  geom_point() +
  geom_line(data=pr, aes(mixed, fit), col='grey') +
  geom_line(data=pr, aes(mixed, upr), lty=2, col='grey') +
  geom_line(data=pr, aes(mixed, lwr), lty=2, col='grey')
#plot2
ggplot(df1, aes(mixed, count)) +
  geom_point(alpha=.2) +
  geom_line(data=pr2, aes(mixed, fit), col='grey') +
  geom_line(data=pr2, aes(mixed, upr), lty=2, col='grey') +
  geom_line(data=pr2, aes(mixed, lwr), lty=2, col='grey')+
  ylim(c(0,90)) #omit acma outliers for plot
  #scale_y_continuous(trans = 'log2')

################################################
#model it in stan with brms
#intercept only
#mbrm0 <- brm(count ~ 1 + (1|species) + (1|leafID),data = df1, family = poisson(),chains = 3, cores = 4, control = list(adapt_delta = .95, max_treedepth=15))
#n.plots (mixed) included in model
mbrm <- brm(
  count ~ mixed + (1|species) + (1|leafID),
  data = df1, family = poisson(),
  chains = 3, cores = 4,
  control = list(adapt_delta = .95, max_treedepth=15))

#Warning message:
#  Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable. Running the chains for more iterations may help. See http://mc-stan.org/misc/warnings.html#bulk-ess 
#see model summary and coefs. kinda just care about 'mixed'
summary(mbrm)
launch_shinystan(mbrm)

#coefficient plot
parnames(mbrm)
mbrm %>% 
  gather_draws(b_Intercept, b_mixed, sd_leafID__Intercept, sd_species__Intercept) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')
mbrm %>% 
  gather_draws( b_mixed) %>%
  rename(par=.variable, value=.value) %>% 
  ggplot(aes(y = fct_rev(par), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient')
mbrm %>% 
  gather_draws(b_Intercept, b_mixed, sd_leafID__Intercept, sd_species__Intercept) %>%
  median_qi(.width=c(.9, .95))

#compare models
#mbrm0 <- add_criterion(mbrm0, c("loo", "waic", 'kfold'))
#mbrm <- add_criterion(mbrm, c("loo", "waic", 'kfold'))
#loo_compare(mbrm, mbrm0, criterion = 'kfold')

#check out model fit
brmsim <- add_predicted_draws(df1, mbrm) %>% 
  select(-.chain, -.iteration) %>% 
  group_by(species, .draw) %>% 
  sample_draws(30) 
brmsim %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.5, color='grey')+
  stat_density(data=df1, aes(x=count), geom="line", color='slateblue')+
  facet_wrap(~species, scales = 'free')

#see model prediction
#plot spores in relation to # plots, not? accounting for random variation
#get posterior predictions
xx <- seq(0,120,length.out = 121)
newd <- data.frame(mixed=xx, leafID=1, species="UMCA")
pp1 <- fitted(mbrm, newdata = newd, scale = 'response', re_formula = ~0, summary = F)#the fitted line, not predicted point
pp1m <- apply(pp1, 2, mean)
pp1pi <- apply(pp1, 2, PI, .9)

#plot it!
#plot points averaged within ind
#d <- df1 %>% group_by(species, leafID, mixed) %>% summarise(count=mean(count)) 
#plot(d$mixed, d$count, col=alpha('slateblue', .2), pch=16, ylim=c(0,85))
#or plot raw data
pdf('plots/ubiquity_sporangia_brms.pdf')
plot(df1$mixed, df1$count, col=alpha('slateblue', .2), pch=16, 
     ylim=c(0,85), xlab='# plots (mixed)', ylab='#sporangia')
lines(xx, pp1m, col=alpha('slateblue4', .9))
shade(pp1pi, xx, col=alpha('slateblue4', .1))
dev.off()  

################################################
#model it in stan with rethinking.
head(df1)
dat <- list(spID = as.integer(factor(df1$species)),
            leafID = as.integer(factor(df1$leafID)),
            ubiq = as.integer(df1$mixed),
            count = df1$count)
m1 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- as[spID] + al[leafID] + b*ubiq,
    #adaptive priors
    as[spID] ~ dnorm(abar, sigmaa), 
    al[leafID] ~ dnorm(0, sigmal),
    #fixed priors
    abar ~ dnorm(2, 1.5),
    b ~ dnorm(0, .5),
    sigmaa ~ dexp(1),
    sigmal ~ dexp(1)
  ), data=dat, chains=3, cores=4
)
precis(m1) #super bad sampling

####################################
#try with rstanarm
library(rstanarm)
m1.1 <- stan_glmer(count ~ mixed + (1|species) + (1|leafID),
         data = df1,
         family = poisson(), 
         chains = 3, cores = 4, seed = 1)

parslist <- c('(Intercept)', 'mixed', "Sigma[leafID:(Intercept),(Intercept)]", "Sigma[species:(Intercept),(Intercept)]")
plot(m1.1, pars=parslist)
ci95 <- posterior_interval(m1.1, prob = 0.9, pars=parslist)
round(ci95, 2)
launch_shinystan(m1.1)

#posterior with same data. its good.
p1 <- posterior_predict(m1.1)
plot(dat$count, p1[1,])
plot(NULL, xlim=c(0,100), ylim=c(0, .5))
for(i in 1:200){
  dens(p1[i,], add=T, col=alpha("grey", .1))
}
dens(dat$count, add=T, col="blue")

#posterior with new data
dnew <- data.frame(mixed=xseq, species="UMCA", leafID=1)
p2 <- posterior_predict(m1.1, dnew, re.form=NA)
pm <- apply(p2, 2, mean)
pPI <- apply(p2, 2, PI)
plot(df1$mixed, df1$count, ylim=c(0,85))
lines(xseq, pm)
lines(xseq, pPI[1,], lty=3)
lines(xseq, pPI[2,], lty=3)

#another way to predict over a range with extracted posterior samples.
ex2 <- as.matrix(m1.1)
str(ex2)
dimnames(ex2)
int <- ex2[,1]
mi <- ex2[,2]
sID <- ex2[,330]
ss <- ex2[,331]
plot(NULL, ylim=c(0,85), xlim=c(0,120), xlab='# plots present', ylab='# sporangia')
#plot ignoring random effects
for(i in 1:300){
  curve(exp(int[i]+mi[i]*x), 0, 120, 100, add = T, col=alpha('grey20', .1))
}
points(df1$mixed, df1$count, ylim=c(0,85), pch=16, col=alpha('blue', .2))

#another way to plot using shade!
xx <- seq(0,120,length.out = 1000)
mum <- sapply( xx , function(x) mean(exp(int + mi*x)))
muPI <- sapply( xx, function(x) PI(exp(int + mi*x), .95))
pPI <- sapply( xx, function(x) PI(exp(int + mi*x + rnorm(length(int), 0, ss)), .95))#predictive with random effex
#pdf('plots/ubiquity_sporangia.pdf', 6, 4)
plot(df1$mixed, df1$count, ylim=c(0,85), xlim=c(0,120),  pch=16, col=alpha('slateblue', .2), xlab='# plots present', ylab='# sporangia', bty='n')
shade(muPI, xx)
#shade(pPI, xx)
lines(xx, mum, col='slateblue4', lty=1)
#dev.off()

#look at distribution of the slope parameter "mixed"
parm <- exp(mi) #backtransform?
dens(parm, show.HPDI = .9, show.zero = T)
dens(parm, show.HPDI = .95, add=T)
