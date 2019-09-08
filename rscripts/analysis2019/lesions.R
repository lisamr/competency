##########################################
# RELATIONSHIP BETWEEN SPORANGIA AND LESION SIZE?

setwd("~/Box/Competency project/competency.git")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)
library(scales)
library(ggridges)
library(rethinking)
library(brms)
library(brmstools)

#for turning pars back 
reset <- function(x) par(mfrow=c(1,1))

#read in master file
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just broad leaf treatment assay
df1 <- df %>% filter(leafID<=530, trt=="T", spore_assay=="S") 

#visualize first
#make data tall in regards to counts
dtall2 <- melt(df1, id.vars = c("species", "leafID", "perc_lesion"), measure.vars = c("count1", "count2", "count3"), variable.name = "sample", value.name = "count") %>% 
  filter(!is.na(count)) %>% 
  mutate(perc_lesion=perc_lesion/100) %>% 
  arrange(leafID) 
head(dtall2)

#all counts vs lesions
ggplot(dtall2, aes(perc_lesion, count))+
  geom_point()+
  facet_wrap(~species, scales = "free_y")

#just mean counts vs lesions
df2 <- dtall2 %>% group_by(species, leafID) %>% summarise(lesion=mean(perc_lesion), countm=mean((count))) #log counts or not?
ggplot(df2, aes(lesion, countm))+
  geom_point()+
  facet_wrap(~species, scales = "free_y")
  #geom_smooth(method = 'lm') #probably need to remove ACMA's outlier because it has too much leverage


########################################
#try it with glmer first to get an idea of what you want
m0 <- glmer(count ~ perc_lesion + (1|leafID) + (perc_lesion|species), data = dtall2, family = 'poisson')
summary(m0)

#thats cool it worked. Not sure how you check out model performance.
plot(m0)
plot(dtall2$count, predict(m0, type='response'))+abline(0,1)

#random slopes for each species
ranef(m0)$species #how do you get the variation across those values
plot(ranef(m0)$species[,2], 1:10)
#library(lattice)
dotplot(ranef(m0, postVar = TRUE))#plots prediction intervals of the random slopes and intercepts

####################################
#try to model it with brms

m4 <- brm(
  count ~ perc_lesion + (1|leafID) + (perc_lesion|species),
  data = dtall2, family = poisson(),
  chains = 3, cores = 4, iter=3000,
  control = list(adapt_delta = .99, max_treedepth=15))
#saveRDS(m4, 'output/lesion/brmsmodel_m4.rds')#save model for easy loading in future
#stancode(m4) #see stancode!
#Warning messages:
#1: There were 3 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See
#http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
#3: Examine the pairs() plot to diagnose sampling problems
#plots pop-level effects only
plot(m4, pars = "^b_")

#see model summary and output
sum4 <- summary(m4);sum4
#check out random slopes. values are the entire coef, which includes fixed + rand effect.
#warning 1: 'forest' is deprecated.Use 'tidybayes' instead.
pdf('plots/lesions/forest.pdf',width = 6,4)
forest(m4, pars='perc_lesion', grouping='species', level=.9, sort = F, digits = 1, av_name = 'mean slope')+
  geom_vline(xintercept = 0, lty=2) +
  geom_point(color='slateblue') 
dev.off()
#hmm...don't match up with reported ranef from brms. that's because ranef is the difference from the population level effect. they do match up if you do:
#coef(m4)
REm4 <- ranef(m4, probs = c(.05, .95))
REm4$species[,,2]

#check out model predictions
#model fit
pp1 <- predict(m4, summary = F)
dens(dtall2$count, col='white')
sapply(1:200, function(x) dens(pp1[x,], add=T, col=alpha('grey', .1)))
dens(dtall2$count, col='blue', add=T)
pp_check(m4) #easy function to do the same

#predict model fit for each species
xx <- seq(0,1,length.out = 100)
newd <- data.frame(perc_lesion=xx, leafID=1, species=rep(unique(dtall2$species), each=length(xx)))
pp2 <- predict(m4, newdata = newd, re_formula = ~(1|species), summary = F)
pp3 <- fitted(m4, newdata = newd, scale = 'response', re_formula = ~(1|species), summary = F)#the fitted line, not predicted points
dim(pp3)
head(pp3)

plotmodfit <- function(i, pp, ...){
  sp <- as.vector(unique(newd$species))
  use <- as.integer(factor(newd$species))==i
  plot(dtall2$perc_lesion[dtall2$species==sp[i]], dtall2$count[dtall2$species==sp[i]], ...,
       xlab='lesion', ylab='# sporangia', main=sp[i], 
       pch=16, col=alpha('slateblue', .5), xlim=c(0,1))
  apply(pp[,use], 2, median) %>% lines(xx, .)
  apply(pp[,use], 2, PI, .9) %>% shade(., xx)
}
#wierd. looks like the same shitty predictions from rethinking package.
pdf('plots/lesions/predictions_brms.pdf',8, 6)
par(mfrow=c(2,5))
sapply(1:10, function(x) plotmodfit(x, pp3))
reset()
dev.off()

###

























########################
#statistical models!
#model with rethinking isn't producing what I think it should
#figure out priors. modeling counts ~ species + lesion 

#prior for intercept. mass should be under 20 or so, which is the maxish spores produced regardless of lesion area.
a <- rnorm(5000, 3, 1.5)
a2 <- exp(a)
dens(a2, xlim=c(0,100), xlab="sp intercept")
mean(a2)

#prior for slope. probably shouldnt go down (or much at least)
b <- rnorm(5000, 1, .5)
sp <- seq(0,1, length.out = 100)
lam <- sapply(1:length(sp), function(x) exp(a + b*sp[x]))
plot(NULL, xlim=c(0,1), ylim=c(0,500), xlab="lesion", ylab="#spores")
lapply(1:200, function(x) lines(sp, lam[x,], col=alpha(1, .3)))

#prior for rho, the correlation between int and slopes. should be centered around zero.
R <- rlkjcorr(5000, 2, 2) #array
dens(R[,1,2], xlab='correlation')
###################################
###################################
dat2 <- list(
  count=dtall2$count,
  species=as.integer(factor(dtall2$species)),
  ID=as.integer(factor(dtall2$leafID)),
  lesion=dtall2$perc_lesion)
#no leafID
m2 <- ulam(
  alist(
    #likelihood and linear model
    count ~ dpois(lambda),
    log(lambda) <- a[species] + g[ID] + b[species]*lesion,
    #varying effects
    g[ID] ~ dnorm(0, sigmag),
    c(a, b)[species] ~ multi_normal(c(abar, bbar), rho, sigma),
    #fixed priors
    sigmag ~ dexp(1),
    abar ~ dnorm(3, 1.5),
    bbar ~ dnorm(1, .5),
    rho ~ dlkjcorr(2),
    sigma ~ dexp(1)
  ), 
  control=list(adapt_delta=0.99),
  data=dat2 , chains=4 , cores=4 )

#HUGE NOTE: MODEL WAS SAMPLING VERRRRY POORLY WHEN MU FOR G[ID] HAD A VALUE OTHER THAN ZERO. THIS IS BECAUSE THE MODEL HAD IDENTIFIABILITY ISSUES. SET IT TO ZERO!
#not sure why it's coming up with different estimates on the random slopes than the glmer model.
traceplot(m2) 

postcheck(m2)
reset()
parsm2 <- m2@pars
sp <- c('ACMA', 'ARME', 'CEOL', 'HEAR', 'LIDE', 'QUAG', 'QUCH', 'QUPA', 'TODI', 'UMCA')
parnames <- c(paste0("b_", sp), paste0("a_", sp), "sigmag", "abar",   "bbar",   rep("rho", 4), rep("sigma", 2))
precis(m2, pars = parsm2[-1], depth = 3) %>% plot(labels=parnames)
####################################
#check out the rethinking model!
#model predictions with observed data overlayed
pred <- sim(m2)
spnames <- unique(df2$species) %>% as.vector()
dens(pred[1,], col=alpha('black', .1))
lapply(1:500, function(x)dens(pred[x,], add = T, col=alpha('black', .05)))
dens(dat2$count, add = T, col='blue')#observed data

fplot1 <- function(i){
  sp <- dat2$species==i
  plot(dat2$lesion[sp], dat2$count[sp], main=spnames[i])
}
par(mfrow=c(2,5))
lapply(1:10, fplot1)
reset()

#get posterior predictions
ex <- extract.samples(m2)
str(ex)
#need to predict values with my own function using the structure of the model.
# log(lambda) <- a[species] + g[ID] + b[species]*lesion,
#g[ID] ~ dnorm(0, sigmag)
#c(a, b)[species] ~ multi_normal(c(abar, bbar), rho, sigma),
ex$a[,1] 
str(ex)
lesionseq <- seq(0,1,length.out = 1000)

#function for predicting posterior
fpost <- function(sp, link=T){
  postraw <- matrix(NA, ncol=length(lesionseq), nrow = nrow(ex$a))
  gbar <- 0
  #gbar <- mean(rnorm(nrow(ex$g), 0, ex$sigmag))
  for(i in 1:length(lesionseq)){
    lesion <- lesionseq[i]
    mean(ex$g)
    #"link" function. will not include intra-ind. variation
    if(link==T){
      postraw[,i] <- with(ex, (a[,sp] + gbar + b[,sp]*lesion)) %>% exp 
      #"sim" function. includes intra-ind variation
    }else{
      postraw[,i] <- with(ex, (a[,sp] + rnorm(nrow(g), 0, sigmag) + b[,sp]*lesion)) %>% exp 
    }
  }
  return(postraw)
}

#function for plotting posterior predictions against data
fplot2 <- function(sp, ...){
  #mean fit
  postraw <- fpost(sp, link=T)
  postm <- apply(postraw, 2, mean)
  postPI <- apply(postraw, 2, PI)
  spcol <- dat2$species==sp
  plot(dat2$lesion[spcol], dat2$count[spcol], pch=16, col=alpha('black', .5), xlim=c(0,1), main=spnames[sp], ...)
  lines(lesionseq, postm)
  shade(postPI, lesionseq)
  #predicted fit
  postraw <- fpost(sp, link=F)
  postm <- apply(postraw, 2, mean)
  postPI <- apply(postraw, 2, PI)
  lines(lesionseq, postm)
  shade(postPI, lesionseq)
}
par(mfrow=c(2,5))
lapply(1:10, fplot2)
reset()

########################
#try it another way.

