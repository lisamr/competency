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
#try it with glmer first
m0 <- glmer(count ~ perc_lesion + (1|leafID) + (perc_lesion|species), data = dtall2, family = 'poisson')
summary(m0)

#thats cool it worked. Not sure how you check out model performance.
plot(m0)
plot(dtall2$count, predict(m0, type='response'))+abline(0,1)

#how do you do counterfactual plots again with lmer obj??
#unique(dtall2$species)
#pd <- data.frame(perc_lesion=1:100, species=rep(unique(dtall2$species), each=100))
#pm0 <-bootMer(m0, function(x) predict(m0, re.form=NA, newdata=pd, type='response'), nsim = 100)
#str(pm0)
#plot(pd$perc_lesion, pm0$t0)
#for(i in 1:nrow(pm0$t)){
#  lines(pd$perc_lesion, pm0$t[i,])}
#plot(pm0$t[1,], type='l')

########################################
#statistical models!
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
#model it?!? having issues. maybe make simpler by removing 3 subsamples (leafID) to troubleshoot and then add it back in?
dtall3 <- dtall2 %>% filter(!duplicated(leafID))
dat <- list(
  count=dtall3$count,
  species=as.integer(factor(dtall3$species)),
  #ID=dtall3$leafID,
  lesion=dtall3$perc_lesion)

#no leafID
m1 <- ulam(
  alist(
    #likelihood and linear model
    count ~ dpois(lambda),
    log(lambda) <- a[species] + b[species]*lesion,
    #varying effects
    c(a, b)[species] ~ multi_normal(c(abar, bbar), rho, sigma),
    #fixed priors
    abar ~ dnorm(3, 1.5),
    bbar ~ dnorm(1, .5),
    rho ~ dlkjcorr(2),
    sigma ~ dexp(1)
    ), 
  data=dat , chains=4 , cores=4 )
traceplot(m1)
dev.off()
precis(m1, depth = 3) %>% plot #sampled well!

########################
#put leafID back in?
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
traceplot(m2) 

postcheck(m2)
reset()
parsm2 <- m2@pars
precis(m2, pars = parsm2[-1], depth = 3) %>% plot 

####################################
####################################
#check out the model!
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
