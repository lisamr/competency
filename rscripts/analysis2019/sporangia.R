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

#I get convergence issues and I forget how to read these model outputs
#m1 <- glmer(count ~ species + (1|leafID), data = dtall, family = "poisson")
#summary(m1)
#m2 <- glmer.nb(count ~ species + (1|leafID), data = dtall)
#summary(m2)

########################################
#BAYESIAN STATS?
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
  ), data=dat2, chains=3
)
traceplot(m4) #the alphas dont look amazing
par(mfrow=c(1,1))
stancode(m4)#check out stan code
precis(m4, depth = 2, pars=c('a', 'sigma_b')) #%>% plot

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

