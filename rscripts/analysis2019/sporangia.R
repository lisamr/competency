setwd("~/Box/Competency project/competency.git/data2019")
rm(list = ls())#clear environment
library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)
library(scales)

#read in master file
df <- read.csv(file = 'master_tall.csv')
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

#there's a fair amount of zeros. I wonder if it's warranted to do a zero inflated model. anyways, viz omitting <1.
df1 %>% filter(spore_assay=="S", !countm<1, trt=="T") %>% 
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
m1 <- glmer(count ~ species + (1|leafID), data = dtall, family = "poisson")
summary(m1)
m2 <- glmer.nb(count ~ species + (1|leafID), data = dtall)
summary(m2)

########################################
#BAYESIAN STATS?
#goal: create a mixture model (zero-inflated) with a poisson likelihood and leafID as random intercept
library(rethinking)

#troubleshoot by using only 1 count. not mixed model.
df3 <- df2 %>% select(species, count1) %>% filter(!is.na(count1)) %>% rename(count=count1)
df3 <- df3[-26,] #for now remove acma outlier

#analyze data with stan
dat <- list(
  count=df3$count,
  species=as.integer(factor(df3$species))#species to numeric index. make sure index is consecutive.
)

#model it
m3 <- ulam(
  alist(
    count ~ dpois(lambda),
    log(lambda) <- a[species], #likelihood
    a[species] ~ dnorm(3, 1.5) #fixed effect
  ), data=dat, chains=3
)

traceplot(m3) #looks good
dev.off()

#check out results
spnames <- levels(factor(df3$species))
precis(m3, depth = 2) %>% plot

#check out posterior means 
post <- extract.samples(m3)
post <- post$a
post2 <- exp(post)
postm <- apply(post2, 2, mean)

#check out how well the model fits the data. predict new data points from model and contrast with real data.
postpred <- sim(m3, data = data.frame(count= postm, species=1:10))

preddf <- as.data.frame(postpred) 
names(preddf) <- spnames
preddf <- melt(preddf)
names(preddf) <- c("species", "count")

#fit doesn't look that great. let's use a better model (and more data)
preddf2 <- rbind.data.frame(cbind(preddf, type="pred"), cbind(df3, type="obs"))
ggplot(preddf2, aes(count, group=type, fill=type))+
  geom_density(alpha=.5)+
  facet_wrap(~species, scales = 'free')
###################################
#mixed model with all 3 counts. leafID is a random effect

#analyze data with stan
#make sure index variables are consecutive integers
dat2 <- list(
  species=as.integer(factor(dtall$species)),
  leafID=as.integer(factor(dtall$leafID)),
  count=dtall$count
)

#figure out prior for the species intercept
curve( dlnorm( x, 3, 1.5) , from=0 , to=100 , n=200, xlab = "mean # spores (lambda)")
curve( dlnorm( x, -2, 1.5) , from=0 , to=50 , n=200, xlab = "deviation around the mean")

#model it
m4 <- ulam(
  alist(
    count ~ dpois(lambda), #likelihood
    log(lambda) <- a[species] + b[leafID], 
    #priors
    a[species] ~ dnorm(3, 1.5), 
    #adaptive priors
    b[leafID] ~ dnorm(b_bar, sigma_b),
    b_bar ~ dnorm(0, 1.5),
    sigma_b ~ dexp(1)
  ), data=dat2, chains=3
)
traceplot(m4) #the alphas dont look amazing
par(mfrow=c(1,1))
stancode(m4)#check out stan code

precis(m4, depth = 2, pars = c("a", "b_bar", "sigma_b")) #it did not sample the species well at all! I don't think my model is right. 

#check out posterior means 
ex <- extract.samples(m4)
beta <- exp(ex$b)
apply(beta, 2, mean) #the error between samples?? doesn't seem like much at all.

#check out how well the model fits the data. predict new data points from model and contrast with real data.
pred <- sim(m4) #951 samples (cols) simulated 1500 times (rows)
#all of the samples
dens(pred[1,], col=alpha('black', .1))
lapply(1:500, function(x)dens(pred[x,], add = T, col=alpha('black', .1)))
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
spp <- unique(colnames(pred))
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

#understand what the model is saying...
precis(m4, depth = 2, pars = c("a", "b_bar", "sigma_b")) %>% plot(labels=c(spp, "b_bar", "sigma_b"))

####################
#get back transformed mean and sd values of counts
####################

#contrasts between species
ex$a %>% head #posterior of species effects
post_a <- ex$a
pairs <- combn(1:10,2)
f <- function(x){
  spdiff <-  post_a[,pairs[1,x]]-post_a[,pairs[2,x]]
  HPDI(spdiff, .95)
}
diffs <- sapply(1:ncol(pairs), f)
#put contrasts into a readable dataframe
diffsdf <- rbind(pairs, diffs) %>% t %>% as.data.frame
diffsdf <- diffsdf %>% 
  mutate(sp1=spp[diffsdf[,1]],
         sp2=spp[diffsdf[,2]],
         sig=sign(diffsdf[,3])==sign(diffsdf[,4]))
diffsdf
write.csv(diffsdf, row.names = F, "~/Box/Competency project/competency.git/output/sporangiacontrasts.csv")

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
  facet_wrap(~species, scales = "free_x")

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

