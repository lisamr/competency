#goal: run through some mock spore data to analyze with stan. 
#approach: do with glm first to get an idea of what to expect and see if that agrees with the stan model. 

library(dplyr)
library(ggplot2)
library(rethinking)

#simulate spore data
n_spp <- 10 #number of spp
L <- rgamma(n_spp, shape = .5, scale = 20) #distribution of lambda across species 
hist(L);L
n_ind <- 32 #number of individuals per species
spores <- rpois(n_ind*n_spp, L) #assume poisson for starters (nbinom later)

#visualize data
dat <- data.frame(species=rep(1:n_spp, n_ind), species2=rep(letters[1:n_spp], n_ind), lam=L, ind=rep(1:n_ind, each=n_spp), spores) #assemble data frame
ggplot(dat, aes(species2, spores))+
  geom_boxplot()

#analyze data with glm
head(dat)
m1 <- glm(spores~species2-1, data=dat, family = "poisson")
summary(m1)
plot(m1)
exp(coef(m1))
log(L)

#analyze data with stan
dat2 <- list(
  S=dat$spores,
  plant=dat$species
)

#figure out prior for the species intercept. mean number of spores is L
L
curve( dlnorm( x, 4, 1.5) , from=0 , to=100 , n=200, xlab = "mean # spores (lambda)")

#model it
m2 <- ulam(
  alist(
    S ~ dpois(lambda),
    log(lambda) <- a[plant],
    a[plant] ~ dnorm(4, 1.5)
    ), data=dat2, chains=3
  )

#check out results
precis(m2, depth = 2) %>% plot

#check out posterior
post <- extract.samples(m2)
post <- post$a

post2 <- exp(post)
postm <- apply(post2, 2, mean)
posth <- apply(post2, 2, HPDI, .95)
rbind(L, posth)

#how does missing values change analysis? (unbalanced design)
datm <- dat %>%
  filter(!(species %in% c(1,2) & ind %in% (1:20)))
datm %>% group_by(species) %>% tally
#analyze data with stan
dat3 <- list(
  S=datm$spores,
  plant=datm$species
)
#model it
m3 <- ulam(
  alist(
    S ~ dpois(lambda),
    log(lambda) <- a[plant],
    a[plant] ~ dnorm(4, 1.5)
  ), data=dat3, chains=3
)
#traceplot(m3)
#dev.off()

precis(m3, depth = 2)
precis(m2, depth = 2)
plot(coeftab(m3, m2))

#pairwise contrasts?

####################################################
####################################################
#negative binomial data

#simulate spore data
n_spp <- 10 #number of spp
L <- rgamma(n_spp, shape = .5, scale = 20) #distribution of lambda across species 
hist(L);L
n_ind <- 32 #number of individuals per species
spores2 <- rgampois(n_ind*n_spp, L, 2) #mu=lambda, scale=phi(dispersion of means)
spores2

#visualize data
dat_nbin <- data.frame(species=rep(1:n_spp, n_ind), species2=rep(letters[1:n_spp], n_ind), lam=L, ind=rep(1:n_ind, each=n_spp), spores2) #assemble data frame
ggplot(dat_nbin, aes(species2, spores))+
  geom_boxplot()
ggplot(dat_nbin, aes(species2, spores))+
  geom_violin()

