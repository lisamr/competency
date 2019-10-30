rm(list=ls()) #clean env

#simulate generating process of data
#params
abar <- 1.5 # average spore intercept
bbar <- .1 # average lesion slope
sigma_a <- .15 # std dev in intercepts
sigma_b <- 0.5 # std dev in slopes
rho <- .6 # correlation between intercepts and slopes

#set correlations between the int and slope
Mu <- c(abar, bbar)
sigmas <- c(sigma_a,sigma_b) # standard deviations 14.5
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
# now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

#generate int and slope for each species
library(MASS)
set.seed(5) # used to replicate example
N_species <- 10
vary_effects <- mvrnorm( N_species , Mu , Sigma )
a <- vary_effects[,1] #species specific int
b <- vary_effects[,2] #species specific slopes

#simulate observations
set.seed(22) 
N_ind <- 32 #n reps or individuals per species
sp_ID <- rep( 1:N_species , each=N_ind )
lesion <- rbeta(length(sp_ID), 1,5) #random lesion size bounded between 0 and 1
lambda <- a[sp_ID] + b[sp_ID]*lesion
spores <- rpois(length(sp_ID), exp(lambda))#generate response
d <- data.frame( species=sp_ID , lesion=lesion , count=spores)
head(d)

#plot data
ggplot(d, aes(lesion, count))+
  geom_point()+
  facet_wrap(~species, scales = "free_y")

#run model
#set prior first
f1 <- bf(count ~ lesion + (lesion|species), family = poisson())
get_prior(f1, d)
prior1 <- c(
  set_prior("normal(2, 1.5)", class="Intercept"),
  set_prior("normal(0, .5)", class="b"),
  set_prior("exponential(1)", class = "sd") )

m1 <- brm(
  f1, data=d, prior = prior1,
  chains = 3, cores = 4, iter=3000,
  control = list(adapt_delta = .95, max_treedepth=15))

#check out model
summary(m1) #close. intercept is off a bit.
parnames(m1)
#population level effects only. traceplots look good in both models.
plot(m1, pars = "^b_")

#simulate prediction data using same dataset to see fit
sims<- add_predicted_draws(d, m1) %>% 
  group_by(species, .draw) %>% 
  sample_draws(50) 
sims %>% ggplot() +
  geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.2, color='grey')+
  stat_density(data=d, aes(x=count), geom="line", color='steelblue') +
  facet_wrap(~species, scales = 'free')

#coef plot
coefs <- m1 %>% 
  spread_draws(b_lesion, r_species[species,lesion]) %>%
  mutate(slope=r_species + b_lesion) %>%
  rename(par=lesion, value=slope)
#point estimates for the random slopes
rand_slopes <- coef(m1)$species[,,2]
rand_slopes[,1]

#plots
coefs %>% filter(par=='lesion') %>% 
  ggplot(aes(y = interaction(par, species), x = value)) +
  geom_halfeyeh(.width = .9, size=.5) +
  geom_vline(xintercept=0, lty=2, color='grey50') +
  labs(y='coefficient') +
  geom_point(data=data.frame(species=1:10, b), aes(b, species), color='steelblue')#'real' slopes

#predict model fit for each species
xx <- seq(0,1,length.out = 100)
newd <- data.frame(lesion=xx, species=rep(1:10, each=length(xx)))
pp1 <- fitted(m1, newdata = newd,scale = 'response', re_formula = ~(lesion|species), summary = F)
pp1.1 <- fitted(m1, newdata = newd,scale = 'response', re_formula = ~(1|species), summary = F)#wrong syntax!!!!
pp2 <- predict(m1, newdata = newd,scale = 'response', re_formula = ~(lesion|species), summary = F)

plotmodfit <- function(i, pp, ...){
  sp <- 1:10
  use <- newd$species==i
  beta.est <- rand_slopes[,1]
  plot(d$lesion[d$species==sp[i]], d$count[d$species==sp[i]], ...,
       xlab='lesion', ylab='# sporangia', 
       main=paste(sp[i], round(beta.est[i], 2)), 
       pch=16, col=alpha('steelblue', .5), xlim=c(0,1))
  apply(pp[,use], 2, median) %>% lines(xx, .)
  apply(pp[,use], 2, PI, .9) %>% shade(., xx)
}
#looks right
par(mfrow=c(2,5))
sapply(1:10, function(x) plotmodfit(x, pp1))

