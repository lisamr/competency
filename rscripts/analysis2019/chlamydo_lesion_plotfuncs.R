#plotting functions for sporangia ~ lesion models
#simulate prediction data using same dataset to see fit
psims <- function(data, model) {
  #sample from posterior
  sims<- add_predicted_draws(data, model) %>% 
    group_by(species, .draw) %>% 
    sample_draws(30) 
  #plot draws against data
  sims %>% ggplot() +
    geom_density(aes(x=.prediction, group=.draw),lwd=.1, alpha=.2, color='grey')+
    stat_density(data=data, aes(x=count), geom="line", color='steelblue') +
    facet_wrap(~species, scales = 'free')
}

#coef plot
pcoef <- function(model){
  #coef plot
  coefs <- model %>% 
    spread_draws(b_lesion, r_species[species,term]) %>%
    mutate(slope=r_species + b_lesion) %>% 
    filter(term=="lesion")
  
  coeftable <- coefs %>% median_hdci(slope, .width = .9) 
  
  #plot
  PLOT <- coefs %>% 
    ggplot(aes(y = fct_rev(species), x = slope)) +
    geom_halfeyeh(.width = .9, size=.5, point_interval = median_hdcih) +
    geom_point(data = coeftable, aes(slope, species), color='steelblue', size=2) +
    geom_vline(xintercept=0, lty=2, color='grey50') +
    labs(y='coefficient') 
  return(list(coeftable, PLOT))
}

#predict model fit for each species
ppredict <- function(data, model, prow=2, pcol=5){
  #create new data to predict to and predict
  xx <- seq(0,1.3,length.out = 100)
  newd <- data.frame(lesion=xx, leafID=1, species=rep(unique(data$species), each=length(xx)))
  head(newd)
  newd <- newd %>% mutate(detached = ifelse(species %in% c('PIPO', 'PSME', 'SESE', 'LIDED', 'UMCAD'), 1, 0))
  pp <- fitted(model, newdata = newd, re_formula = ~(lesion|species), summary = F, scale = 'response')
  #pp <- pp*102.5 #convert to whole leaf sample
  
  #define plot function
  plotmodfit <- function(i, pp, ...){
    sp <- as.vector(unique(newd$species))
    use <- as.integer(factor(newd$species, levels = sp))==i
    plot(data$lesion[data$species==sp[i]], data$count[data$species==sp[i]], ...,
         xlab='Lesion area cm2', ylab='# sporangia', main=sp[i], 
         pch=16, col=alpha('slateblue', .5), xlim=c(0,1.3))
    apply(pp[,use], 2, median) %>% lines(xx, .)
    apply(pp[,use], 2, HPDI, .9) %>% shade(., xx)
  }
  #plot it
  nspecies <- length(unique(data$species))
  par(mfrow=c(prow,pcol))
  sapply(1:nspecies, function(x) plotmodfit(x, pp))
  reset()
}
