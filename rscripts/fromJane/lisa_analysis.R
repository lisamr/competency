# Analysis of Lisa Rosenthal's spore data

# TO DO:
# 
# simulate more data from the model and re-fit. See what is recovered. 
# sensitivity analysis
# plot sampled data from the bayesian fit against true data

# Questions:
# Q1 - Do species have an effect? 
# Q2 - How to handle individual variation (repeated measures)
# Q3 - Model to describe the data?
# Q4 - Bayesian analysis?
# Q5 - Power analysis / how much data we need to rank or group the species?

##############################################################################################################


# 0. Setup ----

library(dplyr)
library(rstan)
library(tidyr)
library(forcats)
library(lme4)
library(plotly)


wide = read.csv("~/Desktop/Competency project/data/MASTERMERGED_wide.csv")

wide.treatment <- filter(wide, trt=="T")

#    Modificiations to data for preliminary analysis ----
# standardize repetitions
# remove species with only one individual
# replace na's with median value
wide.treatment = filter(wide.treatment, species !="CEOL") # CEOL has only one individual -- removing it for preliminary
wide.treatment = droplevels(wide.treatment)
wide.treatment %>% group_by(species) %>% summarize(n_distinct(ind)) # all others have three individuals
# most individuals have six observations but a couple 12 -- making them all 6 for preliminary:
wide.treatment %>% group_by(species, ind) %>% summarize(n())
wide.treatment = wide.treatment %>% group_by(species, ind) %>% 
  arrange(is.na(countS)) %>% #get rid of NAs first if need to get rid of anything
  mutate(index = 1:n(), countS.ind.mean = mean(countS, na.rm = T)) %>% filter(index <=6) #cuts 18
wide.treatment = wide.treatment %>% arrange(species, ind) %>% 
  mutate(ind2 = ind+3*as.numeric(species)) %>% ungroup()
         
# replace NAs with median if necessary:
# wide.treatment = wide.treatment %>% group_by(species, ind) %>% mutate(counts = replace_na(countS, median(countS, na.rm = T)))

# 1. EDA ----

# boxplots for individuals group by species also showing data points
ggplot(wide.treatment, aes(y = countS, x = species, group = interaction(species,ind))) +
         geom_boxplot() + theme_bw() +
         geom_point(aes(y = countS, group = interaction(species,ind), color = as.factor(ind), shape = as.factor(ind)))
ggplotly(.Last.value) #enables zoom in for species near zero

# point distributions for individuals, colored by species
ggplot(wide.treatment, aes(y = countS, x = ind2)) + geom_point(aes(color = species)) + theme_bw()


# 2. Anova/GLM 
#    One-way ANOVA (bad) (Q1) ----
# Use individual means - deals with non-independence of repeated measures of individuals.
# Other assumptions of ANOVA (normality of residuals, homoskedasticity, aren't met)
# We shouldn't pay attention to it.
lm1 = lm(countS.ind.mean ~ -1 + species, data = wide.treatment) # removed intercept for easily interpretable estiamtes
summary(lm1)
aov1 = aov(lm1); aov1
plot(aov1)

#    Anova with one random effect and one fixed effect (Q2) ----
# http://www.maths.bath.ac.uk/~jjf23/mixchange/rbd.html
# https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
# Extending the linear model with R (faraway) chapter 8 and 9.2

# We have repeated measures of each individual, 
  # Species is fixed and individual is a random effect (not a level/type chosen by the researcher)
  # (Standard 2-way anova would assume two fixed effects)

  # lisa.mmod <- glmer(countS ~ species + (1|ind), wide.treatment, family = "poisson") bad, groups individuals into just 1,2,3 groups
  # random effect of individual nested in species
  # Assumes constant variance for random effects

  # Normal - note range of response var definitely not normal
  lisa.mmod <- lmer(countS ~ -1 + species + (1|species:ind), wide.treatment) 
  summary(lisa.mmod)
  anova(lisa.mmod) # can check the species effect with an anova
  lisa.mmod.ci = data.frame(confint(lisa.mmod)[-c(1,2),]); names(lisa.mmod.ci) = c("lower_2.5","upper_97.5")
  lisa.mmod.ci$effect = as.factor(rownames(lisa.mmod.ci))
  ggplot(lisa.mmod.ci) + geom_segment(aes(x = lower_2.5, xend = upper_97.5, y = effect, yend = effect)) + 
                         geom_point(data = data.frame(summary(lisa.mmod)$coefficients),
                                    aes(x = Estimate, y = 1:12)) + theme_bw()
  
  # Poisson - for counts, but not enough variance
  lisa.gmmod <- glmer(countS ~ species + (1|species:ind), wide.treatment, family = "poisson") 
  summary(lisa.gmmod)
  #confint(lisa.gmmod) fails?

  # Neg.binom
  lisa.gmmod.nb <- glmer.nb(countS ~ species + (1|species:ind), wide.treatment) 
  #glmer(countS ~ species + (1|ind2), wide.treatment, family = "poisson") 
  summary(lisa.gmmod.nb)
  lisa.gmmod.nb.ci = confint(lisa.gmmod.nb)

# 3. Poisson lack of fit (Q3) ----
  
  #These distributions are way more dispersed.
  
  # lefts examine points less likely than the following cutoff using group means as poisson parameters
  cutoff1 = 1.0e-10
  # for individual distributions
  wide.treatment %>% group_by(species, ind) %>% 
    mutate(
      ind.mean = round(mean(countS)),
      dpois = dpois(countS, lambda = ind.mean),
      major_outlier = dpois < cutoff1) %>% 
    select(ind.mean, countS, dpois, major_outlier) %>%
    ggplot(aes(x = fct_reorder(interaction(species, ind), as.numeric(species)), y = countS, 
               color = species, size = major_outlier)) + 
    geom_point(alpha = .5) + theme(axis.text.x = element_text(angle = 90))
  
  
  # for species distributions
  wide.treatment %>% group_by(species) %>%
    mutate(species.mean = mean(countS)) %>%
    group_by(species, ind) %>%
    summarize(ind.mean = round(mean(countS)),
              species.mean = round(first(species.mean)),
              dpois = dpois(round(ind.mean), lambda = species.mean),
              major_outlier = dpois < cutoff1
    ) %>% 
    select(species.mean, ind.mean, dpois, major_outlier) %>%
    ggplot(aes(x = species, y = ind.mean, color = species, size = major_outlier)) + 
    geom_point(alpha = .5) +  theme(axis.text.x = element_text(angle = 90))  

# 4. Should we use negative binomial (uses bayesian output below)? (Q4) ----
  convert_to_nbinom <- function(lambda, phi) {
    prob = phi/(1+lambda)
    size = lambda*phi/(1+lambda - phi)
    return(list(size = size, prob = prob))
  }
  
  # individual
  lambda = lisa_summary_individual[grepl("lambda", rownames(lisa_summary_individual)), 1]
  phi = lisa_summary_individual[grepl("phi", rownames(lisa_summary_individual)), 1]
  sp = convert_to_nbinom(as.vector(lambda), phi)
  with(sp, {assign("mean.ind", size*(1-prob)/prob, pos = 1)})
  with(sp, {assign("var.ind", size*(1-prob)/prob^2, pos = 1)})
  
  dnbinom(wide.treatment$countS, size = sp$size, p = sp$prob) # Better, but huge variances
  sum(.Last.value < cutoff1)
  
  # species
  lambda = lisa_summary_species[grepl("lambda", rownames(lisa_summary_species)), 1]
  phi = lisa_summary_species[grepl("phi", rownames(lisa_summary_species)), 1]
  sp = convert_to_nbinom(as.vector(lambda), phi)
  with(sp, {assign("mean.spec", size*(1-prob)/prob, pos = 1)})
  with(sp, {assign("var.spec", size*(1-prob)/prob^2, pos = 1)})
  
  dnbinom(wide.treatment$countS, size = sp$size, p = sp$prob) # Better, fewer huge variances and means are closer in
  sum(.Last.value < cutoff1)
  
  # compare
  par(mfrow= c(1,1))
  plot(mean.ind, type = "l"); points(rep(mean.spec, each = 6), type = "l", col = "red"); points(wide.treatment$countS)
  
# 5. Bayesian (Q4) ####

Sys.setenv(USE_CXX14 = 1)
rstan_options(auto_write = T)

#    model with individual effects drawn from species level means ----

# data
lisa_data = list(
  S = n_distinct(wide.treatment$species),
  I = 3,
  J = 6,
  Y = wide.treatment$countS # data arranged by species then individual 
)

# fit
lisa_model = "~/Documents/DSI/consultations/lisa_model.stan"
stanc(lisa_model)
lisa_fit = stan(file = lisa_model, model_name = "lisa", data = lisa_data, 
                 control = list(max_treedepth = 15), verbose = TRUE, iter = 1000, chains = 2)

# results
print(lisa_fit)
plot(lisa_fit, pars = c("alpha_s", "alpha0")) # means
plot(lisa_fit, pars = "beta") # individual
plot(lisa_fit, pars = c("sigma_s", "sigma0", "phi")) # dispersion
lisa_summary = summary(lisa_fit, pars = c("alpha_s", "phi", "sigma_s", "beta","lambda"), probs = c(0.1, 0.9))$summary
lisa_extracted = rstan::extract(lisa_fit)

# negative binomial fit
lambda = lisa_summary[grepl("lambda", rownames(lisa_summary)), 1]
phi = lisa_summary[grepl("phi", rownames(lisa_summary)), 1]
sp = convert_to_nbinom(as.vector(lambda), phi)
with(sp, {assign("mean.lisa", size*(1-prob)/prob, pos = 1)})  # should equal lambda
with(sp, {assign("var.lisa", size*(1-prob)/prob^2, pos = 1)}) # some very large variances

dnbinom(wide.treatment$countS, size = sp$size, p = sp$prob) # Better, fewer huge variances and means are closer in
sum(.Last.value < cutoff1)

# compare to poisson
lisa_summary = summary(lisa_fit_individual_poisson, pars = c("alpha", "beta","lambda"), probs = c(0.1, 0.9))$summary
lambda = lisa_summary[grepl("lambda", rownames(lisa_summary)), 1]
plot(wide.treatment$countS)
points(mean.lisa, type = "l", col = "black", lwd = 2)
points(lambda, type = "l", col = "green", lwd = 2)
points(rep(1:216, each = 1000), c(lisa_extracted[["lambda"]]), pch = ".")

#    (preliminary) model with individual effects ----

# data
lisa_data_individual = list(
  S = n_distinct(wide.treatment$species),
  I = 3,
  J = 6,
  Y = wide.treatment$countS # data arranged by species then individual 
)

# fit
lisa_model_individual = "~/Documents/DSI/consultations/lisa_individual.stan"
stanc(lisa_model_individual)
# lisa_fit_individual_poisson -- stored result before switching to neg binom
lisa_fit_individual = stan(file = lisa_model_individual, model_name = "lisa_individual", 
                           data = lisa_data_individual, 
                           control = list(max_treedepth = 15), verbose = TRUE, iter = 1000, chains = 2)

# results
print(lisa_fit_individual)
#plot(lisa_fit_individual_poisson, pars = "alpha") # species
plot(lisa_fit_individual, pars = "alpha") # species
#plot(lisa_fit_individual_poisson, pars = "beta") # species
plot(lisa_fit_individual, pars = "beta") # individual
plot(lisa_fit_individual, pars = "phi") # individual
lisa_summary_individual = summary(lisa_fit_individual, pars = c("phi", "alpha","beta","lambda"), probs = c(0.1, 0.9))$summary

lisa_extracted_individual = rstan::extract(lisa_fit_individual)


#    (preliminary) model with only species effects ----

# data
lisa_data_species = list(
  S = n_distinct(wide.treatment$species),
  I = 3,
  Y = round((wide.treatment %>% group_by(species, ind) %>% summarize(mean = mean(countS)))$mean)# data arranged by species then individual 
)

lisa_model_species = "~/Documents/DSI/consultations/lisa_species.stan"
stanc(lisa_model_species)
lisa_fit_species= stan(file = lisa_model_species, model_name = "lisa_species", data = lisa_data_species, 
                           control = list(max_treedepth = 15), verbose = TRUE, iter = 1000, chains = 2)

# results
print(lisa_fit_species)
plot(lisa_fit_species, pars = c("alpha","alpha0")) # species
lisa_summary_species = summary(lisa_fit_species, pars = c("phi", "alpha", "lambda"), probs = c(0.1, 0.9))$summary
lisa_extracted_species = rstan::extract(lisa_fit_species)

#    (bonus) investigate poisson lack of fit using stan output ----
lisa_summary_individual_poisson = summary(lisa_fit_individual_poisson, pars = c("lambda"))$summary
lambda = lisa_summary_individual_poisson[grepl("lambda", rownames(lisa_summary_individual_poisson)), 1]
dpois(wide.treatment$countS, lambda) #not so good
sum(.Last.value < cutoff1)

lisa_summary_species_poisson = summary(lisa_fit_individual_poisson, pars = c("lambda"))$summary
lambda = lisa_summary_individual_poisson[grepl("lambda", rownames(lisa_summary_individual_poisson)), 1]
dpois(wide.treatment$countS, lambda) #not so good
sum(.Last.value < cutoff1)


# 6. Power using simulated data ----

# How would the Bayesian results look with different levels of data and spread?
# How much data to be need to recover true parameters?
