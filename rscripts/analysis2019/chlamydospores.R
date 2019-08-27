#contrast chlamydospore counts, similar to what you did for sporangia. will have more zeros and might need to do a zeroinflated model?

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
#map the plants
df1 %>% filter(spore_assay=="S", trt=="T") %>%
  ggplot(aes(easting, northing)) +
  geom_point(aes(color=species))+
  scale_color_viridis_d()
#plot chlamydospore counts on a map for select species
df1 %>% filter(spore_assay=="C", trt=="T", species=="CEOL") %>%
  ggplot(aes(easting, northing)) +
  geom_point(aes(color=count1))+
  scale_color_viridis_c()

#let's try to see how chlamydo counts differ across species
#2 choices: estimate # of spores by dividing counts by proportion sampled OR (zero inflated to handle zeros and gamma distribution?) OR use a count dist with offsets (poisson)

#start with averages
df1<- df1 %>% mutate(count_est=count1/prop) %>% filter(spore_assay=="C", !is.na(count_est), trt=="T")

#plot it. well they aren't all that different!
df1 %>% 
  ggplot( aes(species, log(count_est)))+
  geom_boxplot()+
  geom_point(alpha=.4)

#there's a fair amount of zeros. I wonder if it's warranted to do a zero inflated model. anyways, viz omitting <1.
df1 %>% filter( !count_est<1) %>% 
  ggplot( aes(species, log(count_est)))+
  geom_boxplot()+
  geom_point(alpha=.4)


#chlamydos vs lesion size?
df1 %>% 
  ggplot( aes(count_est, perc_lesion))+
  geom_point(alpha=.4)+
  facet_wrap(~species, scales = 'free_x')

