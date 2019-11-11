#script to collate all the important values into tables that will go into the manuscript
library(dplyr)
library(ggplot2)
setwd("~/Box/Competency project/competency.git")

#get lesion size summarie
df <- read.csv(file = 'data2019/master_tall.csv')
#analyze just trt, sporangia
lesion_sum <- df %>% 
  group_by(species) %>% 
  filter(trt=="T", is.na(omit)) %>% 
  summarise(mean_les=mean(perc_lesion, na.rm = T), sd_les=sd(perc_lesion, na.rm = T)) %>% 
  mutate_if(is.numeric, signif, 3)
lesion_sum

#species intercepts (mean, sd, lower, upper)
spor <- read.csv('output/sporangia_both/predictions.csv')
chlam <- read.csv('output/chlamydo_both/chlamydosbroad_predicted_both.csv')

left_join(spor, chlam, by='species') %>% 
  left_join(lesion_sum, by='species') %>% 
  select(species, estimate.x, sd.x, lower.x, upper.x,
         estimate.y,sd.y, lower.y, upper.y,
         mean_les, sd_les) %>% 
  write.csv('output/final/meansummaries.csv', row.names = F)
