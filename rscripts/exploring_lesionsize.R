#Is any of the lesion area data useful? I'm wondering if I should bother taking pictures of the leaves for the redo experiment.

#sporangia analysis.
library(dplyr)
library(ggplot2)


wide <- read.csv("data/MASTERMERGED_wide.csv")
head(wide)
#clean up data. remove CEOL (only 1 ind)
wideT <- filter(wide, trt=="T", !is.na(countS), species!="CEOL") %>% mutate(ind2=interaction(species, ind)) %>% droplevels()
wideT %>% group_by(species, ind) %>% summarise(min(countS), max(countS)) 

#Q1: do plants still sporulate when the are asymptomatic?
#Q2: does lesion size predict sporulation in general?
#Q3: does lesion size predict sporulation within a given species?

#Q1
#calculate average sporulation for leaves that had no lesions
wideT %>% filter(p.areaS==0) # QUAG, QUCH, QUPA only. 
wideT %>% filter(p.areaS==0) %>% summarise(muS=mean(countS), sdS=sd(countS))
wideT %>% filter(p.areaC==0) #most were the oaks and couldn't count chlamydospores.
#consensus: just to shows with data that leaves produce spores without lesions is enough to keep it. But I can do that with the sporangia viability assays.

#Q2
#check out sporangia vs lesion size. when log transforming sporangia, theres an upward trend with wide interval. Still plenty of plants with no lesions and sporulation.
ggplot(wideT, aes(p.areaS, (countS))) +
  geom_point()
ggplot(wideT, aes(p.areaS, log(countS))) +
  geom_point() 

#Q3. It looks like lesion area correlates with sporangia for some species but not for others. I think I should keep the lesion size just so I can address this part of the story.
ggplot(wideT, aes(p.areaS, (countS), color=species)) +
  geom_point() +
  facet_wrap(~species)
ggplot(wideT, aes(p.areaS, log(countS), color=species)) +
  geom_point() +
  facet_wrap(~species)


