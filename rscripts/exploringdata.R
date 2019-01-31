#first stab at looking at all of the (preliminary) data
setwd("~/Desktop/Competency project/")

library(dplyr)
library(ggplot2)
library(reshape2)

master <- read.csv("data/MASTERMERGED.csv")

#only going to look at leaf disc assay. detached leaves didnt work.
master1 <- master %>% filter(assay=="L") %>% select(-notes)
tail(master1)

ggplot(filter(master1, spore=="S", trt=="T"), aes(species, count.mu)) +
  geom_point()+
  geom_boxplot()
ggplot(filter(master1, spore=="C", trt=="T"), aes(species, count.mu)) +
  geom_point()+
  geom_boxplot()
ggplot(filter(master1, trt=="T"), aes(species, p.area)) +
  geom_point()+
  geom_boxplot()

#what does the distribution of spore counts and lesion area look like? zero-inflated poisson for sure. lesion area looks funny.
ggplot(filter(master1, spore=="S"), aes(count.mu)) +
  geom_histogram()+
  facet_wrap(~trt)
ggplot(filter(master1, spore=="C"), aes(count.mu)) +
  geom_histogram()+
  facet_wrap(~trt)
ggplot(filter(master1, species!="ACMA"), aes(p.area)) +
  geom_histogram()+
  facet_wrap(~trt)

#Let's look at data from master_wide
wide <- read.csv("data/MASTERMERGED_wide.csv")
head(wide)
#first some distributions
ggplot(filter(wide, trt=="T"), aes(log(1+countS))) +
  geom_density()+
  facet_wrap(~species)
ggplot(filter(wide, trt=="T"), aes(log(SC_ratio+1))) +
  geom_histogram()

#now relationships.
#SPORANGIA
ggplot(filter(wide, trt=="T"), aes(species, countS)) +
  geom_boxplot() 
ggplot(filter(wide, trt=="T"), aes(species, countS, group=interaction(species, ind))) +
  geom_boxplot() 
#water doesn't seem to have much effect on sporangia counts besides for ACMA, where it def increased the counts. could remove those ones for analysis.
ggplot(filter(wide, trt=="T"), aes(species, countS, group=interaction(species, ind), color=water)) +
  geom_boxplot() +
  geom_point(aes(group=ind), position=position_dodge(.8))
ggplot(filter(wide, trt=="T"), aes(species, S.per.area, group=interaction(species, ind))) +
  geom_boxplot() 
#CHLAMYDOSPORES
ggplot(filter(wide, trt=="T"), aes(species, countC)) +
  geom_boxplot() 
ggplot(filter(wide, trt=="T"), aes(species, countC, group=interaction(species, ind))) +
  geom_boxplot() 
ggplot(filter(wide, trt=="T"), aes(species, C.per.area, group=interaction(species, ind))) +
  geom_boxplot() 
#CHLAMYDO VS SPORANGIA
ggplot(wide, aes(log(countS), log(countC), color=as.factor(species))) +
  geom_point()
ggplot(wide, aes((countS), (countC), color=as.factor(ind))) +
  geom_point()+
  facet_wrap(~species)
ggplot(filter(wide, trt=="T"), aes(species, SC_ratio, group=interaction(species, ind))) +
  geom_boxplot() 
ggplot(filter(wide, trt=="T", species!="RHCA", !(species=="ACMA" & ind==3)), aes(species, SC_ratio, group=interaction(species, ind))) +
  geom_boxplot() 
#LESION AREA
ggplot(filter(wide, trt=="T"), aes(species, avg.area)) +
  geom_boxplot() 
ggplot(filter(wide, trt=="T"), aes(species, avg.area, group=interaction(species, ind))) +
  geom_boxplot() 
#SPORES VS LESIONS. generally the more spores you have, the bigger the lesion.
ggplot(wide, aes(log(1+countS), (areaS), color=as.factor(ind))) +
  geom_point()+
  facet_wrap(~species)
ggplot(wide, aes(log(1+countS), areaS,  color=as.factor(species))) +
  geom_point()+
  geom_smooth(aes(), method = "lm")
ggplot(wide, aes((1+countC), (areaS), color=as.factor(ind))) +
  geom_point()+
  facet_wrap(~species)
ggplot(wide, aes((countC), areaC,  color=as.factor(species))) +
  geom_point()+
  geom_smooth(aes(), method = "lm")


#ACMA3 is a wierd sample. I wonder if the sporangia that I was counting were really PR. look at the notes and also look at the controls. Also, where did you collect it? ACMA and QUKE had samples collected in both davis and berkeley. data isnt entered. still in my notebook in the lab.

#well the controls look fine.
ggplot(filter(wide, species=="ACMA"), 
       aes(trt, countS, group=interaction(trt,ind))) +
  geom_boxplot() 
#notes dont help. controls didnt grow on parp and treatments did. ind 3 might have just been a freak. can do a power analysis with and without the outlier.
master %>% filter(species=="ACMA", ind==3)



