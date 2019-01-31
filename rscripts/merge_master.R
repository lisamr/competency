#merging data for spore counts, lesion area and notes into one master file for competency tests done in June 2018
setwd("~/Desktop/Competency project/")

library(dplyr)
library(ggplot2)

#load data
master <- read.csv("data/attributes.csv", header = T)
lesion <- read.csv("data/lesions.csv", header = T)
spores <- read.csv("data/sporecounts.csv", header = T)
head(master)
head(lesion)
head(spores)

#MERGE DATAFRAMES
#clean up the data, making sure they all have the same # rows
dim(master);dim(lesion);dim(spores) #they don't all have the same number of rows...
#finds which rows in column 1 are not in column 2
fmissing <- function(col1, col2){
  r <- which(col1 %in% col2 == F)
  col1[r]
}
fmissing(lesion$code, spores$code) #has one empty. delete it. 
empty <- which(lesion$code=="")
lesion <- lesion[-empty,]

names(master)
names(spores)
names(lesion)
lesion2 <- select(lesion, code, p.area, manual, Notes_lesions)
#Ok. now merge. can ignore factor length warning. has to do wtih the empty.
master2 <- master %>% 
  left_join(spores, by=c("No", "code", "spore")) %>%
  left_join(lesion2, by=c("code")) %>% 
  mutate(code = as.factor(code),
         notes=paste(notes, Notes_spores, Notes_lesions, sep = ""),
         Notes_spores=NULL,
         Notes_lesions=NULL) %>% 
  select(No:tubes, count.mu:manual, notes)
head(master2)

#things that need to be changed to NA: p.area for VAOV, ch. counts for QUPA 1-3, QUCH 1-3, QUAG 1- 3
master2 <- master2 %>% 
  mutate(p.area = ifelse(species=="VAOV", NA, p.area),
         count.mu = ifelse(species %in% c("QUPA", "QUCH", "QUAG") & spore == "C", NA, count.mu)) 

write.csv(master2, "data/MASTERMERGED.csv", row.names = F)

#############################################
#############################################
#I want a wide version of the master sheet so that I can link spore counts and lesion areas to the same leaf. remember, for most species, 4 leaf discs were cut from the same leaf and each disc had lesion area calculated.
master2 <- read.csv("data/MASTERMERGED.csv")
head(master2)
#just looking at leaf disc experiment
master2 <- filter(master2, assay=="L") %>% select(-notes)

#make wide version of master
tmp1 <- master2 %>%  dcast(., species + ind + trt + leaf ~ spore, var.name = "test", value.var = c("count.mu")) %>% rename(countC=C, countS=S)
tmp2 <-master2 %>%  dcast(., species + ind + trt + leaf ~ spore, var.name = "test", value.var = c("p.area")) %>% rename(p.areaC=C, p.areaS=S)

master3 <- tmp1 %>% 
  left_join(tmp2, by=c("species", "ind", "trt", "leaf")) %>%
  mutate(countS = countS*100, #get absolute # by multiplying concentration by volume of water/LPCB added to leaf disc 
         areaC = p.areaC * 1.227, #convert proportion area to absolute area
         areaS = p.areaS * 1.227,
         avg.area = (areaC + areaS)/2,
         C.per.area = countC/areaC, # #spores per area
         S.per.area = countS/areaS,
         SC_ratio = countS/countC) 
head(master3)
#avg lesion area couldn't be calculated with ones with NA. can't do it with dplyr for some reason.
master3$avg.area <- ifelse(is.na(master3$areaC), master3$areaS, master3$areaC)
#add in water column from master sheet
master4 <- master2 %>% filter(spore=="S") %>% select(species, ind, trt, leaf, water) %>% right_join(master3)
dim(master3); dim(master4)


write.csv(master4, "data/MASTERMERGED_wide.csv", row.names = F)

master4 %>% 
  group_by(species, ind, trt) %>% 
  summarise(n()) %>% View






