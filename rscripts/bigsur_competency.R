#how are IV (relative and total) and # sites present related to each other? I need to arrive at a metric to represent abundance or density. hopefully they all correlate.

library(dplyr)
library(ggplot2)
library(readxl)
Rename <-   function(values_to, index_from, from_column){
    result <- (values_to[match(from_column, index_from)])
    return(result)
  }

##########################################
##########################################
#SUMMARIES OF TREES BY PLOT 
df <- read_excel("/Users/lisarosenthal/Box Sync/big sur/Big Sur Github/bigsuraccessdata/LRosenthal_3_26_2018_Masterfromkerriemail.xlsx")

#data selecting, adding columns and cleaning up
df2 <- df %>% filter(SampleYear %in% c(2006, 2007)) %>% 
  select("BSPlotNumber", "ForestAllianceType", "StemNumber", "TreeNumber", "Species", "DBH", "Status", "CankerTrunkSymptoms":"TwigBasalSproutSymptoms", "Canker", "PramCulturePositive", "PramPCRPositive") %>% 
  filter(!is.na(Status), !Species=="MISSING", !is.na(DBH)) %>% 
  #get disease symptoms out. Only trusting Y. if it has symptoms and p.ram is plot level, counts as diseased (that's what Haas 2011 did)
  mutate_at(.vars = vars(c("CankerTrunkSymptoms":"TwigBasalSproutSymptoms")), .funs = function(x)(Rename(c(0,0,0, 0, 1, 1), c(NA, "N", "P", "p", "y", "Y"), x))) %>% 
  mutate(nsymptoms = rowSums(.[8:14]), P_ramorum_plotlevel = ifelse(PramCulturePositive==T|PramPCRPositive==T,1,0)) %>%
  group_by(BSPlotNumber) %>% 
  mutate(P_ramorum_plotlevel = ifelse(any(P_ramorum_plotlevel)==1,1,0),
         infected = ifelse(P_ramorum_plotlevel ==1 & nsymptoms>0,1,0 )) %>% 
  #basal area (m2), relative 
  mutate(BA = ifelse(is.na(DBH), 0, pi*(DBH/200)^2)) %>% 
  select(-c(CankerTrunkSymptoms:Canker))

#out ouf curiosity, how many trees show symptoms but no disease at plot level? 
tmp <- df2 %>% 
  group_by(BSPlotNumber) %>% 
  #filter(Status=="L") %>% 
  summarise(PR = unique(P_ramorum_plotlevel),
            #PR = ifelse(any(P_ramorum_plotlevel)>0,1,0),
            stems.w.symptoms = sum(nsymptoms>0),
            stems.w.posID = sum(PramCulturePositive==1),
            stems.infected = sum(infected==1)) 
ggplot(tmp, aes(PR, stems.w.symptoms, group=PR))+
  geom_boxplot()
ggplot(tmp, aes(PR, stems.w.posID, group=PR))+
  geom_boxplot() 
ggplot(tmp, aes(PR, stems.infected, group=PR))+
  geom_boxplot() 
tmp %>% group_by(PR) %>% summarise(nplots = n(), sum(stems.w.symptoms), sum(stems.w.posID ))

##########################################
##########################################
#get species summaries: relative IV, total IV, occurence (# plots present)
#IV = (rel.dens + rel.BA)/2
head(df2)
summary(df2)
tmp <- df2 %>% 
  #filter(Status == "L") %>% 
  group_by(BSPlotNumber, ForestAllianceType, Species) %>% 
  summarise(n.stems = n(), sum.BA = sum(BA)) %>% 
  mutate(tot.stems = sum(n.stems), total.BA = sum(sum.BA), 
         IV = (n.stems/tot.stems + sum.BA/total.BA)/2)

species.sum.ME <- tmp %>% filter(ForestAllianceType=="MixedEvergreen") %>% 
  group_by(Species) %>% 
  summarise(IV.mean = mean(IV), IV.sd = sd(IV),IV.se=sd(IV)/sqrt(n()),
            n.stems = sum(n.stems), sum.BA = sum(sum.BA),
            n.plots = n()) %>% 
  mutate(IV.tot = (n.stems/sum(n.stems) + sum.BA/sum(sum.BA))/2 )
species.sum.Red <- tmp %>% filter(ForestAllianceType=="Redwood") %>% 
  group_by(Species) %>% 
  summarise(IV.mean = mean(IV), IV.sd = sd(IV), IV.se=sd(IV)/sqrt(n()),
            n.stems = sum(n.stems), sum.BA = sum(sum.BA),
            n.plots = n()) %>% 
  mutate(IV.tot = (n.stems/sum(n.stems) + sum.BA/sum(sum.BA))/2 )

pairs(species.sum.ME[,-1])
pairs(species.sum.Red[,-1])

reorder <- function(df, var){
  df=as.data.frame(df) 
  ord <- df$Species[order(df[,var], decreasing = T)]
  factor(df$Species, levels = ord)
  }

#I think I will use n.plots to inform which plants will be inoculated
dat=species.sum.Red
dat=species.sum.ME
dat$Species2 <- reorder(dat, "IV.tot")
ggplot(dat, aes(Species2, IV.tot)) +
  geom_col()

dat$Species2 <- reorder(dat, "n.plots")
ggplot(dat %>% filter(n.plots>0), aes(Species2, n.plots)) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
 # geom_point(aes(Species2, IV.tot*200)) + geom_point(aes(Species2, IV.mean*200), color="blue")
ggplot(dat, aes(log(n.plots), log(IV.tot))) +
  geom_point() +
  geom_text(aes(label = Species), hjust=1, vjust=1)
ggplot(dat, aes(log(IV.tot), log(IV.mean)) ) +
  geom_point() +
  geom_errorbar( aes(ymin=log(IV.mean-2*IV.se), ymax=log(IV.mean+2*IV.se) ))+
  geom_text(aes(label = Species), hjust=1, vjust=1)
ggplot(dat, aes((IV.tot), (IV.mean)) ) +
  geom_point() +
  geom_errorbar( aes(ymin=(IV.mean-2*IV.se), ymax=(IV.mean+2*IV.se) ))+
  geom_text(aes(label = Species), hjust=1, vjust=1)

#write.csv(species.sum.Red, "output/redwoodspecies.csv", row.names = F)
#write.csv(species.sum.ME, "output/mixedeverspecies.csv", row.names = F)
##################################################################
##################################################################
#I'm using Kerri's dataset from now on. there were slight differences from allison's dataframe and I'm not sure if it has to do with my post-processing that I did as a first year. 
df <- read.csv("~/Desktop/big sur/Big Sur Github/outputs/allison_firstyr_w.species.csv", header=T)
df2 <- df %>% 
  select("BSPlotNumber", "ForestAllianceType", "StemNumber", "TreeNumber", "Species", "DBH", "Status", "CankerTrunkSymptoms":"TwigBasalSproutSymptoms","P_ramorum_plotlevel") %>% 
  filter(!is.na(Status)) %>% 
  filter(!Species=="MISSING") %>% 
  #get disease symptoms out. Only trusting Y. if it has symptoms and p.ram is plot level, counts as diseased (that's what Haas 2011 did)
  mutate_at(.vars = vars(c("CankerTrunkSymptoms":"TwigBasalSproutSymptoms")), .funs = function(x)(Rename(c(0,0,0,1), c(NA, "N", "P", "Y"), x))) %>% 
  mutate(nsymptoms = rowSums(.[8:12]),
         infected = ifelse(P_ramorum_plotlevel ==1 & nsymptoms>0,1,0 )) %>% 
  #basal area (m2), relative 
  mutate(BA = ifelse(is.na(DBH), 0, pi*(DBH/200)^2)) %>% 
  select(-c(CankerTrunkSymptoms:TwigBasalSproutSymptoms))
df2 %>% group_by(ForestAllianceType, P_ramorum_plotlevel) %>% summarise(n_distinct(BSPlotNumber)) 

