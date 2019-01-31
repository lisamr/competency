setwd("~/Desktop/Competency project")

#I inoculated UMCA leaf discs with 15 ul of 4000-ish sporangia/ml and processed them 5 days later. All 5 leaves were from the same tree. I scraped them with different intensities, which is why I think there's so much variation. The colors are the plates and the samples on those plates should have been processed together. The blue plate was the last one I processed and I think I got the hang of it--don't scrape so hard. 
dat <- read.csv("data/umca_trial.csv", header = T)
(dat)
plot(data=dat, counts_mean~leaf, col=plate)
plot(data=dat, lesion~leaf, col=plate)
plot(data=dat, lesion~counts_mean, col=plate)

#here's the data sheet for my samples. There's a lot. I need to figure out how many tubes I need to label for the sporangia and chlamydospore counts. this will inform how much KOH and LPCB I need to make.
df <- read.csv("data/competencydata.csv", header = T)
head(df)

library(dplyr)
#make summary table to see the overall design
df_sum <- df %>% 
  group_by(species, spore.type, assay.type, trt) %>% 
  summarise("n.ind"=length(unique(ind)), "leaves_per_ind"=length(unique(leaf)), "total leaves/discs"=n.ind * leaves_per_ind)
#there are 16 total species and 1152 leaf discs/leaves to be inoculated and counted. :/
list("n.spp"=length(unique(df_sum$species)), "n.leaves"=sum(df_sum$`total leaves/discs`))

#a table with each species, individual and the number of leaves collected. Probably need to update this tomorrow when I actually cut the leaf discs because it will likely change. I should print that and update it while cutting discs.
n.leaves <- df %>% 
  group_by(species, ind) %>% 
  summarise(length(unique(leaf)))

write.csv(n.leaves, "data/nleavescollected.csv", row.names = F)

#Here I need to know how many of the tubes need to be filled and labeled. 
df %>% 
  group_by(spore.type) %>% 
  count(tubes)

#most of the sporangia will be collected in .2 ml tubes. Label them 1-552. I think it's fine to label just the first and last of each strip. Will save lots of time.
a <- df %>% 
  filter(spore.type=="sporangia",tubes==".2 ml")
a
a$No

#UMCA and LIDE sporangia are those collected in 1.5 ml tubes. Label them 553-576
df %>% 
  filter(spore.type=="sporangia",tubes=="1.5 ml")

#all of the leaves will be added to 2 ml tubes filled with KOH. label them 577-1152. There will be 576 tubes.
b <- df %>% 
  filter(spore.type=="chlamydo")
b
b$No
nrow(b)

######
#making spreadsheet with species, tube numbers for treatments and type
head(df)
df$species.ind <- interaction(df$species, df$ind)
df2 <- df %>% 
  filter(spore_type=="S", assay_code == "L") %>% 
  group_by(species.ind, trt) %>% 
  summarize("rep1"= min(No), "rep2"= min(No)+1,"rep3"= min(No)+2, "rep4"= min(No)+3, "rep5"= min(No)+4, "rep6"= max(No)) %>% 
  arrange(rep1)

df3 <- df %>% 
  filter(spore_type=="C", assay_code == "L") %>% 
  group_by(species.ind, trt) %>% 
  summarize("rep1"= min(No), "rep2"= min(No)+1,"rep3"= min(No)+2, "rep4"= min(No)+3, "rep5"= min(No)+4, "rep6"= max(No)) %>% 
  arrange(rep1)

df4 <- df %>% 
  filter(spore_type=="S", assay_code == "D") %>% 
  group_by(species.ind, trt) %>% 
  summarize("rep1"= min(No), "rep2"= min(No)+1,"rep3"= min(No)+2, "rep4"= min(No)+3, "rep5"= min(No)+4, "rep6"= max(No)) %>% 
  arrange(rep1)

df5 <- df %>% 
  filter(spore_type=="C", assay_code == "D") %>% 
  group_by(species.ind, trt) %>% 
  summarize("rep1"= min(No), "rep2"= min(No)+1,"rep3"= min(No)+2, "rep4"= min(No)+3, "rep5"= min(No)+4, "rep6"= max(No)) %>% 
  arrange(rep1)

write.csv(df2, "output/sporangiatube.leafdiscs_spreadsheet.csv", row.names = F)
write.csv(df3, "output/chlamydotube.leafdiscs_spreadsheet.csv", row.names = F)
write.csv(df4, "output/sporangiatube.detached_spreadsheet.csv", row.names = F)
write.csv(df5, "output/chlamydotube.detached_spreadsheet.csv", row.names = F)
