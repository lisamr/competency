#finalize master file.
#need to convert lesion samples that were manually done from pixels to %lesion
#create linear model between area and lesion to figure out %lesion.

setwd("~/Box/Competency project/competency.git/data2019")
library(dplyr)
library(ggplot2)
rm(list=ls())#clear environment

df <- readxl::read_excel('master_bothassays.xlsx')
gps <- readxl::read_excel('plant_collection2019.xlsx')

#convert random numbers to actual dates
df$date_counted <- as.Date(as.numeric(df$date_counted), origin = as.Date("1899-12-31"))
df$perc_lesion <- as.numeric(df$perc_lesion) #character for some reason
head(df)

#remove samples that were manually done and the conifer assay samples. includes mostly 0 and 100 or near 100. also ones without both total area and perc_lesion
df2 <- df %>% filter(
  !(is.na(perc_lesion)|is.na(npixels)),
  perc_lesion<99, leafID<530)
head(df2)
tail(df2)
#plot the relationship between area and lesion. pretty damn linear.
ggplot(df2, aes(npixels, perc_lesion))+
  geom_point()+
  geom_smooth(method = 'lm')
m1 <- lm(data=df2, perc_lesion~0+npixels) #0 constrains model to pass through 0 at intercept
summary(m1)

#create data to predict lesion size
ndat <- df %>% filter(!is.na(npixels), is.na(perc_lesion))
ndat$perc_lesion <- predict(m1, newdata = ndat)
head(ndat)

#plot points with new data
ggplot(df2, aes(npixels, perc_lesion))+
  geom_point()+
  geom_smooth(method = 'lm') +
  geom_point(data=ndat, aes(npixels, perc_lesion), color='red')

#merge new data with original data
ndat2 <- ndat %>% select(leafID2, perc_lesionp=perc_lesion)
as.numeric(df$leaf_area_cm2)
dfinal <- df %>% 
  left_join(ndat2, by = "leafID2") %>% 
  mutate(
    perc_lesion = case_when(
    is.na(perc_lesion) ~ perc_lesionp,
    T ~ perc_lesion)) %>% 
  mutate(
    perc_lesionp = NULL,
    npixels = NULL,
    leaf_area_cm2 = as.numeric(leaf_area_cm2),
    lesion_area_cm2 = leaf_area_cm2*(perc_lesion/100)) %>% 
  mutate_at(.vars = c("perc_lesion", "leaf_area_cm2", "lesion_area_cm2"), ~round(., 2)) %>% 
  arrange(leafID)

#merge gps data to master file
dfinal <- merge(gps, dfinal, by=c("species", "ind")) %>% arrange(leafID)

write.csv(dfinal, "master_tall.csv", row.names = F)
