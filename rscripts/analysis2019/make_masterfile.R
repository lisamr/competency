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

#remove samples that were manually done. includes mostly 0 and 100 or near 100. also ones without both total area and perc_lesion
df2 <- df %>% filter(
  !(is.na(perc_lesion)|is.na(total_area)),
  perc_lesion<99)
head(df2)

#plot the relationship between area and lesion. pretty damn linear.
ggplot(df2, aes(total_area, perc_lesion))+
  geom_point()+
  geom_smooth(method = 'lm')
m1 <- lm(data=df2, perc_lesion~0+total_area) #0 constrains model to pass through 0 at intercept
summary(m1)

#create data to predict lesion size
ndat <- df %>% filter(!is.na(total_area), is.na(perc_lesion))
ndat$perc_lesion <- predict(m1, newdata = ndat)
head(ndat)

#plot points with new data
ggplot(df2, aes(total_area, perc_lesion))+
  geom_point()+
  geom_smooth(method = 'lm') +
  geom_point(data=ndat, aes(total_area, perc_lesion), color='red')

#merge new data with original data
dfinal <- rbind(ndat, df) %>%  
  filter(!duplicated(leafID2)) %>% 
  arrange(leafID)
#View(dfinal)

#merge gps data to master file
dfinal <- merge(gps, dfinal, by=c("species", "ind")) %>% arrange(leafID)

write.csv(dfinal, "master_tall.csv", row.names = F)
