setwd("~/Desktop/Competency project/data")

############################
#Inoculum viability
############################

#read in data and view it
df <- read.csv("inoculum_viability.csv", header = T)
head(df)
df <- df[,-c(7,8)] #there where extra hidden columns in the spreadsheet that I'm removing

#making boxplots to compare the groups, which are when-who-assays
#first load the packages
library(ggplot2)
library(dplyr)

#make combination groups to compare with boxplots
df$whenwho <- interaction(df$When, df$Who)

#plot it
ggplot(aes(y = CFUs, x = whenwho), data = df) + 
  geom_boxplot()+
  facet_wrap(~assay)

#plot just before and after
ggplot(aes(y = CFUs, x = When), data = df) + 
  geom_boxplot()+
  facet_wrap(~assay)

#there's a decent difference between assays. I think that has to do with the amount of time before plating and agitation that the leaf disc inoculum had. I plated all of these at the same time, despite the fact that the leaf disc inoculum was made about 2ish hours before the detached leaves. I looked at the inoculum under the scope and noticed more zoospores had released in the leaf discs, so I think that would explain the difference bewtween assays. It doesn't look like there's much of a difference between Lisa and Sebastian in the detached leaves, but kinda in the leaf discs. I really can't tell for sure though with only 3 counts of CFUs. Let's call it fine.

#as for the second plot (before and after), it looks good for detached leaves. There is a difference between before and after for leaf disc, likely due to the zoospores encysting. Hopefully it didn't play a major part, but at least we randomized the plates. The "after" solution was made first though, so it had the longest time for the zoospores to encyst. 

#counts for the leaf disc inoculum. These counts were from 5 ul drops. I made it and then split it into 2 aliquots, one for me and the other for sebastian. I assumed they were equal and did not make sure they were the same. I hope it was okay, but looking at the leaf disc box plots makes me worry a little bit. 
v <- c(26, 19, 22, 21, 21, 29)

#standardizing to sporangia/mL
v2 <- v/5*1000
se <- function(x){
  sd(x)/sqrt(length(x))
}

#range for leaf disc inoculum: 3990-5210, mu=4600
CI <- 2*se(v2)
mu <- mean(v2)
mu + CI
mu - CI

#counts for the detached leaf dip inoculations. These counts were from 5 ul drops and I counted them seperately for me and sebastian.
L <- c(31, 28, 21, 29, 25)
S <- c(30, 31, 28, 22, 28)

#standardizing to sporangia/mL
L2 <- L/5*1000 
S2 <- S/5*1000

#for Lisa detached leaves: 4660-6010, mu=5360
CI <- 2*se(L2)
mu <- mean(L2)
mu + CI
mu - CI

#for Lisa detached leaves: 4930-6180, mu=5560
CI <- 2*se(S2)
mu <- mean(S2)
mu + CI
mu - CI

