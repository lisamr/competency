setwd("~/Desktop/Competency project")
df <- read.csv("data/spore data recording sheet JUNE 2018.csv", header = T)
head(df)

library(dplyr)
library(ggplot2)

df2 <- df %>% 
  select(No, code, ind, species, trt, rep, sporangia.counts, X.2, X.3, X.4) %>% 
  filter(!is.na(sporangia.counts), trt=="T")

head(df2)
df2$meansp <- apply(df2[,7:10], 1, function(x) mean(x, na.rm=T))

df2$speciesind <- interaction(df2$species, df2$ind)

#plot it
df2$speciesind <- factor(df2$speciesind, levels = c("ACMA.1", "ACMA.2", "ACMA.3", "UMCA.1", "UMCA.2", "UMCA.3", "QUKE.1", "QUKE.2", "QUKE.3"))

ggplot(aes(y = meansp, x = speciesind), data = df2) + 
  geom_boxplot()
ggplot(aes(y = meansp, x = species), data = df2) + 
  geom_boxplot()

#remove outlier ind
df3 <- df2 %>% 
  filter(!speciesind=="ACMA.3", !speciesind=="QUKE.3")
ggplot(aes(y = meansp, x = speciesind), data = df3) + 
  geom_boxplot()
ggplot(aes(y = meansp, x = species), data = df3) + 
  geom_boxplot()

#getting mean and se of the data. kinda crappy
tapply(df2$meansp, df2$speciesind, function(x) list(mean=mean(x), se=2*sd(x)/length(x)))
x <- mean(tapply(df2$meansp, df2$speciesind, function(x) relse=(sd(x)/sqrt(length(x)))/mean(x)))


#simulating stuff
####################################
####################################
rnorm(10, mean(df2$meansp), sd(df2$meansp))

f <- function(n){
  test <- tapply(df2$meansp, df2$speciesind, function(x) rnorm(n, mean(df2$meansp), sd(df2$meansp)))
  test2 <- sapply(test, function(x) list(relse=(sd(x)/sqrt(length(x)))/mean(x)))
  test3 <- mean(sapply(test2, mean))
  points(n, test3)
  
}

f(6)
test <- tapply(df2$meansp, df2$speciesind, function(x) rnorm(6, mean(df2$meansp), sd(df2$meansp)))
test2 <- sapply(test, function(x) list(relse=(sd(x)/sqrt(length(x)))/mean(x)))
test3 <- mean(sapply(test2, mean))
plot(6, test3, xlim=c(0, 101), ylim=c(0, 2))
points(6,x, col="red")

