#turn the contrast table into a matrix
library(tidyverse)
df <- read.csv('output/sporangia_both/contrasts.csv')

head(df)

broad <- df %>% filter(assay == 'broad') %>% droplevels()
broadmat <- xtabs(broad[,5] ~ broad[,6] + broad[,7]) %>% as.matrix()
broadmatsig <- xtabs(broad[,8] ~ broad[,6] + broad[,7])

conif <- df %>% filter(assay == 'conif') %>% droplevels()
conifmat <- xtabs(conif[,5] ~ conif[,6] + conif[,7])
conifmatsig <- xtabs(conif[,8] ~ conif[,6] + conif[,7])

write.table(broadmat, 'output/sporangia_both/contrasts_matrix_broad.csv', sep = ',')
write.table(conifmat, 'output/sporangia_both/contrasts_matrix_conif.csv', sep = ',')
write.table(broadmatsig, 'output/sporangia_both/contrasts_matrix_broadsig.csv', sep = ',')
write.table(conifmatsig, 'output/sporangia_both/contrasts_matrix_conifsig.csv', sep = ',')