library(dplyr)

#writing a script to create field data sheets for plant collection

#trip 1: collect leaves from 32 individuals for each species. broad-leaved plants
broads <- c("UMCA", "QUAG", "LIDE", "ARME", "QUPA", "QUCH", "TODI", "HEAR", "CEOL", "ACMA")
ninds <- 32

#make a table with species, individual, easting, northing, general location, notes. 
df <- expand.grid(Species=broads, Ind=1:32, Easting=NA, Northing=NA, General_location=NA, Notes=NA)
df <- df %>% arrange(Species)
head(df)

write.csv(df, file="data_sheets/collection_broadlvs.csv", row.names = F, na = "") #specify NA as blank

#trip 2: collect leaves from 32 individuals of each conifer species and 10 individuals of UMCA and LIDE
conifers <- c("PIPO", "SESE", "PSME")

#make the table. mind the difference in # individuals
df2 <- expand.grid(Species=conifers, Ind=1:32, Easting=NA, Northing=NA, General_location=NA, Notes=NA)
df3 <- expand.grid(Species=c("UMCA", "LIDE"), Ind=1:10, Easting=NA, Northing=NA, General_location=NA, Notes=NA)
df4 <- rbind(df2, df3) %>% arrange(Species) 
head(df4)

write.csv(df4, file="data_sheets/collection_conifers.csv", row.names = F, na = "") #specify NA as blank 
