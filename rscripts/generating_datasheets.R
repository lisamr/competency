library(dplyr)
rm(list=ls()) #clear environment

#writing a script to create field data sheets for plant collection

#trip 1: collect leaves from 32 individuals for each species. broad-leaved plants
broads <- c("UMCA", "QUAG", "LIDE", "ARME", "QUPA", "QUCH", "TODI", "HEAR", "CEOL", "ACMA")
broads <- sort(broads) #put in alphabetical order
ninds <- 32

#make a table with species, individual, easting, northing, general location, notes. 
df <- expand.grid(species=broads, ind=1:ninds, easting=NA, northing=NA, general_location=NA, notes=NA)
df <- df %>% arrange(species)
head(df)

write.csv(df, file="data_sheets2019/collection_broadlvs.csv", row.names = F, na = "") #specify NA as blank

#trip 2: collect leaves from 32 individuals of each conifer species and 10 individuals of UMCA and LIDE
conifers <- c("PIPO", "SESE", "PSME")

#make the table. mind the difference in # individuals
df2 <- expand.grid(species=conifers, ind=1:ninds, easting=NA, northing=NA, general_location=NA, notes=NA)
df3 <- expand.grid(species=c("UMCA", "LIDE"), ind=1:10, easting=NA,northing=NA, general_location=NA, notes=NA)
df4 <- rbind(df2, df3) %>% arrange(species) 
head(df4)

write.csv(df4, file="data_sheets2019/collection_conifers.csv", row.names = F, na = "") #specify NA as blank 

############################################################
############################################################
#master data sheet for processing
df5 <- expand.grid(leafID=NA, leafID2=NA, species=broads, ind=1:ninds, trt=c("T", "C"), spore_assay=c("S", "C"))
#create empty columns
colvector <- c("processor", "symptoms", "sp_viable", "inf_viable", "counter", "date_counted", paste0("count", 1:4), paste0("vol", 1:4), "lesion_size", "notes") #cols to be filled with NA
df5[colvector] <- NA #create those col and fill with NA
head(df5)

#Oh, actually, I'm going to remove the samples where chlamydospore counting isnt possible. That is going to be quch, qupa and quag. 
df5 <- df5 %>% filter(!(species %in% c("QUCH", "QUPA", "QUAG") & spore_assay == "C"))

#assign unique IDs
df5 <- df5 %>% 
  mutate(leafID = group_indices(., species, trt, ind), 
         leafID2 = paste0(leafID, spore_assay)) %>% 
  arrange(leafID)
head(df5)
length(df5$leafID2)#cut it down from 1280 to 1088. that's nice! unfortunately it only cuts out chlamydos and it doesn't reduce the number of petri plates. 

View(df5)
#master data sheet for processing detached leaf assays
df6a <- expand.grid(leafID=NA, leafID2=NA, species=conifers, ind=1:ninds, trt=c("T", "C"), spore_assay=c("S", "C"))
df6b <- expand.grid(leafID=NA, leafID2=NA, species=c("UMCA", "LIDE"), ind=1:10, trt=c("T", "C"), spore_assay=c("S", "C"))
df6 <- rbind(df6a, df6b)
df6[colvector] <- NA
head(df6)

#assign unique IDs. making ID continuous from the leaf discs. 
df6 <- df6 %>% 
  mutate(leafID = group_indices(., species, trt, ind) + max(df5$leafID), 
         leafID2 = paste0(leafID, spore_assay)) %>% 
  arrange(leafID)
head(df6)



#########
#export master files to csv
write.csv(df5, file="data_sheets2019/master_broadlvs.csv", row.names = F, na = "") #specify NA as blank 
write.csv(df6, file="data_sheets2019/master_conifers.csv", row.names = F, na = "") #specify NA as blank 

