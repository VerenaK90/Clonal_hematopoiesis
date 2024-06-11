#### Analyze SNVs in human brain regions
source("./Settings.R")

### read in the SNV data 

SNV.data.lieber <- read.xlsx("./Published_data/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s2.xlsx", sheet = 4)
SNV.data.yale <- read.xlsx("Published_data/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s2.xlsx", sheet = 2)
SNV.data.harvard <- read.xlsx("Published_data/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s2.xlsx", sheet = 3)

### meta.data

donor.info.lieber <- read.xlsx("./Published_data/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s1.xlsx", sheet = 6, startRow = 7)
donor.info.lieber <- donor.info.lieber[!is.na(donor.info.lieber$Case.ID),]
donor.info.yale <- read.xlsx("./Published_data/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s1.xlsx", sheet = 2, startRow = 7)
donor.info.yale <- donor.info.yale[!is.na(donor.info.yale$Case.ID),]
donor.info.harvard <- read.xlsx("./Published_data/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s1.xlsx", sheet = 4, startRow = 7)
donor.info.harvard <- donor.info.harvard[!is.na(donor.info.harvard$Case.ID),]
sample.info.lieber <- read.xlsx("./Published_data/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s1.xlsx", sheet = 7)
sample.info.lieber$Age <- sapply(sample.info.lieber$Case.ID, function(x){
  donor.info.lieber[donor.info.lieber$Case.ID == x,]$Age
})
sample.info.lieber$Sex <- sapply(sample.info.lieber$Case.ID, function(x){
  donor.info.lieber[donor.info.lieber$Case.ID == x,]$Sex
})
sample.info.lieber$Phenotype <- sapply(sample.info.lieber$Case.ID, function(x){
  donor.info.lieber[donor.info.lieber$Case.ID == x,]$Phenotype
})
sample.info.lieber$Center <- "LIBD"
sample.info.lieber$Location <- sample.info.lieber$Origin.of.Tissue
sample.info.lieber$CellType <- sample.info.lieber$Prep.Method

sample.info.yale <- read.xlsx("./Nextcloud/Hematopoiesis/Blood_oxford/Mutation_analysis/Brain/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s1.xlsx", sheet = 3)
sample.info.yale$Age <- sapply(sample.info.yale$Case.ID, function(x){
  donor.info.yale[donor.info.yale$Case.ID == x,]$`Age.(years)`
})
sample.info.yale$Sex <- sapply(sample.info.yale$Case.ID, function(x){
  donor.info.yale[donor.info.yale$Case.ID == x,]$Sex
})
sample.info.yale$Phenotype <- sapply(sample.info.yale$Case.ID, function(x){
  donor.info.yale[donor.info.yale$Case.ID == x,]$Phenotype
})
sample.info.yale$Center <- "Yale"
sample.info.yale$Location <- sapply(sample.info.yale$Sample.ID, function(x){
  strsplit(x, split="-")[[1]][2]
})
sample.info.yale$CellType <- sapply(sample.info.yale$Sample.ID, function(x){
  strsplit(x, split="-")[[1]][3]
})
sample.info.harvard <- read.xlsx("./Nextcloud/Hematopoiesis/Blood_oxford/Mutation_analysis/Brain/Bae_et_al/science.abm6222_tables_s1_s2_s6_and_s7/science.abm6222_table_s1.xlsx", sheet = 5)
sample.info.harvard$Age <- sapply(sample.info.harvard$Case.ID, function(x){
  donor.info.harvard[donor.info.harvard$Case.ID == x,]$`Age.(years)`
})
sample.info.harvard$Sex <- sapply(sample.info.harvard$Case.ID, function(x){
  donor.info.harvard[donor.info.harvard$Case.ID == x,]$Sex
})
sample.info.harvard$Phenotype <- sapply(sample.info.harvard$Case.ID, function(x){
  donor.info.harvard[donor.info.harvard$Case.ID == x,]$Phenotype
})
sample.info.harvard$Center <- "Harvard"
sample.info.harvard$Location <- sample.info.harvard$Origin.of.Tissue
sample.info.harvard$CellType <- "bulk"

sample.info.bae.et.al <- rbind(sample.info.lieber[,c("Case.ID", "Sample.ID", "Age", "Sex", "Phenotype", "Average.coverage", "Center", "Location", "CellType")],
                           sample.info.yale[,c("Case.ID",  "Sample.ID","Age", "Sex", "Phenotype", "Average.coverage", "Center", "Location", "CellType")],
                           sample.info.harvard[,c("Case.ID",  "Sample.ID","Age", "Sex", "Phenotype", "Average.coverage", "Center", "Location", "CellType")])

# harmonize location annotations
sample.info.bae.et.al$Location[sample.info.bae.et.al$Location == "Frontal Cortex - generic"] <- "Generic FC"
sample.info.bae.et.al$Region <- sample.info.bae.et.al$Location
sample.info.bae.et.al$Region[sample.info.bae.et.al$Location %in% c("CX", "DLPFC", "Generic PFC", "PFC",
                                                    "Generic FC", "BA10", "BA46", "BA9/46",
                                                    "BA10 (Right)", "BA9/10", "BA9")] <- "Cortex"
sample.info.bae.et.al$Region[sample.info.bae.et.al$Location=="STR"] <- "Striatum"

### collect SNVs

snvs <- list()

for(i in unique(SNV.data.lieber$Individual)){
  
  for(tissue in c("DLPFC", "Hippocampus")){
    
    to.store <- SNV.data.lieber[SNV.data.lieber$Individual==i,]
    to.store$Depth <- sapply(to.store[,paste0(tissue, ".REF;ALT")], function(x){
      sum(as.numeric(strsplit(x, split=";")[[1]]))
    } )
    to.store$VAF <- sapply(to.store[,paste0(tissue, ".REF;ALT")], function(x){
      x <- as.numeric(strsplit(x, split=";")[[1]])
      x[2]/sum(x)
    } )
    to.store$varCounts <- to.store$VAF * to.store$Depth
    to.store <- to.store[to.store$VAF>0,]

    snvs[[paste(i, tissue, sep="_")]] <- to.store[,-which(colnames(to.store) %in% c("DLPFC.VAF", "DLPFC.REF;ALT", 
                                                            "Hippocampus.VAF", "Hippocampus.REF;ALT"))]
    
  
  }
  
}

for(i in unique(SNV.data.harvard$Individual)){

  to.store <- SNV.data.harvard[SNV.data.harvard$Individual==i,]
  to.store$Depth <- sapply(to.store[, "Bulk.REF;ALT"], function(x){
    sum(as.numeric(strsplit(x, split=";")[[1]]))
  } )
  to.store$VAF <- sapply(to.store[, "Bulk.REF;ALT"], function(x){
    x <- as.numeric(strsplit(x, split=";")[[1]])
    x[2]/sum(x)
  } )
  to.store$varCounts <- to.store$VAF * to.store$Depth
  to.store <- to.store[to.store$VAF>0,]
  
  snvs[[i]] <- to.store[,-which(colnames(to.store) %in% c("Bulk.VAF", "Bulk.REF;ALT"))]
  
 }

for(i in unique(SNV.data.yale$Individual)){
  
  for(tissue in c("STR-Bulk", "STR-INT", "STR-OLI", "STR-MSN", "STR-ASTMIG",
                  "CX-Bulk", "CX-INT", "CX-OLI", "CX-PYR", "CX-ASTMIG")){
    
    to.store <- SNV.data.yale[SNV.data.yale$Individual==i,]
    to.store$Depth <- sapply(to.store[,paste0(tissue, ".REF;ALT")], function(x){
      sum(as.numeric(strsplit(x, split=";")[[1]]))
    } )
    to.store$VAF <- sapply(to.store[,paste0(tissue, ".REF;ALT")], function(x){
      x <- as.numeric(strsplit(x, split=";")[[1]])
      x[2]/sum(x)
    } )
    to.store$varCounts <- to.store$VAF * to.store$Depth
    to.store <- to.store[!is.na(to.store$VAF) & to.store$VAF>0,]
    
    snvs[[paste(i, tissue, sep="-")]] <- to.store[,-which(colnames(to.store) %in% c("STR-Bulk.VAF", "STR-INT.VAF", "STR-OLI.VAF", "STR-MSN.VAF", 
                                                                                          "STR-ASTMIG.VAF", "CX-Bulk.VAF", "CX-INT.VAF", "CX-OLI.VAF", 
                                                                                    "CX-PYR.VAF", "CX-ASTMIG.VAF",
                                                                                    "STR-Bulk.REF;ALT", "STR-INT.REF;ALT", "STR-OLI.REF;ALT", 
                                                                                    "STR-MSN.REF;ALT", "STR-ASTMIG.REF;ALT", "CX-Bulk.REF;ALT", 
                                                                                    "CX-INT.REF;ALT", "CX-OLI.REF;ALT", "CX-PYR.REF;ALT", "CX-ASTMIG.REF;ALT"))]
    
    }
  
}

save(snvs, file="./RData/Bae_et_al/SNVs_brain.RData")
