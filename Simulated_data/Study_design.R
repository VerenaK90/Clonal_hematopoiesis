source("./Metadata/Settings.R")
######################### ######################### ######################### ######################### ######################### 
## Design the simulation study for the ROC curve. To this end, use the VAFs from the simulated trees. Overall, this will generate a snvs-object, storing all snvs and a sample information sheet

study.directory <- "Simulated_data/"
sample.info <- read.xlsx("Metadata/Supplementary Table.xlsx", sheet = 1, startRow = 3)
rownames(sample.info) <- sample.info$SampleID

## simulate the VAFs for each sample at 90x coverage

depth <- 90
snvs <- list()
for(i in sample.info$SampleID){
  set.seed(which(sample.info$SampleID==i)*368)
  # load the simulated data
  load(paste0("./RData/Simulated_data/",
              sample.info[i,]$Tree_simulation_file))
  # access the VAFs
  true.vafs <- vafs[[as.character(sample.info[i,]$Age)]]
  # transform them into simulated VAFs
  sim.vafs <- as.data.frame(t(sapply(true.vafs, function(x){
    cov <- rpois(1, depth)
    vaf <- rbinom(1, depth, x)/cov
    c(VAF=vaf, Depth=cov, varCounts = vaf*cov)
  })))
  snvs[[i]] <- sim.vafs[sim.vafs$VAF>0,]
}

save(snvs, file="RData/Simulated_data/SNVs_90x.RData")

# 270x coverage
depth <- 270
snvs <- list()

for(i in sample.info$SampleID){
  set.seed(which(sample.info$SampleID==i)*368*2)

  load(paste0("./RData/Simulated_data/",
              sample.info[i,]$Tree_simulation_file))
  true.vafs <- vafs[[as.character(sample.info[i,]$Age)]]
  sim.vafs <- as.data.frame(t(sapply(true.vafs, function(x){
    cov <- rpois(1, depth)
    vaf <- rbinom(1, depth, x)/cov
    c(VAF=vaf, Depth=cov, varCounts = vaf*cov)
  })))
  snvs[[i]] <- sim.vafs[sim.vafs$VAF>0,]
}

save(snvs, file="RData/Simulated_data/SNVs_270x.RData")


# 30x coverage
depth <- 30
snvs <- list()

for(i in sample.info$SampleID){
  set.seed(which(sample.info$SampleID==i)*368*2)
  load(paste0("./RData/Simulated_data/",
              sample.info[i,]$Tree_simulation_file))
  true.vafs <- vafs[[as.character(sample.info[i,]$Age)]]
  sim.vafs <- as.data.frame(t(sapply(true.vafs, function(x){
    cov <- rpois(1, depth)
    vaf <- rbinom(1, depth, x)/cov
    c(VAF=vaf, Depth=cov, varCounts = vaf*cov)
  })))
  snvs[[i]] <- sim.vafs[sim.vafs$VAF>0,]
}

save(snvs, file="RData/Simulated_data/SNVs_30x.RData")
