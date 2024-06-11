rm(list=ls())
####### load libraries
library(deSolve)
library(openxlsx)
library(SCIFER)
###### load functions

sample <- "AN02255" # specify sample ID

sample.info <- read.xlsx("MetaData/Supplementary Tables.xlsx", startRow = 7, sheet = 8)
sample.info$Age <- as.numeric(sample.info$Age)
sample.info$Average.coverage <- as.numeric(sample.info$Average.coverage)
sample.info$Sample.ID <- replace(sample.info$Sample.ID, sample.info$Sample.ID == "M3663M ", "M3663M")

donor <- sample.info[sample.info$Sample.ID==sample,]$Case.ID
age <- sample.info[sample.info$Case.ID==donor & sample.info$Sample.ID==sample,]$Age*365
depth <- round(sample.info[sample.info$Case.ID==donor & sample.info$Sample.ID==sample,"Average.coverage"])
if(depth < 150){
min.vaf <- 0.05
min.clone.size = 0.05
min.prior.size=0.01

}else{

min.vaf <- 0.02
min.clone.size = 0.01
min.prior.size=0.001
}

load("./RData/Bae_et_al/SNVs_brain.RData")

if(sample %in% names(snvs)){
	snvs <- list(snvs[[sample]])
}else{
	snvs <- list(snvs[[donor]])
}

use.sensitivity = F
source("~/Parameter_estimation/Bayesian_fit.R")
