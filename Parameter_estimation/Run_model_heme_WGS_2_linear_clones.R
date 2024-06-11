rm(list=ls())
####### load libraries
library(SCIFER)
library(data.table)
library(openxlsx)
###### load functions

patient.id <- 'A1' # specify patient ID here

sample.info <- read.xlsx("MetaData/Supplementary Tables.xlsx", sheet = 2, startRow = 6)
rownames(sample.info) <- sample.info$Paper_ID

load("./RData/WGS_heme/SNVs.RData")

depth=270
min.vaf <- 0.02
min.clone.size = 0.01
min.prior.size=0.001

use.sensitivity = F

snvs <- list(snvs.cd34.deep[[patient.id]])

mother.daughter <- matrix(c(1, 2, 2,3), byrow=T, ncol=2) # clonal topology according to linear evolution

age <- sample.info[patient.id,]$Age*365

source("~/Parameter_estimation/Bayesian_fit_multiclone.R")
