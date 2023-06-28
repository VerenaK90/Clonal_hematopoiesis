rm(list=ls())
####### load libraries
library(deSolve)
library(openxlsx)
library(DriftAndSelection)
###### load functions

patient.id='STS1' 

sample.info <- read.xlsx("Metadata/Sample_information_simulated_data.xlsx", sheet = 1)
rownames(sample.info) <- sample.info$SampleID
age <- sample.info[patient.id,]$Age*365

depth=270

min.vaf <- 0.01
min.clone.size = 0.01
min.prior.size=0.001

load("RData/Simulated_data/SNVs_270x.RData")
snvs <- list(snvs[[patient.id]])

use.sensitivity = F
source("Parameter_estimation/Bayesian_fit.R")
