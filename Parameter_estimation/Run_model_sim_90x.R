rm(list=ls())
####### load libraries
library(deSolve)
library(openxlsx)
library(FLORENCE)
###### load functions

patient.id='STS1' 

sample.info <- read.xlsx("Metadata/Sample_information_simulated_data.xlsx", sheet = 1)
rownames(sample.info) <- sample.info$SampleID
age <- sample.info[patient.id,]$Age*365

depth=90

load("RData/Simulated_data/SNVs_90x.RData")
snvs <- list(snvs[[patient.id]])

use.sensitivity = F
source("Parameter_estimation/Bayesian_fit.R")
