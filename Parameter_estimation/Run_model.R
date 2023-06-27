rm(list=ls())
####### load libraries
library(deSolve)
library(openxlsx)
library(DriftAndSelection)
###### load functions

patient.id <- 'N1' # specify patient ID here
sort <- "CD34" # specify cell sort here

sample.info <- read.xlsx("~/Metadata/Sample_info.xlsx")
rownames(sample.info) <- sample.info$Sample
age <- sample.info[patient.id,]$Age*365

depth=sample.info[patient.id,paste("Depth", sort, sep=".")]

load("./pyABC/Blood/SNVs_indels.RData")
snvs <- match.arg(sort, c("snvs.cd34", "snvs.mnc", "snvs.mnc_minus_t", "snvs.pb"))
snvs <- list(snvs.cd34[[patient.id]])

use.sensitivity = F

source("~/Parameter_estimation/Bayesian_fit.R")
