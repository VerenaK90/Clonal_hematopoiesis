rm(list=ls())
####### load libraries
library(deSolve)
library(openxlsx)
library(SCIFER)
###### load functions

patient.id <- 'N1' # specify patient ID here
sort <- "CD34" # specify cell sort here, should be either CD34, MNC, MNC_minus_T, PB_gran

sample.info <- read.xlsx("MetaData/Supplementary Tables.xlsx", sheet = 2, startRow = 6)
rownames(sample.info) <- sample.info$Paper_ID

load("./RData/WGS_heme/SNVs.RData")

if(sort=="CD34"){
  depth <- sample.info[patient.id,"Coverage.WGS.CD34+.1"]
  snvs <- snvs.cd34
}else if(sort=="MNC"){
  depth <- sample.info[patient.id,"Coverage.WGS.BM.MNC"]
  snvs <- snvs.mnc
}else if(sort=="MNC_minus_T"){
  depth <- sample.info[patient.id,"Coverage.WGS.BM.MNC-T"]
  snvs <- snvs.mnc_minus_t
}else if(sort=="PB_gran"){
  depth <- sample.info[patient.id,"Coverage.WGS.PB.granulocytes"]
  snvs <- snvs.pb_gran
}else if(sort=="CD34_deep"){
  depth <- 270
  snvs <- snvs.cd34.deep
  min.vaf <- 0.02
  min.prior.size <- 0.001
  min.clone.size <- 0.01
}

age <- sample.info[patient.id,]$Age*365

snvs <- list(snvs[[patient.id]])

use.sensitivity = F

source("~/Parameter_estimation/Bayesian_fit_no_size_compensation.R")
