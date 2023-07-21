rm(list=ls())
####### load libraries
library(deSolve)
library(openxlsx)
library(FLORENCE)
###### load functions

patient.id <- 'AX001' # specify patient ID here

age <- 63*365 # age in days

depth=10000 # sequencing depth (arbitrary, as we run in single-cell mode)

load("./Mitchell_et_al/SNVs.RData")

snvs <- list(data.frame(VAF=snvs[[patient.id]], Depth=100, varCounts=snvs[[patient.id]]*100)) # set arbitrary Depth of 100x (doesn't matter, as the information is later only used to filter variants with very small or very high coverage, what was, however, already done for this data).

use.sensitivity <- F

ncells = 361 # number of sequenced cells

seq.type="sc" # parameter estimation mode: single-cell

source("~/Parameter_estimation/Bayesian_fit.R")
