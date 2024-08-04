####### load libraries
library(SCIFER)
library(data.table)
library(openxlsx)
###### load functions

patient.id <- 'KX004' # specify patient

age <- 77*365 # age in days

depth=10000 # sequencing depth (arbitrary, as we run in single-cell mode)

load("./RData/Mitchell_et_al/SNVs.RData")

snvs <- list(data.frame(VAF=snvs[[patient.id]], Depth=100, varCounts=snvs[[patient.id]]*100))


min.vaf <- 0.02
min.clone.size = 0.01
min.prior.size=0.001

use.sensitivity = F
seq.type="sc"
ncells=451 # number of analyzed cells
mother.daughter <- matrix(c(1, 2, 1,3), byrow=T, ncol=2) # branched clonal topology

source("~/Parameter_estimation/Bayesian_fit_multiclone.R")

