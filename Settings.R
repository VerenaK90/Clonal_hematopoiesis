############################################################################################################################################
### load libraries
library(openxlsx)
library(wesanderson)
library(bedr)
library(ggplot2); theme_set(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=8, color="black"),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.75),
                                  axis.ticks = element_line(color = "black", size=0.75),
                                  axis.text = element_text(size=8, color="black")))
library(ggpubr)
library(ggformula)
library(ggsci)
library(deSolve)
library(HDInterval)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(ggbeeswarm)
library(Hmisc)
library(ggbeeswarm)
library(lemon)
library(data.table)
library(SCIFER)

##############################################################################################################################################
## set directories

custom.script.directory <- "./"
analysis.directory <- "./"
dir.create(paste(analysis.directory, "Figures", sep="/"))
meta.data <- "./MetaData/"
rdata.directory <- "./RData/"

############################################################################################################################################
### read in sample information

sample.info <- read.xlsx("Metadata/Supplementary Tables.xlsx", sheet = 2, startRow = 6)
rownames(sample.info) <- sample.info$Paper_ID
sample.info <- sample.info[order(sample.info$Age),]

sample.info$CHIP.mutation.associated.with.fit <- sapply(sample.info$Paper_ID, function(x){
  if(grepl("N", x)){
    "no driver"
  }else if(grepl("A", x)){
    "ASXL1"
  }else if(grepl("D", x)){
    "DNMT3A"
  }else if(grepl("T", x)){
    "TET2"
  }else if(grepl("U", x)){
    "unknown driver"
  }
})

samples <- rownames(sample.info)
normal.samples <- samples[sample.info$CHIP.mutation.associated.with.fit=="no driver"]

chip.samples.unknown.driver <- samples[sample.info$CHIP.mutation.associated.with.fit=="unknown driver"]

chip.samples <- samples[sample.info$CHIP.mutation.associated.with.fit!="no driver" &
                          sample.info$CHIP.mutation.associated.with.fit != "unknown driver"]

############################################################################################################################################
### read in sample information for published data

sample.info.published.data <- read.xlsx(paste0(meta.data, "Sample_info_published_data.xlsx"), sheet = 1)

############################################################################################################################################
### define color schemes

sample.color <- brewer.pal(4, "Dark2")
names(sample.color) <- c("MNC", "PB_gran", "CD34", "MNC_minus_T")  

CHIP.color <- c("darkgrey", "black", wes_palette("Darjeeling1", 4, type = "discrete"))
names(CHIP.color) <- c("no driver", "unknown driver", "TET2", "DNMT3A", "ASXL1", "SF3B1")

model.colors <- c(neutral = '#7891F9', selection = '#F79292')


############################################################################################################################################
## custom function to adjust hdi interval size for plotting density distributions; courtesy to https://stackoverflow.com/questions/59471972/shade-an-area-under-density-curve-to-mark-the-highest-density-interval-hdi

hdi_custWidth <- function(...) {
  dots <- list(...)
  quantiles <- dots[[2]]
  hdi_width <- quantiles[[length(quantiles)]] # uses the last entry if its a vector which should be the biggest one; better pass a single double < 1.0
  if (is.na(hdi_width)) hdi_width <- .89 # happens is quantiles = 1L
  message(paste0('HDI credible interval width = ', hdi_width))
  HDInterval::hdi(dots[[1]], credMass = hdi_width)
}

