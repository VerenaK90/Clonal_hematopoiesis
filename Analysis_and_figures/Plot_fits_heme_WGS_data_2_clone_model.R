###### This scripts plots the output from the population genetics model for 2 clones. It produces
###### - For each patient the model fit vs data and the posterior probabilities of the parameters
###### - A summary of the 80%HDI estimates across the population
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### libraries and functions

library(cdata)
library(scales)
library(openxlsx)
library(ggridges)

source("./Settings.R")

## mutation data
load("RData/WGS_heme/SNVs.RData")

putative.drivers <- read.xlsx("MetaData/Supplementary Tables.xlsx", startRow = 8, sheet = 7)
sample.info <- read.xlsx("MetaData/Supplementary Tables.xlsx", sheet = 2, startRow = 6)
rownames(sample.info) <- sample.info$Paper_ID
patient.ids <- sample.info$Paper_ID[sample.info$`Coverage.WGS.CD34+.2`==180]

############################################################################################################################################
####### Set parameters as used for model fits
depth=270
min.vaf <- 0.02
min.clone.size = 0.01
min.prior.size=0.001

use.sensitivity <- F
sample.color["CD34+"] <- sample.color["CD34"]
seq.type <- "bulk"

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Plot per patient and compare parameter across cohort

# collect highest density intervals of the parameters:
# .. overall ..
parameters.2cm <- c()
# .. neutral fits only ..
neutral.parameters.2cm <- c()
# .. selection fits supporting the first selected clone ... 
selected.parameters.2cm.1clone <- c()
# .. selection fits supporting the second selected clone 
selected.parameters.2cm.2clones <- c()

# collect model support for selection of the first and second clone (%)
model.support.selection.2_clones <- matrix(NA, nrow = 2, ncol = length(patient.ids)*2,
                                dimnames = list(c("Clone 1", "Clone 2"), 
                                                c(paste0(patient.ids, "_branched"),
                                                  paste0(patient.ids, "_linear"))))

## store model fits in a list
plotlist.model.vs.data.2cm <- list()


for(patient.id in patient.ids){
  
  snvs <- list(snvs.cd34.deep[[patient.id]])
  age <- sample.info[patient.id,]$Age*365
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                            colnames(putative.drivers)[grepl("CD34", colnames(putative.drivers)) &
                                                                         grepl("deep", colnames(putative.drivers))])]
  
  driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(patient.id, "_"), colnames(driver.information))])]
  if(any(grepl(paste0(patient.id, "_"), colnames(driver.information)))){
    driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", type="snvs", mutationcaller="mpileup", sample.col.mpileup = colnames(driver.information)[grep(paste0(patient.id, "_"), colnames(driver.information))])
    driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", type="snvs", mutationcaller="mpileup", sample.col.mpileup = colnames(driver.information)[grep(paste0(patient.id, "_"), colnames(driver.information))])
    driver.information$nvar <- driver.information$VAF*driver.information$Depth
    driver.information <- driver.information[driver.information$nvar>2 &
                                               driver.information$Depth > 1/3*depth,,drop=F]
    driver.information$lower <- driver.information$VAF-1.96*sqrt(driver.information$VAF*(1-driver.information$VAF)/driver.information$Depth) ## Wald approximation for 95% binomial CI
    driver.information$upper <- driver.information$VAF+1.96*sqrt(driver.information$VAF*(1-driver.information$VAF)/driver.information$Depth) ## Wald approximation for 95% binomial CI
    driver.information$lower[driver.information$lower<0] <- 0
    driver.information$upper[driver.information$upper>1] <- 1
    driver.information <- driver.information[!is.na(driver.information$VAF),]
  }else{
    driver.information <- data.frame()
  }
  
  ## subset on drivers in "TET2", "DNMT3A", "ASXL1" and, in T2, KMT2D - the others are putative drivers with unknown consequences
  driver.information <- driver.information[driver.information$GENE %in% c("DNMT3A", "TET2", "ASXL1") |
                                             (driver.information$GENE == "KMT2D" & patient.id == "T2"),]
  
  ## for each sample we fit 2-clone models with phylogenetic topologies according to linear and branched evolution
  for(mode in c("branched", "linear")){
    ##########################################################################################################
    ## read in fits and analyze the parameters
    
    if(mode == "branched"){
      
      fits <- read.csv(paste0(analysis.directory, "Model_fits_2_clones/WGS_heme/", patient.id, "/Model_fit_branched.csv"))
      # mother-daughter relationship for branched evolution: both selected clones (2, 3) stem from the normal cells (1)
      mother.daughter <- matrix(c(1, 2, 1,3), byrow=T, ncol=2)
      
    }else if(mode =="linear"){
      
      fits <- read.csv(paste0(analysis.directory, "/Model_fits_2_clones/WGS_heme/", patient.id, "/Model_fit_linear.csv"))
      # mother-daughter relationship for linear evolution: the second selected clone (3) stems from the first selected clone (2), which stems from normal cells (1)
      mother.daughter <- matrix(c(1, 2, 2,3), byrow=T, ncol=2)
      
    }
    
    # compute selective advantage from size of the selected clone
    fits$s1 <- 1 - log(10^fits$par_size1*10^fits$par_N)/(10^fits$par_lambda_ss*(age - fits$par_ts1))
    fits$s2 <- 1 - log(10^fits$par_size2*10^fits$par_N)/(10^fits$par_lambda_ss*(age - fits$par_ts2))
    
    fits$mutations_per_year <- fits$par_mu*10^fits$par_lambda_ss*365
    
    # compute clone sizes
    clone_sizes <- apply(fits, 1, function(x){
      res <- SCIFER:::.compute_actual_size(t.s = c(0, as.numeric(x["par_ts1"]), as.numeric(x["par_ts2"])), mother.daughter = mother.daughter, 
                                           N = 10^as.numeric(x["par_N"]), lambda.ss = 10^as.numeric(x["par_lambda_ss"]), 
                                           delta.ss = 10^as.numeric(x["par_lambda_ss"]),
                                           lambda.exp = 1, delta.exp = as.numeric(x["par_delta_exp"]), 
                                           size = c(10^as.numeric(x["par_size1"]), 10^as.numeric(x["par_size2"])), t.end = age)
      res[-1]/10^as.numeric(x["par_N"])
    })
    
    fits$clone_size_1 <- clone_sizes[1,]
    fits$clone_size_2 <- clone_sizes[2,]
    fits$clone_size_1_log10 <- log10(clone_sizes[1,])
    fits$clone_size_2_log10 <- log10(clone_sizes[2,])
    fits$age_of_clone_1 <- (age - fits$par_ts1)/365
    fits$age_of_clone_2 <- (age - fits$par_ts2)/365
    
    # to compute the average growth per year from the final clone size
    # to do this, go from 1 cell to the final clone size with exponential growth
    fits$growth_per_year1 <- apply(fits, 1, function(x){
      # effective selective advantage
      seff = 1 - log(as.numeric(x["clone_size_1"])*10^as.numeric(x["par_N"]))/(10^as.numeric(x["par_lambda_ss"])*(age - as.numeric(x["par_ts1"])))
      # effective growth rate per year
      exp(10^as.numeric(x["par_lambda_ss"])*(1 - seff)*365) - 1
    })
    fits$growth_per_year2 <- apply(fits, 1, function(x){
      # effective selective advantage
      seff = 1 - log(as.numeric(x["clone_size_2"])*10^as.numeric(x["par_N"]))/(10^as.numeric(x["par_lambda_ss"])*(age - as.numeric(x["par_ts2"])))
      # effective growth rate per year
      exp(10^as.numeric(x["par_lambda_ss"])*(1 - seff)*365) - 1
    })
    
    
    ## sort clones by size
    order.fits <- t(apply(fits, 1, function(x){
      order(x[c("clone_size_1", "clone_size_2")], decreasing = T)
    }))
    
    fits[,c("par_size1", "par_size2")] <- t(apply(cbind(fits[,c("par_size1", "par_size2")], order.fits), 1, function(x){
      x[c(1,2)][x[c(3,4)]]
    }))
    fits[,c("s1", "s2")] <- t(apply(cbind(fits[,c("s1", "s2")], order.fits), 1, function(x){
      x[c(1,2)][x[c(3,4)]]
    }))
    fits[,c("age_of_clone_1", "age_of_clone_2")] <- t(apply(cbind(fits[,c("age_of_clone_1", "age_of_clone_2")], order.fits), 1, function(x){
      x[c(1,2)][x[c(3,4)]]
    }))
    fits[,c("clone_size_1", "clone_size_2")] <- t(apply(cbind(fits[,c("clone_size_1", "clone_size_2")], order.fits), 1, function(x){
      x[c(1,2)][x[c(3,4)]]
    }))
    fits[,c("par_ts1", "par_ts2")] <- t(apply(cbind(fits[,c("par_ts1", "par_ts2")], order.fits), 1, function(x){
      x[c(1,2)][x[c(3,4)]]
    }))
    fits[,c("growth_per_year1", "growth_per_year2")] <- t(apply(cbind(fits[,c("growth_per_year1", "growth_per_year2")], order.fits), 1, function(x){
      x[c(1,2)][x[c(3,4)]]
    }))
    
    # store the model support for the first and second selected clone, where first and second are defined based on size
    model.support.selection.2_clones["Clone 1", paste(patient.id, mode, sep="_")] <- sum(fits$clone_size_1 >= 2*min.vaf )/nrow(fits)*100
    model.support.selection.2_clones["Clone 2", paste(patient.id, mode, sep="_")] <- sum(fits$clone_size_2 >= 2*min.vaf)/nrow(fits)*100
    
    ##########################################################################################################
    ## plot model parameters
    
    ## compute Nxtau
    
    fits$N_tau <- 10^fits$par_N / (365 * 10^fits$par_lambda_ss)
    
    modelResult <- list()
    
    to.plot <- fits[,c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "N_tau", "par_offset", "par_ts1",
                       "par_ts2", "s1", "s2", "clone_size_1", "clone_size_2", "growth_per_year1", "growth_per_year2",
                       "mutations_per_year")]
    
    # transform scaling
    to.plot$par_ts1 <- to.plot$par_ts1/365
    to.plot$par_ts2 <- to.plot$par_ts2/365
    to.plot$par_lambda_ss <- 10^to.plot$par_lambda_ss*365
    
    to.plot <- melt(to.plot)
    
    pdf(paste0(analysis.directory, "Model_fits_2_clones/WGS_heme/", patient.id, "/Parameter_estimates_", mode, ".pdf"), width = 6, height = 6)

    p <- ggplot(to.plot, aes(x=value, y= 0, fill = stat(quantile))) +
      geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
      #geom_density(fill="grey") +
      facet_wrap("variable", scales="free") +
      scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none")+
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")
    
    print(p)
    
    ## plot again with log-scale
    p <- ggplot(to.plot[!to.plot$variable %in% c("clone_size_1", "clone_size_2"),], aes(x=value, y= 0, fill = stat(quantile))) +
      geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
      #geom_density(fill="grey") +
      facet_wrap("variable", scales="free", nrow = 4, ncol = 4) +
      scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none")+
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability") + scale_x_log10()
    
    print(p)
    
    ## compute highest density intervals for all fits ... 
    hdinterval <- as.data.frame(t(hdi(fits, credMass=0.8)))
    hdinterval$Parameter <- rownames(hdinterval)
    hdinterval$Median <- apply(fits, 2, function(x){median(as.numeric(x))})
    hdinterval$Paper_ID <- patient.id
    hdinterval$Mode <- mode
    
    parameters.2cm <- rbind(parameters.2cm, hdinterval)
    
    ## ... for all fits supporting a neutral model ....
    hdinterval.neutral <- as.data.frame(t(hdi(fits[fits$clone_size_1 < 2*min.vaf &
                                                     fits$clone_size_2 < 2*min.vaf,], credMass=0.8)))
    hdinterval.neutral$Parameter <- rownames(hdinterval.neutral)
    hdinterval.neutral$Median <- apply(fits, 2, function(x){median(as.numeric(x))})
    hdinterval.neutral$Paper_ID <- patient.id
    hdinterval.neutral$Mode <- mode
    
    neutral.parameters.2cm <- rbind(neutral.parameters.2cm, hdinterval.neutral)
    
    ## ... and the selection parameters associated with situations where there's evidence for 1 selected clone or 2 selected clones
    fits.selection.1clone <- fits[(fits$clone_size_1 >= 2*min.vaf | fits$clone_size_2 >= 2*min.vaf) & 
                                    !(fits$clone_size_1 >= 2*min.vaf & fits$clone_size_2 >= 2*min.vaf),] ## fits supporting only 1 selected clone
    fits.selection.2clones <- fits[fits$clone_size_1 >= 2*min.vaf & fits$clone_size_2 >= 2*min.vaf,] ## fits supporting two selected clones
    
    hdinterval.selection.1clone <- as.data.frame(t(hdi(fits.selection.1clone, credMass=0.8)))
    hdinterval.selection.1clone$Parameter <- rownames(hdinterval.selection.1clone)
    hdinterval.selection.1clone$Median <- apply(fits.selection.1clone, 2, function(x){median(as.numeric(x))})
    hdinterval.selection.1clone$Paper_ID <- patient.id
    hdinterval.selection.1clone$Mode <- mode
    selected.parameters.2cm.1clone <- rbind(selected.parameters.2cm.1clone, hdinterval.selection.1clone)
    
    hdinterval.selection.2clones <- as.data.frame(t(hdi(fits.selection.2clones, credMass=0.8)))
    hdinterval.selection.2clones$Parameter <- rownames(hdinterval.selection.2clones)
    hdinterval.selection.2clones$Median <- apply(fits.selection.2clones, 2, function(x){median(as.numeric(x))})
    hdinterval.selection.2clones$Paper_ID <- patient.id
    hdinterval.selection.2clones$Mode <- mode
    selected.parameters.2cm.2clones <- rbind(selected.parameters.2cm.2clones, hdinterval.selection.2clones)
    
    ### 2D-correlations
    
    # physiological parameters
    
    p1 <- ggplot(fits, aes(x=par_mu, y=par_N)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_continuous(name="mu", limits=c(0.1, 10)) + scale_y_continuous(name="N (log10)", limits = c(2.5, 8))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    p2 <- ggplot(fits, aes(x=par_lambda_ss, y=par_N)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_continuous(name="lambda (log10)", limits=c(-3, -1)) + scale_y_continuous(name="N (log10)", limits = c(2.5, 8))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    
    p3 <- ggplot(fits, aes(x=par_mu, y=par_lambda_ss)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_continuous(name="mu", limits=c(0.1, 10)) + scale_y_continuous(name="lambda (log10)", limits = c(-3, -1))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    
    p4 <- ggplot(fits, aes(x=10^par_N/10^par_lambda_ss, y=par_mu)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_log10(name="N/lambda", limits=c(10^2.5/10^-1, 10^11)) + scale_y_continuous(name="mu", limits = c(0.1, 10))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    
    p5 <- ggplot(fits, aes(x=par_mu, y=par_delta_exp)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_continuous(name="mu", limits=c(0.1, 10)) + scale_y_continuous(name="delta", limits = c(0, 0.75))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    p6 <- ggplot(fits, aes(x=par_N, y=par_delta_exp)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_continuous(name="N (log10)", limits=c(2.5, 8)) + scale_y_continuous(name="delta", limits = c(0, 0.75))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    p7 <- ggplot(fits, aes(x=par_lambda_ss, y=par_delta_exp)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_continuous(name="lambda (1/d; log10)", limits=c(-3, -1)) + scale_y_continuous(name="delta", limits = c(0, 0.75))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    p8 <- ggplot(fits, aes(x=par_mu*10^par_N, y=par_delta_exp)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_log10(name="mu*N", limits=c(0.1*10^2.5, 10*10^8)) + scale_y_continuous(name="delta", limits = c(0, 0.75))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    p9 <- ggplot(fits, aes(x=10^par_N, y=clone_size_1)) +
      geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_log10(name="N", limits=c(10^2.5, 10^8)) + scale_y_continuous(name="size1", limits = c(0.01, 1))+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
    
    print(ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow=3, ncol=3))
    
    
    
    parms.of.interest <- c("par_N", "par_lambda_ss", "par_delta_exp", "par_mu",
                           "clone_size_1_log10", "clone_size_2_log10", "par_ts1", "par_ts2")
    
    parameters.to.plot <- data.frame(parms_i = rep(parms.of.interest, each = length(parms.of.interest)),
                                     parms_j = rep(parms.of.interest, length(parms.of.interest)),
                                     i = rep(1:length(parms.of.interest), each = length(parms.of.interest)),
                                     j = rep(1:length(parms.of.interest), length(parms.of.interest)))
    
    plot.design <- data.frame(Parameter = parms.of.interest,
                              Min = c(2.5, -3, 0, 0.1, -3, -3, 0, 0),
                              Max = c(8, -1, 0.75, 10, 0, 0, age, age),
                              Label = c("N (log10)", "lambda_ss (log10)", "delta/lambda",
                                        "mu (1/div)", "size of clone 1 (log10)", "size of clone2 (log10)",
                                        "age of clone1 (years)", "age of clone2 (years)"))
    
    
    # fill up
    p <- list()
    
    for(ii in 1:length(parms.of.interest)){
      for(j in 1:length(parms.of.interest)){
        if(ii == j){
          p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits[!is.infinite(fits[,parms.of.interest[ii]]) &
                                                                     !is.infinite(fits[,parms.of.interest[[j]]]),], aes(x=.data[[parms.of.interest[ii]]])) +
            geom_histogram() + scale_x_continuous(limits=unlist(plot.design[ii,c("Min", "Max")]),
                                                  name = unlist(plot.design[ii, "Label"])) + 
            scale_y_continuous() +
            theme(axis.text.y=element_blank(),  #remove y axis labels
                  axis.ticks.y=element_blank(),  #remove y axis ticks
                  axis.title.y = element_blank()
            )
        }else if(ii < j){
          p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits[!is.infinite(fits[,parms.of.interest[ii]]) &
                                                                     !is.infinite(fits[,parms.of.interest[[j]]]),], aes(x=.data[[parms.of.interest[ii]]],
                                                                                                                        y=.data[[parms.of.interest[j]]])) +
            geom_point() + scale_x_continuous(limits=unlist(plot.design[ii,c("Min", "Max")]),
                                              name = unlist(plot.design[ii, "Label"])) + 
            scale_y_continuous(limits=unlist(plot.design[j,c("Min", "Max")]),
                               name = unlist(plot.design[j, "Label"]))
        }else{
          p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits[!is.infinite(fits[,parms.of.interest[ii]]) &
                                                                     !is.infinite(fits[,parms.of.interest[[j]]]),], aes(x=.data[[parms.of.interest[ii]]],
                                                                                                                        y=.data[[parms.of.interest[j]]])) +
            geom_density2d() + scale_x_continuous(limits=unlist(plot.design[ii,c("Min", "Max")]),
                                                  name = unlist(plot.design[ii, "Label"])) + 
            scale_y_continuous(limits=unlist(plot.design[j,c("Min", "Max")]),
                               name = unlist(plot.design[j, "Label"]))
        }
        if(j!=length(parms.of.interest)){
          p[[length(parms.of.interest)*(ii-1) + j]] <-  p[[length(parms.of.interest)*(ii-1) + j]] + theme(axis.text.x=element_blank(),
                                                                                                          axis.ticks.x=element_blank(),
                                                                                                          axis.title.x = element_blank())
        }
        if(ii!=1){
          p[[length(parms.of.interest)*(ii-1) + j]] <- p[[length(parms.of.interest)*(ii-1) + j]] + theme(axis.text.y=element_blank(),
                                                                                                         axis.ticks.y=element_blank(),
                                                                                                         axis.title.y = element_blank())
        }
      }
    }
    
    print(egg::ggarrange(plots = p, nrow=length(parms.of.interest), byrow=F)) 
    
    p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
      facet_wrap(~Parameter, scales="free")
    
    print(p)
    
    
    dev.off()
    
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
    ###### Plot fits for neutral and selected case
   
    source(paste0(custom.script.directory, "/Parameter_estimation/Bayesian_fit_multiclone.R"))
    
    # simulate the trajectories for 100 model fits
    if(!file.exists(paste0(analysis.directory, "Model_fits_2_clones/WGS_heme/", patient.id, "/Sim_trajectories_", mode, ".RData"))){
      sim <- matrix(0, nrow=100, ncol=length(mySumStatData$mutation.count[[1]]))
      
      for(j in 1:100){
        print(j)
        
        parms <- list(mu=fits$par_mu[j], N=fits$par_N[j], delta_exp = fits$par_delta_exp[j], lambda_ss=fits$par_lambda_ss[j],
                      offset=fits$par_offset[j], ts1=fits$par_ts1[j], ts2 = fits$par_ts2[j],
                      size1 = fits$par_size1[j], size2 = fits$par_size2[j])
        
        model <- myModel(parms)
        sim[j,] <- model$modelResult[[1]]
        
      }
      
      max.pred <- apply(sim, 2, quantile, p=0.975)
      min.pred <- apply(sim, 2, quantile, p=0.025)
      
      data.vs.prediction <- data.frame(VAF=rep(vafs.of.interest), mean.data=mySumStatData$mutation.count[[1]],
                                       sd.data = mySumStatData$sampled.sd[[1]],
                                       min.model = min.pred, max.model=max.pred,
                                       Age=mySumStatData$age/365)
      
      save(sim, data.vs.prediction, file=paste0(analysis.directory, "Model_fits_2_clones/WGS_heme/", patient.id, "/Sim_trajectories_", mode, ".RData"))
      
    }else{
      load(paste0(analysis.directory, "Model_fits_2_clones/WGS_heme/", patient.id, "/Sim_trajectories_", mode, ".RData"))
    }
    
    to.plot <- data.vs.prediction
    
    grDevices::pdf(paste0(analysis.directory, "/Model_fit2_2_clones/WGS_heme/", patient.id, "/Model_fit_", mode, ".pdf"), width=3, height=2.5)
    
    max.y <- max(to.plot$max.model)
    
    p <- ggplot(data=to.plot, aes(x=VAF, y=mean.data, ymin=mean.data-sd.data, ymax=mean.data+sd.data)) +
      geom_ribbon(data=to.plot, aes(x=VAF, y=mean.data, ymin=min.model,ymax=max.model), alpha=1, fill="violet") + 
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + scale_y_continuous(name="Cumulative # of mutations") +
      scale_x_continuous(limits=c(0, 0.6))
    
    if(nrow(driver.information)>0){
      p <- p + geom_segment(data = driver.information, aes(x = as.numeric(VAF), y = -250, xend = as.numeric(VAF), yend = 0),
                            arrow = arrow(length = unit(0.5, "cm")), col="firebrick", inherit.aes = F) +
        geom_text(data = driver.information, aes( x=as.numeric(VAF), y=-100), label=driver.information[,'GENE'],  color="firebrick",
                  size=5 ,  fontface="italic", inherit.aes = F)
    }
    
    print(p)
    
    
    
    xlimits <- c(0,1/min.vaf)
    
    p <- ggplot(data=to.plot, aes(x=1/VAF, y=mean.data, ymin=mean.data-sd.data, ymax=mean.data+sd.data)) +
      geom_ribbon(data=to.plot, aes(x=1/VAF, y=mean.data, ymin=min.model,ymax=max.model), alpha=1,
                  fill="violet") +
      geom_pointrange(lwd=0.25, shape=1, fatten=1) +
      scale_y_continuous(name="Cumulative # of mutations") +
      theme(aspect.ratio = min.vaf/0.05, legend.position = "bottom") + guides( col = guide_legend(title="Mutation", ncol = 4)) +
      scale_x_continuous( breaks = c(5, 10, 20, 50, 100), labels = c("0.2", "0.1", "0.05", "0.02", "0.01"), name = "Variant allele frequency") +
      coord_cartesian(ylim=c(0, max.y*1.05), xlim=xlimits)
    
    if(nrow(driver.information[driver.information$upper>=1/max(xlimits),])>0){
      
      driver.information$y <- seq(max(to.plot$max.model/10), max(to.plot$max.model), length.out = nrow(driver.information))
      
      p <- p + geom_point(data = driver.information[driver.information$upper>=1/max(xlimits),], aes(x = 1/as.numeric(VAF), y = y, col = paste(GENE, AAchange)),  inherit.aes = F) +
        geom_errorbarh(data =  driver.information[driver.information$upper>=1/max(xlimits),], aes(xmin = 1/as.numeric(upper), y = y,
                                                                                                  xmax = 1/as.numeric(lower), col =  paste(GENE, AAchange)), height = max.y/50,  inherit.aes = F) +
        theme(legend.text=element_text(size=6))
      
    }
    
    
    if(model.support.selection.2_clones[1,paste(patient.id, mode, sep="_")] > 15){
      p <- p + geom_ribbon(data = data.frame(x = unlist(hdinterval.selection.1[hdinterval.selection.1$Parameter=="clone_size_1",c("lower", "upper")]),
                                             ymin = c(0,0),
                                             ymax = c(max.y, max.y)), aes(x=2/x, ymin=ymin, ymax = ymax), inherit.aes = F, fill="darkblue", alpha = 0.5)
      
    }
    
    if(model.support.selection.2_clones[2,paste(patient.id, mode, sep="_")] > 15){
      p <- p + geom_ribbon(data = data.frame(x = unlist(hdinterval.selection.2[hdinterval.selection.2$Parameter=="clone_size_2",c("lower", "upper")]),
                                             ymin = c(0,0),
                                             ymax = c(max.y, max.y)), aes(x=2/x, ymin=ymin, ymax = ymax), inherit.aes = F, fill="darkgreen", alpha = 0.5)
    }
    
    print(p)
    
    plotlist.model.vs.data.2cm[[paste(patient.id, mode, sep="_")]] <- p
    
    dev.off()
    
  }
}

####################################################################################################################################################
## Figures 6g,i, Extended Data. Fig. 8e,f: plot the model fits 

pdf(paste0(analysis.directory, "/Figures/Figure_6gi_S8ef.pdf"), width=8, height=6)

ggarrange(plotlist=plotlist.model.vs.data[c("A2_linear", "T3_linear", "U6_branched",
                                            "A1_branched", "D1_branched", "D2_linear",
                                            "D4_branched", "U2_branched", "U3_branched")],
  nrow=3, ncol=3)

dev.off()

####################################################################################################################################################
## Figure 6h plot the model support for the second clone 

to.plot <- melt(t(model.support.selection.2_clones), value.name = "P_selection")
colnames(to.plot)[c(1,2)] <- c("Fit_ID", "Clone")
to.plot$ID <- gsub("_.*", "", to.plot$Fit_ID)
## focus on the second clone and take the topology (linear/branched) with maximum support for second clone per ID
to.plot <- to.plot[to.plot$Clone == "Clone 2",]
to.plot <- to.plot[order(to.plot$P_selection, decreasing = T),]
to.plot <- to.plot[order(to.plot$ID),]
to.plot <- to.plot[!duplicated(to.plot$ID),]
to.plot$P_neutral <- 100 - to.plot$P_selection
to.plot <- melt(to.plot, value.name = "Posterior probability", id.vars = c("Fit_ID", "Clone", "ID"),
                measure.vars = c("P_selection", "P_neutral"), variable.name = "Type")
to.plot$Clone_size <- apply(to.plot, 1, function(x){
  mode <- gsub(".*_", "", x[1])
  res <- selected.parameters.2cm.2clones[selected.parameters.2cm.2clones$Paper_ID == x[3] &
                                       selected.parameters.2cm.2clones$Mode == mode &
                                       selected.parameters.2cm.2clones$Parameter=="true_size_2" &
                                       x[4] == "P_selection",]$Median
  if(length(res)==0){
    return(0)
  }else{
    return(res)
  }
})
to.plot$Type <- factor(to.plot$Type, levels=c("P_neutral", "P_selection"))
to.plot$Clone_size[to.plot$Clone_size==0] <- NA


pdf(paste0(analysis.directory, "/Figures/Figure_6h.pdf"), width=6, height=3.5)

to.plot.selected <- to.plot[ !is.na(to.plot$`Posterior probability`) & !grepl("N", to.plot$ID),]

ggplot(to.plot.selected,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") + 
  scale_x_discrete(breaks = levels(to.plot.selected$ID), limits = levels(to.plot.selected$ID)) +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, selection associated with 2nd CH driver")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

dev.off()

