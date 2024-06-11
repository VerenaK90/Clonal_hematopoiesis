###### This scripts plots the output from the population genetics model for a one-clone-model without size compensation. It produces
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
patient.ids <- sample.info$Paper_ID[grepl("T", sample.info$Paper_ID)]

############################################################################################################################################
####### Set parameters as used for model fits
depth <- 90
min.vaf <- 0.05
min.clone.size <- 0.05
min.prior.size <- 0.01

use.sensitivity <- F
sample.color["CD34+"] <- sample.color["CD34"]
seq.type <- "bulk"

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Plot per patient and compare parameter across cohort

# collect highest density intervals of the parameters:
# .. overall ..
parameters.nsc <- data.frame()
# .. neutral fits only ..
neutral.parameters.nsc <- data.frame()
# .. selection fits supporting the selected clone 
selected.parameters.nsc <- data.frame()

# collect model support for selection (%)
model.support.selection.nsc <- matrix(NA, nrow=1, ncol=length(patient.ids), dimnames = list(c("CD34_wo_size_compensation"),
                                                                                      patient.ids))
plotlist.model.vs.data <- list()

for(patient.id in patient.ids){
  print("....")
  print(patient.id)
  age <- sample.info[patient.id, ]$Age*365
  
  snvs <- list(snvs.cd34[[patient.id]])
  
  ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
  ###### load observed data
  
  directory <- paste0(analysis.directory, "/Model_fits/WGS_heme/", paste(patient.id, sort, "nsc", sep="_"))
  if(!file.exists(paste0(directory,  "/Model_fit.csv"))){next}
  fits <- read.csv(paste0(directory,  "/Model_fit.csv"))

  sort.type = "CD34+"
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                            colnames(putative.drivers)[grep(sort.type, colnames(putative.drivers))])]

  driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(patient.id, "_"), colnames(driver.information))])]
  driver.information <- driver.information[,colnames(driver.information)[grep( "_CD38", colnames(driver.information), invert = T)]]
  driver.information <- driver.information[,colnames(driver.information)[grep( "deep", colnames(driver.information), invert = T)]]
  
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
  }else{
    driver.information <- data.frame()
  }
  
  
  ## subset on drivers in "TET2", "DNMT3A", "ASXL1" and, in T2, KMT2D - the others are putative drivers with unknown consequences
  driver.information <- driver.information[driver.information$GENE %in% c("DNMT3A", "TET2", "ASXL1") |
                                             (driver.information$GENE == "KMT2D" & patient.id == "T2"),]
  
  ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
  #### Print parameter estimates
  
  pdf(paste0(directory, "/Parameter_estimates.pdf"), width=6, height=6)
  
  fits$par_t_s_absolute <- apply(fits, 1, function(x){
    min.t.s <- 0
    max.t.s <- age - log(min.prior.size*10^as.numeric(x["par_N"]))/10^as.numeric(x["par_lambda_ss"])
    
    t.s <- min.t.s +as.numeric(x["par_t_s"])*(max.t.s - min.t.s)
  })
  
  
  fits$par_s_absolute <- apply(fits, 1, function(x){
    ts <- as.numeric(x["par_t_s_absolute"])
    
    min.s <- (10^as.numeric(x["par_lambda_ss"])*(age-ts) - log(10*10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*
                                                                                                 (age-ts))
    if(min.s < 0){
      min.s <- 0
    }
    max.s <- (10^as.numeric(x["par_lambda_ss"])*(age-ts) - log(min.prior.size*10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*(age-ts))
    
    s <- min.s + as.numeric(x["par_s"])*(max.s-min.s)
    
  })
  
  
  fits$size_of_clone <- exp(10^fits$par_lambda_ss*(1-fits$par_s_absolute)*(age - fits$par_t_s_absolute))/(10^fits$par_N + exp(10^fits$par_lambda_ss*(1-fits$par_s_absolute)*(age - fits$par_t_s_absolute)))
  
  fits$size_of_clone[fits$size_of_clone>1] <- 1
  
  fits$age_of_clone <- (age - fits$par_t_s_absolute)/365
  
  ## how long did the clone grow until reaching 2% of the stem cell compartment?
  
  fits$age_of_clone_at_4pct <- log(0.04*10^fits$par_N)/(10^fits$par_lambda_ss*(1-fits$par_s_absolute))/365
  
  fits$growth_per_year <- exp(10^fits$par_lambda_ss*(1-fits$par_s_absolute)*365)-1
  
  fits$mutations_per_year <- fits$par_mu*10^fits$par_lambda_ss*365
  
  ## compute Nxtau
  
  fits$N_tau <- 10^fits$par_N / (365 * 10^fits$par_lambda_ss)
  
  to.plot <- fits[,c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "N_tau", "par_offset", "par_t_s_absolute",
                     "par_s_absolute", "size_of_clone", "growth_per_year",
                     "mutations_per_year", "age_of_clone")]
  
  to.plot$par_t_s_absolute <- to.plot$par_t_s_absolute/365
  to.plot$par_lambda_ss <- 10^to.plot$par_lambda_ss*365
  
  to.plot <- melt(to.plot)
  
  p <- ggplot(to.plot, aes(x=value, y= 0, fill = stat(quantile))) +
    geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
    #geom_density(fill="grey") +
    facet_wrap("variable", scales="free") +
    scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none")+
    scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")
  
  print(p)
  
  ## plot again with log-scale
  p <- ggplot(to.plot, aes(x=value, y= 0, fill = stat(quantile))) +
    geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
    #geom_density(fill="grey") +
    facet_wrap("variable", scales="free", nrow = 3, ncol = 4) +
    scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none")+
    scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability") + scale_x_log10()
  
  print(p)
  
  hdinterval <- as.data.frame(t(hdi(fits, credMass=0.8)))
  hdinterval$Parameter <- rownames(hdinterval)
  hdinterval$Median <- apply(fits, 2, function(x){median(as.numeric(x))})
  hdinterval$Paper_ID <- patient.id
  hdinterval$Sort <- "CD34"
  hdinterval$Depth <- depth
  
  
  ### 2D-correlations
  
  # physiolgical paramteres
  
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
  
  print(ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow=3, ncol=3))
  
  
  
  parms.of.interest <- c("par_N", "par_lambda_ss", "par_delta_exp", "par_mu", "age_of_clone", "par_s_absolute",
                         "size_of_clone")
  
  parameters.to.plot <- data.frame(parms_i = rep(parms.of.interest, each = length(parms.of.interest)),
                                   parms_j = rep(parms.of.interest, length(parms.of.interest)),
                                   i = rep(1:length(parms.of.interest), each = length(parms.of.interest)),
                                   j = rep(1:length(parms.of.interest), length(parms.of.interest)))
  
  plot.design <- data.frame(Parameter = parms.of.interest,
                            Min = c(2.5, -3, 0, 0.1, 0, 0, 0),
                            Max = c(8, -1, 0.75, 10, age, 0.99, 1),
                            Label = c("N (log10)", "lambda_ss (log10)", "delta/lambda",
                                      "mu (1/div)", "age of clone (years)", "s", "size of clone"))
  
  
  # fill up
  p <- list()
  
  for(ii in 1:length(parms.of.interest)){
    for(j in 1:length(parms.of.interest)){
      if(ii == j){
        p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits, aes(x=.data[[parms.of.interest[ii]]])) +
          geom_histogram() + scale_x_continuous(limits=unlist(plot.design[ii,c("Min", "Max")]),
                                                name = unlist(plot.design[ii, "Label"])) + 
          scale_y_continuous() +
          theme(axis.text.y=element_blank(),  #remove y axis labels
                axis.ticks.y=element_blank(),  #remove y axis ticks
                axis.title.y = element_blank()
          )
      }else if(ii < j){
        p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits, aes(x=.data[[parms.of.interest[ii]]],
                                                                      y=.data[[parms.of.interest[j]]])) +
          geom_point() + scale_x_continuous(limits=unlist(plot.design[ii,c("Min", "Max")]),
                                            name = unlist(plot.design[ii, "Label"])) + 
          scale_y_continuous(limits=unlist(plot.design[j,c("Min", "Max")]),
                             name = unlist(plot.design[j, "Label"]))
      }else{
        p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits, aes(x=.data[[parms.of.interest[ii]]],
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
  

  ## for selection associated values, store only the parameters associated with a clone >= 2*min.vaf
  
  fits.selection <- fits[fits$size_of_clone >= 2*min.vaf, ]
  fits.neutral <- fits[fits$size_of_clone < 2*min.vaf, ]
  
  if(nrow(fits.selection)>0){
    size.of.selected.clone <- median(fits$size_of_clone[fits$size_of_clone>=min.vaf])
    model.support.selection.nsc[1, patient.id] <- nrow(fits.selection)/10
    
  }else{
    model.support.selection.nsc[1,patient.id] <- 0
  }
  
  
  if(nrow(fits.selection)==0){
    
    parameters.nsc <- rbind(parameters.nsc, hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                 "mutations_per_year"),])
    
    neutral.parameters.nsc <- rbind(neutral.parameters.nsc,hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                                "mutations_per_year"),])
    
    p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
      facet_wrap(~Parameter, scales="free")
    
    print(p)
    
    
    dev.off()
    
  }else{
    hdinterval.selection <- as.data.frame(t(hdi(fits.selection, credMass=0.8)))
    hdinterval.selection$Parameter <- rownames(hdinterval.selection)
    hdinterval.selection$Median <- apply(fits.selection, 2, function(x){median(as.numeric(x))})
    hdinterval.selection$Paper_ID <- patient.id
    hdinterval.selection$Sort <- "CD34"
    hdinterval.selection$Depth <- depth
    
    
    p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
      facet_wrap(~Parameter, scales="free") + ggtitle("All parameters")
    
    print(p)
    
    p <- ggplot(data=hdinterval.selection, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
      facet_wrap(~Parameter, scales="free") + ggtitle("Selection")
    
    print(p)
    
    
    parameters.nsc <- rbind(parameters.nsc, hdinterval)
    
    selected.parameters.nsc <- rbind(selected.parameters.nsc, hdinterval.selection)
    
    if(nrow(fits.neutral)>0){
      hdi.neutral <- as.data.frame(t(hdi(fits.neutral, credMass = 0.8)))
      hdi.neutral$Parameter <- rownames(hdi.neutral)
      hdi.neutral$Median <- apply(fits.neutral, 2, function(x){median(as.numeric(x))})
      hdi.neutral$Paper_ID <- patient.id
      hdi.neutral$Sort <- "CD34"
      hdi.neutral$Depth <- depth
      
      neutral.parameters.nsc <- rbind(neutral.parameters.nsc, hdi.neutral[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                                    "mutations_per_year"),])
      
      p <- ggplot(data=hdi.neutral, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free") + ggtitle("Neutral")
      
      print(p)
      
    }
    
    dev.off()
  }
  
  
  ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
  ###### Plot fits for neutral and selected case
  
  source(paste0(custom.script.directory, "/Parameter_estimation/Bayesian_fit_no_size_compensation.R"))
  
  if( !file.exists(paste0(directory, "Sim_trajectories.RData") )){
    sim <- matrix(0, nrow=100, ncol=length(mySumStatData$mutation.count[[1]]))
    
    for(j in 1:100){
      print(j)
      
      parms <- list(mu=fits$par_mu[j], N=fits$par_N[j], delta_exp = fits$par_delta_exp[j], lambda_ss=fits$par_lambda_ss[j],
                    offset=fits$par_offset[j], t_s=fits$par_t_s[j], s=fits$par_s[j])
      
      model <- myModel(parms)
      sim[j,] <- model$modelResult[[1]]
      
    }
    
    max.pred <- apply(sim, 2, quantile, p=0.975)
    min.pred <- apply(sim, 2, quantile, p=0.025)
    
    best.fit <- which.min(fits$distance)
    
    best.parms <- list(mu=fits$par_mu[best.fit], N=fits$par_N[best.fit], delta_exp = fits$par_delta_exp[best.fit], lambda_ss=fits$par_lambda_ss[best.fit],
                       offset=fits$par_offset[best.fit], t_s=fits$par_t_s[best.fit], s=fits$par_s[best.fit])
    
    best.fit <- myModel(best.parms)
    best.fit <- best.fit$modelResult
    
    data.vs.prediction <- data.frame(VAF=rep(vafs.of.interest), mean.data=mySumStatData$mutation.count[[1]],
                                     sd.data = mySumStatData$sampled.sd[[1]],
                                     min.model = min.pred, max.model=max.pred,
                                     best.fit.model = best.fit[[1]],
                                     Age=mySumStatData$age/365)
    
    save(sim, data.vs.prediction, file=paste0(directory, "/Sim_trajectories.RData"))
    
  }else{
    load(paste0(directory, "/Sim_trajectories.RData"))
  }
  
  to.plot <- data.vs.prediction
  
  grDevices::pdf(paste0(directory, "/Model_fit.pdf"), width=3, height=2.5)
  
  max.y <- max(to.plot$max.model)
  
  xlimits <- c(0,1/min.vaf)
  
  as.ratio <- ifelse(min.vaf==0.05, 1, min.vaf/0.05)
  
  p <- ggplot(data=to.plot, aes(x=1/VAF, y=mean.data, ymin=mean.data-sd.data, ymax=mean.data+sd.data)) +
    geom_ribbon(data=to.plot, aes(x=1/VAF, y=mean.data, ymin=min.model,ymax=max.model), alpha=1,
                fill=sample.color[sort.type]) +
    ggtitle(paste(sample.info[patient.id,]$Paper_ID, sort.type))+
    geom_pointrange(lwd=0.25, shape=1, fatten=1) +
    scale_y_continuous(name="Cumulative # of mutations") +
    theme(aspect.ratio = as.ratio) +
    scale_x_continuous( breaks = c(5, 10, 20, 50, 100), labels = c("0.2", "0.1", "0.05", "0.02", "0.01"), name = "Variant allele frequency") +
    coord_cartesian(ylim=c(0, max.y*1.05), xlim=xlimits)
  
  if(nrow(driver.information[driver.information$upper>=1/max(xlimits),])>0){
    
    driver.information$y <- sapply(driver.information$VAF, function(x){
      to.plot[min(which(to.plot$VAF >= x)),]$max.model
    })
    
    p <- p + geom_point(data = driver.information[driver.information$upper>=1/max(xlimits),], aes(x = 1/as.numeric(VAF), y = y), col="firebrick", inherit.aes = F) +
      geom_errorbarh(data =  driver.information[driver.information$upper>=1/max(xlimits),], aes(xmin = 1/as.numeric(upper), y = y,
                                                                                                xmax = 1/as.numeric(lower)), height = max.y/50, col="firebrick", inherit.aes = F) +
      geom_text(data = driver.information[driver.information$upper>=1/max(xlimits),], aes( x=1/as.numeric(VAF), y=y*1.05), label=driver.information[driver.information$upper>=1/max(xlimits),'GENE'],  color="firebrick",
                size=6 ,  fontface="italic", inherit.aes = F)
  }
  
  if(nrow(fits.selection)>0){
    
    p <- p + geom_ribbon(data = data.frame(x = unlist(hdinterval.selection[hdinterval.selection$Parameter=="size_of_clone",c("lower", "upper")]),
                                           ymin = c(0,0),
                                           ymax = c(max.y, max.y)), aes(x=2/x, ymin=ymin, ymax = ymax), inherit.aes = F, fill="grey", alpha = 0.5)
  }
  
  print(p)
  plotlist.model.vs.data[[paste(patient.id, sort.type)]] <- p
  
  dev.off()
  
  
  
}

####################################################################################################################################################
## Extended Data. Fig. 7d: plot the model fit

pdf(paste0(analysis.directory, "/Figures/Figure_S7d.pdf"), width=3.5, height=3.5)

print(plotlist.model.vs.data[["T2"]])

dev.off()

