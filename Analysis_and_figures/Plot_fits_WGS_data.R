###### This scripts plots the output from the population genetics model. It produces
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
load("RData/WGS_data/SNVs.RData")

putative.drivers <- read.xlsx("MetaData/Supplementary Tables.xlsx", startRow = 8, sheet = 6)
sample.info <- read.xlsx("MetaData/Supplementary Tables.xlsx", sheet = 2, startRow = 6)
rownames(sample.info) <- sample.info$Paper_ID
patient.ids <- sample.info$Paper_ID
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

############################################################################################################################################
####### Set parameters as used for model fits
use.sensitivity <- F
sample.color["CD34+"] <- sample.color["CD34"]
seq.type <- "bulk"

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Plot per patient and compare parameter across cohort

# collect highest density intervals of the parameters:
# .. overall ..
parameters <- data.frame()
# .. neutral fist only ..
neutral.parameters <- data.frame()
# .. selection fits only
selected.parameters <- data.frame()

# collect model support for selection (%)
model.support.selection <- matrix(NA, nrow=4, ncol=length(patient.ids), dimnames = list(c("MNC", "CD34", "MNC_minus_T", "PB_gran"),
                                                                                      patient.ids))

plotlist.model.vs.data <- list()

for(patient.id in patient.ids){
  print(patient.id)
  age <- sample.info[patient.id, ]$Age*365
 
  cell.sorts <- c("CD34", "MNC", "MNC_minus_T", "PB_gran")
  
  for(tissue in cell.sorts){
    
    # specify coverage and snvs of this sample
    if(tissue == "CD34"){
      depth <- sample.info[patient.id,]$`Coverage.WGS.CD34+`
      if(depth==0){next}
      snvs <- list(snvs.cd34[[patient.id]])
    }else if(tissue == "MNC"){
      depth <- sample.info[patient.id,]$Coverage.WGS.BM.MNC
      if(depth==0){next}
      snvs <- list(snvs.mnc[[patient.id]])
    }else if (tissue =="MNC_minus_T"){
      depth <- sample.info[patient.id,]$`Coverage.WGS.BM.MNC–T`
      if(depth==0){next}
      snvs <- list(snvs.mnc_minus_t[[patient.id]])
    }else {
      depth <- sample.info[patient.id,]$Coverage.WGS.PB.granulocytes
      if(depth==0){next}
      snvs <- list(snvs.pb_gran[[patient.id]])
    }
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    ###### load observed data
    
    directory <- paste0(analysis.directory, "/Model_fits/", patient.id, "/Model_fit/")
    if(!file.exists(paste0(directory, "/", tissue, "/Model_fit.csv"))){next}
    fits <- read.csv(paste0(directory, "/", tissue, "/Model_fit.csv"))

    if(tissue == "CD34"){
      tissue.type = "CD34+"
    }else if(tissue == "MNC_minus_T"){
      tissue.type = "MNC_minus_T"
    }else if(tissue=="MNC"){
      tissue.type = "MNC"
    }else if(tissue=="PB_gran"){
      tissue.type="PB_gran"
    }

    ## get the drivers associated with this case
    driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(putative.drivers)[grep(tissue.type, colnames(putative.drivers))])]
    driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(patient.id, "_"), colnames(driver.information))])]
    driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", type="snvs", mutationcaller="mpileup", sample.col.mpileup = colnames(driver.information)[grep(paste0(patient.id, "_"), colnames(driver.information))])
    driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", type="snvs", mutationcaller="mpileup", sample.col.mpileup = colnames(driver.information)[grep(paste0(patient.id, "_"), colnames(driver.information))])
    driver.information$nvar <- driver.information$VAF*driver.information$Depth
    driver.information <- driver.information[driver.information$nvar>2,,drop=F]
    driver.information$lower <- driver.information$VAF-1.96*sqrt(driver.information$VAF*(1-driver.information$VAF)/driver.information$Depth) ## Wald approximation for 95% binomial CI
    driver.information$upper <- driver.information$VAF+1.96*sqrt(driver.information$VAF*(1-driver.information$VAF)/driver.information$Depth) ## Wald approximation for 95% binomial CI
    driver.information$lower[driver.information$lower<0] <- 0
    driver.information$upper[driver.information$upper>1] <- 1
    ## subset on drivers in "TET2", "DNMT3A", "ASXL1" - the others are putative drivers with unknown consequences
    driver.information <- driver.information[driver.information$GENE %in% c("DNMT3A", "TET2", "ASXL1"),]
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    #### Print parameter estimates
    
    pdf(paste0(directory, tissue, "/Parameter_estimates.pdf"), width=5, height=5)
    
    # the model converts the prior distribution for t_s and s into absolute values to match clones within the limits of min.prior.size and N. We here convert the parameter estimates accordingly.
    fits$par_t_s_absolute <- apply(fits, 1, function(x){
      min.t.s <- 0
      max.t.s <- age - log(0.01*10^as.numeric(x["par_N"]))/10^as.numeric(x["par_lambda_ss"])
      
      t.s <- min.t.s +as.numeric(x["par_t_s"])*(max.t.s - min.t.s)
    })
  
    fits$par_s_absolute <- apply(fits, 1, function(x){
      ts <- as.numeric(x["par_t_s_absolute"])

      min.s <- (10^as.numeric(x["par_lambda_ss"])*(age-ts) - log(10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*
                                                                                                                (age-ts))
      if(min.s < 0){
        min.s <- 0
      }
      max.s <- (10^as.numeric(x["par_lambda_ss"])*(age-ts) - log(0.01*10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*(age-ts))
      
      s <- min.s + as.numeric(x["par_s"])*(max.s-min.s)
      
    })
    
    fits$size_of_clone <- exp(10^fits$par_lambda_ss*(1-fits$par_s_absolute)*(age - fits$par_t_s_absolute))/10^fits$par_N

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
      facet_wrap("variable", scales="free") +
      scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none")+
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")
    
    print(p)
    
    ## plot again with log-scale 
    p <- ggplot(to.plot, aes(x=value, y= 0, fill = stat(quantile))) + 
      geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
      facet_wrap("variable", scales="free", nrow = 3, ncol = 4) +
      scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none")+
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability") + scale_x_log10()
    
    print(p)
    
    hdinterval <- as.data.frame(t(hdi(fits, credMass=0.8)))
    hdinterval$Parameter <- rownames(hdinterval) 
    hdinterval$Median <- apply(fits, 2, function(x){median(as.numeric(x))})
    hdinterval$Sample <- patient.id
    hdinterval$Tissue <- tissue
    
    ### 2D-correlations
    
    ## parameters to plot
    parameter.names <- c("par_N",  "par_lambda_ss", "par_mu",
                         "age_of_clone", "growth_per_year")
    
    ## specify the variables I want to plot
    meas_vars <- parameter.names
    
    ## a data frame of all combinations of its arguments
    
    controlTable <- data.frame(expand.grid(meas_vars, meas_vars, stringsAsFactors = F))
    
    ## rename the columns
    colnames(controlTable) <- c("x", "y")
    
    ## add the key column
    controlTable <- cbind(data.frame(par_key = paste(controlTable[[1]], controlTable[[2]]), stringsAsFactors = F), controlTable)
    
    ## create the new data frame
    to.plot <- rowrecs_to_blocks(fits, controlTable)
    
    ## re-arrange with facet_grid
    splt <- strsplit(to.plot$par_key, split=" ", fixed=TRUE)
    to.plot$xv <- vapply(splt, function(si) si[[1]], character(1))
    to.plot$yv <- vapply(splt, function(si) si[[2]], character(1))
    
    to.plot$xv <- factor(as.character(to.plot$xv), meas_vars)
    to.plot$yv <- factor(as.character(to.plot$yv), meas_vars)
    
    
    to.plot$xaxis <- F
    to.plot$yaxis <- F
    to.plot$xaxis[to.plot$yv == to.plot$xv[sqrt(length(unique(to.plot$par_key)))]] <- T
    to.plot$yaxis[to.plot$xv==to.plot$xv[1]] <- T
    to.plot$topm <- F
    to.plot$rightm <- F
    to.plot$topm[to.plot$yv == to.plot$xv[1]] <- T
    to.plot$rightm[to.plot$xv==to.plot$xv[sqrt(length(unique(to.plot$par_key)))]] <- T
    
    p <- list()
    
    ## introduce an artificial top row and right column
    
    for(i in 1:(sqrt(length(unique(to.plot$par_key)))+1)){
      p[[length(p)+1]] <- ggplot(data.frame()) + geom_point()+
        theme_bw() + theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"),
                           panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
      
    }
    
    for(i in unique(to.plot$par_key)){
      
      
      tmp <- to.plot[to.plot$par_key==i,]
      
      if(tmp$xv[1]==tmp$yv[1]){
        p[[length(p)+1]] <- ggplot(tmp, aes(x=x)) + 
          geom_histogram() + scale_x_continuous(name=tmp$xv[1]) + scale_y_continuous(name=tmp$yv[1])+
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
          theme(legend.position = "none")
        
      }else{

        p[[length(p)+1]] <- ggplot(tmp, aes(x=x, y=y)) + 
          geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) + 
          scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
          scale_x_continuous(name=tmp$xv[1]) + scale_y_continuous(name=tmp$yv[1])+
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
        
      }
      
      ## top-row and right column: adjust margins differently
      if(tmp$rightm[1] & tmp$topm[1]){
        p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
      }else if(tmp$rightm[1]){
        p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
      }else if(tmp$topm[1]){
        p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
      }else{
        p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
      }
      
      if(tmp$xaxis[1]==F){
        p[[length(p)]] <-  p[[length(p)]] +  theme(axis.title.x = element_blank(),
                                                   axis.text.x = element_blank())
      }
      
      if(tmp$yaxis[1]==F){
        p[[length(p)]] <-  p[[length(p)]] +  theme(axis.title.y = element_blank(),
                                                   axis.text.y = element_blank())
      }
      
      if(tmp$xv[1] %in% c("age_of_clone", "growth_per_year")){
        p[[length(p)]] <-  p[[length(p)]] +  scale_x_log10(name=tmp$xv[1])
      }
      
      if(tmp$yv[1] %in% c("age_of_clone", "growth_per_year")){
        p[[length(p)]] <-  p[[length(p)]] +  scale_y_log10(name=tmp$yv[1])
      }
      
      if(tmp$rightm[1]){
        p[[length(p)+1]] <- ggplot(data.frame()) + geom_point()+
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
          theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"),
                legend.position = "none")
        
      }
      
    }
    
    print(ggarrange(plotlist=p, nrow=6, ncol=6, align="hv"))
    
    
    ## for selection associated values, store only the parameters associated with a clone >= 0.05
    fits.selection <- fits[fits$size_of_clone >= 0.1, ]
    fits.neutral <- fits[fits$size_of_clone < 0.1, ]
    
    if(nrow(fits.selection)>0){
      model.support.selection[tissue, patient.id] <- nrow(fits.selection)/10
    }else{
      model.support.selection[tissue,patient.id] <- 0
    }
    
    if(nrow(fits.selection)==0){
      
      parameters <- rbind(parameters, hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                   "mutations_per_year"),])
      
      neutral.parameters <- rbind(neutral.parameters,hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                                  "mutations_per_year"),])
      
      p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free") 
      
      print(p)
      
      dev.off()
    }else{
      # compute hdi
      hdinterval.selection <- as.data.frame(t(hdi(fits.selection, credMass=0.8)))
      hdinterval.selection$Parameter <- rownames(hdinterval.selection) 
      hdinterval.selection$Median <- apply(fits.selection, 2, function(x){median(as.numeric(x))})
      hdinterval.selection$Sample <- patient.id
      hdinterval.selection$Tissue <- tissue
      
      p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free") 
      
      print(p)
      
      parameters <- rbind(parameters, hdinterval)
      
      selected.parameters <- rbind(selected.parameters, hdinterval.selection)
      
      if(nrow(fits.neutral)>0){
        hdi.neutral <- as.data.frame(t(hdi(fits.neutral, credMass = 0.8)))
        hdi.neutral$Parameter <- rownames(hdi.neutral) 
        hdi.neutral$Median <- apply(fits.neutral, 2, function(x){median(as.numeric(x))})
        hdi.neutral$Sample <- patient.id
        hdi.neutral$Tissue <- tissue
        
        neutral.parameters <- rbind(neutral.parameters, hdi.neutral[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                                      "mutations_per_year"),])
        
      }
      
      dev.off()
      
    }
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    ###### Plot fits for neutral and selected case
    
    source(paste0(custom.script.directory, "/Bayesian_fit.R"))
    
    if( !file.exists(paste0(analysis.directory, patient.id, "/Model_fit/", tissue, "/Sim_trajectories.RData") )){
      
      # simulate for 100 parameter sets
      sim<- matrix(0, nrow=100, ncol=length(mySumStatData$mutation.count[[1]]))
      
      for(j in 1:100){
        print(j)
        
        parms <- list(mu=fits$par_mu[j], N=fits$par_N[j], delta_exp = fits$par_delta_exp[j], lambda_ss=fits$par_lambda_ss[j],
                      offset=fits$par_offset[j], t_s=fits$par_t_s[j], s=fits$par_s[j])
        
        model <- myModel(parms)
        sim[j,] <- model$modelResult[[1]]
        
      }
      
      # compute 95% credible intervals
      max.pred<- apply(sim, 2, quantile, p=0.975)
      min.pred<- apply(sim, 2, quantile, p=0.025)
      
      data.vs.prediction <- data.frame(VAF=rep(vafs.of.interest), mean=mySumStatData$mutation.count[[1]],
                                       sd= mySumStatData$sampled.sd[[1]],
                                       min.model= min.pred, max.model=max.pred, 
                                       Age=mySumStatData$age/365)
      
      save(sim, data.vs.prediction, file=paste0(analysis.directory, patient.id, "/Model_fit/", tissue, "/Sim_trajectories.RData"))
      
    }else{
      load(paste0(analysis.directory, patient.id, "/Model_fit/", tissue, "/Sim_trajectories.RData"))
    }
    
    to.plot <- data.vs.prediction
    
    # plot the model fits
    grDevices::pdf(paste0(directory, tissue, "/Model_fit.pdf"), width=3, height=2.5)
    
    max.y <- max(to.plot$max.model)
    
    p <- ggplot(data=to.plot, aes(x=VAF, y=mean, ymin=mean-sd, ymax=mean+sd)) +
      geom_ribbon(data=to.plot, aes(x=VAF, y=mean.bm, ymin=min.model,ymax=max.model), alpha=1, fill=model.colors["selection"]) + ggtitle(tissue.type)+
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + scale_y_continuous(name="Cumulative # of variants") + 
      scale_x_continuous(limits=c(0, 0.6)) 
    
    if(nrow(driver.information)>0){
      p <- p + geom_segment(data = driver.information, aes(x = as.numeric(VAF), y = -250, xend = as.numeric(VAF), yend = 0),
                            arrow = arrow(length = unit(0.5, "cm")), col="firebrick", inherit.aes = F) +
        geom_text(data = driver.information, aes( x=as.numeric(VAF), y=-100), label=driver.information[,'GENE'],  color="firebrick",
                  size=5 ,  fontface="italic", inherit.aes = F)
    }
    
    print(p)
    
    xlimits <- c(0,20)
    
    p <- ggplot(data=to.plot, aes(x=1/VAF, y=mean, ymin=mean-sd, ymax=mean+sd)) +
      geom_ribbon(data=to.plot, aes(x=1/VAF, y=mean, ymin=min.model,ymax=max.model), alpha=1, 
                  fill=sample.color[tissue.type]) +
      ggtitle(paste(sample.info[patient.id,]$Paper_ID, tissue.type))+
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + 
      scale_y_continuous(name="Cumulative # of mutations") + theme(aspect.ratio = 1) +
      scale_x_continuous(limits=xlimits, breaks = c(5, 10, 20), labels = c("0.2", "0.1", "0.05"), name = "Variant allele frequency") + 
      coord_cartesian(ylim=c(0, max.y))
    
    if(nrow(driver.information[driver.information$upper>=1/max(xlimits),])>0){
      
      driver.information$y <- sapply(driver.information$VAF, function(x){
        to.plot[min(which(to.plot$VAF >= x)),]$max.model.BM
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
    plotlist.model.vs.data[[paste(patient.id, tissue.type)]] <- p
    
    dev.off()
    
    
    
  }

}

####################################################################################################################################################
## classify samples as selected or neutral according to CD34:

clearly.neutral.samples <- intersect(normal.samples, colnames(model.support.selection)[model.support.selection["CD34",]<15])
clearly.selected.samples <- setdiff(colnames(model.support.selection)[model.support.selection["CD34",]>=15],c(chip.samples.unknown.driver, normal.samples))
selection.no.driver <-  intersect(c(chip.samples.unknown.driver, normal.samples), colnames(model.support.selection)[model.support.selection["CD34",]>=15])
neutral.driver <- intersect(colnames(model.support.selection)[model.support.selection["CD34",]<15],c(chip.samples, chip.samples.unknown.driver))
####################################################################################################################################################
## Figures 4a/ 5a / 6a / S4, S6a/b: plot the model fits stratified by type

# clearly neutral:
pdf(paste0(analysis.directory, "/Figures/Figure_4a_S4.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[names(plotlist.model.vs.data) %in% 
                                            paste(clearly.neutral.samples, "CD34+")],
          nrow=6, ncol=6) 
dev.off()


# clearly selected:
pdf(paste0(analysis.directory, "/Figures/Figure_5a_S6_a.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[names(plotlist.model.vs.data) %in% 
                                            paste(clearly.selected.samples, "CD34+")],
          nrow=6, ncol=6) 

dev.off()

# selection for an unknown driver
pdf(paste0(analysis.directory, "/Figures/Figure_6_a.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[names(plotlist.model.vs.data) %in% 
                                            paste(selection.no.driver, "CD34+")],
          nrow=6, ncol=6) 
dev.off()

## for sample U5, plot a zoom in:

pdf(paste0(analysis.directory, "/Figures/Figure_6a_U5_zoom.pdf"), width=3.5, height=3.5)

plotlist.model.vs.data$`U6 CD34` + 
  scale_x_continuous( breaks = c(1/0.5, 1/0.4, 1/0.3, 1/0.2), labels = c("0.5", "0.4", "0.3", "0.2"), name="Variant allele frequency") + 
  scale_y_continuous( breaks=seq(0,15), labels=c("0", "", "", "", "", "5", "", "", "", "", "10", "", "", "", "", "15"), name = "Number of SSNVs") + 
  coord_cartesian(ylim=c(0, 15), xlim=c(1/0.5, 1/0.2))

dev.off()

## estimate time point of clone emergence: ~5-6 SSNVs

n.div <- 6/selected.parameters[selected.parameters$Sample=="U6" & selected.parameters$Tissue=="CD34" &
                                 selected.parameters$Parameter=="par_mu" ,c("lower", "Median", "upper")]

# driver but no evidence for selection

pdf(paste0(analysis.directory, "/Figures/Figure_S6_b_right.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[names(plotlist.model.vs.data) %in% 
                                            paste(neutral.driver, "CD34+")],
          nrow=6, ncol=6) 

dev.off()

####################################################################################################################################################
## Fig. 4b, 5b, 6b, S6b, S7: plot the posterior probability for the neutral and the selection model for each sample

to.plot <- melt(t(model.support.selection), value.name = "P_selection")
colnames(to.plot)[c(1,2)] <- c("Patient", "Sample")
to.plot$P_neutral <- 100 - to.plot$P_selection
to.plot <- melt(to.plot, value.name = "Posterior probability")
to.plot$Sample <- as.character(to.plot$Sample)
to.plot$Clone_size <- apply(to.plot, 1, function(x){
  res <- selected.parameters[selected.parameters$Sample==x[1] & selected.parameters$Parameter=="size_of_clone" &
                        selected.parameters$Tissue==x[2] & x[3]=="P_selection",]$Median
  if(length(res)==0){
    return(0)}else{
      return(res)
    }
})

to.plot$Sample <- factor(to.plot$Sample, levels = c("CD34", "MNC_minus_T", "MNC", "PB_gran"))
to.plot$variable <- factor(to.plot$variable, levels=c("P_neutral", "P_selection"))
to.plot$Clone_size[to.plot$Clone_size==0] <- NA
to.plot$ID <- sample.info[as.character(to.plot$Patient),]$Paper_ID
to.plot$CHIP_driver <- sample.info[as.character(to.plot$Patient),]$CHIP.mutation.associated.with.fit

## clearly neutral samples

pdf(paste0(analysis.directory, "/Figures/Figures_4b.pdf"), width=6, height=6)

to.plot.neutral <- to.plot[to.plot$Sample=="CD34" & !is.na(to.plot$`Posterior probability`) & to.plot$Patient %in% clearly.neutral.samples,]
to.plot.neutral$ID <- factor(to.plot.neutral$ID, levels=unique(to.plot.neutral$ID)[order(sapply(unique(to.plot.neutral$ID), function(x){
  sample.info[sample.info$Paper_ID==x,]$Age}))])

ggplot(to.plot.neutral,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34 neutral")+ theme( strip.background = element_blank() )+ geom_hline(yintercept = 15, linetype=2)

dev.off()

## clearly selected samples, associated with CH driver

pdf(paste0(analysis.directory, "/Figures/Figures_5b.pdf"), width=6, height=6)

to.plot.selected <- to.plot[to.plot$Sample=="CD34" & !is.na(to.plot$`Posterior probability`) & to.plot$Patient %in% clearly.selected.samples,]
ggplot(to.plot.selected,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, selection associated with CH driver")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

## neutral samples with evidence for selection

pdf(paste0(analysis.directory, "/Figures/Figures_6b_left.pdf"), width=6, height=6)

ggplot(to.plot[to.plot$Sample=="CD34" & !is.na(to.plot$`Posterior probability`) & to.plot$Patient %in% selection.no.driver,],
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34 selection for an unkonwn driver")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

dev.off()

## samples with driver but no evidence for selection

pdf(paste0(analysis.directory, "/Figures/Figures_S6b_left.pdf"), width=6, height=6)

ggplot(to.plot[to.plot$Sample=="CD34" & !is.na(to.plot$`Posterior probability`) & to.plot$Patient %in% neutral.driver,],
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34 no evidence for selection despite driver")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

dev.off()

## compare with other tissues
pdf(paste0(analysis.directory, "/Figures/Figures_7_posteriors.pdf"), width=6, height=6)

to.plot <- to.plot[!is.na(to.plot$`Posterior probability`),]
to.plot <- to.plot[to.plot$ID %in% names(table(to.plot$ID))[table(to.plot$ID)==2*4],] # select samples with all 4 tissues sequenced

ggplot(to.plot, aes(x=Sample, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col=NA) +
  facet_rep_wrap(~ID) +   geom_hline(yintercept = 15, linetype=2)+
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) + 
  theme( strip.background = element_blank() )

dev.off()


####################################################################################################################################################
## Figure 4c: plot raw data w/o fits for normal blood

to.plot <- data.frame(VAF= c(), MutationCount=c(), Sample=c(), Age=c())
vafs.of.interest <- seq(0.05, 1, 0.01)

for(i in clearly.neutral.samples){
  tmp <- snvs.cd34[[i]][snvs.cd34[[i]]$varCounts >=3 & snvs.cd34[[i]]$Depth >=10 &
                          snvs.cd34[[i]]$Depth <=300,]
  
  M <- sapply(vafs.of.interest, function(x){
    sum(tmp$VAF >=x)
  })  
  
  M <- M - min(M)
  age <- sample.info[sample.info$SAMPLE==i,]$Age
  
  to.plot <- rbind(to.plot, data.frame(VAF=vafs.of.interest,
                                       MutationCount=M,
                                       Sample=i,
                                       Age=age))
}

pdf(paste0(analysis.directory, "/Figures/Figure_4c.pdf"), width=3.5, height=3.5)

ggplot(to.plot, aes(x=1/VAF, y=MutationCount, col=Age, group=Sample)) + geom_point() + geom_line() +
  scale_color_gradientn(colors=hcl.colors(n=7, palette="Zissou 1"), limits=c(25, 80)) +
  scale_x_continuous(breaks = c(5, 10, 20), labels = c("0.2", "0.1", "0.05"), name="Variant allele frequency") +
  theme(aspect.ratio = 1)+ ggtitle("CD34")

dev.off()

####################################################################################################################################################
## Figure 4d, e, S7a: compare the physiological parameters in normal samples estimates across the cohort

pdf(paste0(analysis.directory, "/Figures/Figure_4d_e_S7a.pdf"), width=5, height=w)

## Stem cell number
## CD34+
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_N" & neutral.parameters$Tissue=="CD34" ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.neutral.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)

## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_N" &
                                                                        neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_N"&
                                                                        neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Tissue=="CD34")))

ggplot(to.plot, 
       aes(x=Patient.ID, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Stem cell number (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  + ggtitle("CD34") 


## Fig. S7a: compare to other tissues

to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_N"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.neutral.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$Tissue <- factor(to.plot$Tissue, levels=c("CD34", "MNC_minus_T", "MNC", "PB_gran"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed

ggplot(to.plot, 
       aes(x=Tissue, y=Median, ymin=lower, ymax=upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Stem cell number (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")


## Division rate

## CD34+
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss" & neutral.parameters$Tissue=="CD34",], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.neutral.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)

## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_lambda_ss" &
                                                                        neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_lambda_ss"&
                                                                        neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Tissue=="CD34")))

ggplot(to.plot[ to.plot$Tissue=="CD34",,drop=F], 
       aes(x=Patient.ID, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = 365*10^min, ymax = 365*10^max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  + ggtitle("CD34") 


## compare to other tissues (Fig. S7a)

to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.neutral.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$Tissue <- factor(to.plot$Tissue, levels=c("CD34", "MNC_minus_T", "MNC", "PB_gran"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed

ggplot(to.plot, 
       aes(x=Tissue, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")


## N times tau
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="N_tau" & neutral.parameters$Tissue=="CD34" ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.neutral.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)

## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="N_tau" &
                                                                        neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="N_tau"&
                                                                        neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Tissue=="CD34")))

ggplot(to.plot[ to.plot$Tissue=="CD34",,drop=F], 
       aes(x=Patient.ID, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="N x Tau") +
  theme(aspect.ratio = 1)+ expand_limits(y = 1)  + ggtitle("CD34") 


## Mutation rate

## CD34+ 

to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_mu" & neutral.parameters$Tissue=="CD34",], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.neutral.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)

## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_mu" &
                                                                        neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                              max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_mu"&
                                                                      neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                              x=c(0, sum( to.plot$Tissue=="CD34")))

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=Patient.ID, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  + ggtitle("CD34") 


## Fig. S7a: compare to other tissues

to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_mu"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.neutral.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$Tissue <- factor(to.plot$Tissue, levels=c("CD34", "MNC_minus_T", "MNC", "PB_gran"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed

ggplot(to.plot, 
       aes(x=Tissue, y=Median, ymin=lower, ymax=upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")


dev.off()


####################################################################################################################################################
## Fig. S6c: physiological parameters in a case with neutral dynamics but known CH driver 

pdf(paste0(analysis.directory, "/Figures/Figure_S6c.pdf"), width=5, height=3.5)

## Stem cell number

to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_N"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% neutral.driver,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))

## take medians of lower and upper quantiles of the neutral samples for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_N" &
                                                                        neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_N"&
                                                                        neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Tissue=="CD34")))

ggplot(to.plot, 
       aes(x=Patient.ID, y=Median, ymin=lower, ymax=upper)) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Stem cell number (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) 


## Mutation rate

to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_mu"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% neutral.driver,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))

## take medians of lower and upper quantiles of the neutral samples for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_mu" &
                                                                        neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                neutral.parameters$Parameter=="par_mu"&
                                                                neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Tissue=="CD34")))

ggplot(to.plot, 
       aes(x=Patient.ID, y=Median, ymin=lower, ymax=upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) 


## Division rate, compare between all sorts

to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% neutral.driver,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))

## take medians of lower and upper quantiles of the neutral samples for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_lambda_ss" &
                                                                        neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                        neutral.parameters$Parameter=="par_lambda_ss"&
                                                                        neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Tissue=="CD34")))

ggplot(to.plot, 
       aes(x=Patient.ID, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = 365*10^min, ymax = 365*10^max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) 

dev.off()

####################################################################################################################################################
## Figure 5c: compare the estimated clone size with the driver VAF

to.plot <- selected.parameters[selected.parameters$Parameter=="size_of_clone",]
to.plot$CHIP.mutation <- apply(to.plot, 1, function(x){
  as.character(sample.info[as.character(x["Sample"]),"CHIP.mutation.associated.with.fit"])
})
to.plot <- to.plot[to.plot$CHIP.mutation!="no driver" & !is.na(to.plot$CHIP.mutation),]

# get ML estimate of the driver VAF and compute 95% CI based on Binomial distr.
to.plot$Driver.mean <- apply(to.plot, 1, function(x){
  if(x["Tissue"]=="CD34"){
    type <- "CD34+"
  }else if(x["Tissue"]=="MNC_minus_T"){
    type <- "MNC_minus_T"
  }else if(x["Tissue"]=="MNC"){
    type <- "MNC"
  }else if(x["Tissue"]=="PB_gran"){
    type <- "PB_gran"
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                            colnames(putative.drivers)[grep(type, colnames(putative.drivers))])]
  if(type=="CD34+"){
    driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                                colnames(driver.information)[grep( "_CD38", colnames(driver.information), invert = T)])]
  }
  driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(x["Sample"], "_"), colnames(driver.information))])]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = ncol(driver.information))
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Sample"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  
  driver.information <- driver.information[which.max(driver.information$VAF),]
  vaf=as.numeric(driver.information["VAF"])
  if(is.na(vaf)){
    vaf <- 0
  }
  return(vaf)
})
to.plot$Driver.min <- apply(to.plot, 1, function(x){
  if(x["Tissue"]=="CD34"){
    type <- "CD34+"
  }else if(x["Tissue"]=="MNC_minus_T"){
    type <- "MNC_minus_T"
  }else if(x["Tissue"]=="MNC"){
    type <- "MNC"
  }else if(x["Tissue"]=="PB_gran"){
    type <- "PB_gran"
  }
  driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                            colnames(putative.drivers)[grep(type, colnames(putative.drivers))])]
  if(type=="CD34+"){
    driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                                colnames(driver.information)[grep( "_CD38", colnames(driver.information), invert = T)])]
  }
  driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(x["Sample"], "_"), colnames(driver.information))])]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = ncol(driver.information))
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = colnames(driver.information)[grep(paste0(x["Sample"], "_"), colnames(driver.information))])
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Sample"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  
  driver.information <- driver.information[which.max(driver.information$VAF),]
  vaf=as.numeric(driver.information["VAF"])
  DP=as.numeric(driver.information["Depth"])
  
  if(length(vaf)==0 || is.na(vaf)){
    vaf <- 0
    DP <- 0
  }
  lower=vaf-1.96*sqrt(vaf*(1-vaf)/DP) ## Wald approximation for 95% binomial CI
  return(lower)
})
to.plot$Driver.max <- apply(to.plot, 1, function(x){
  if(x["Tissue"]=="CD34"){
    type <- "CD34+"
  }else if(x["Tissue"]=="MNC_minus_T"){
    type <- "MNC_minus_T"
  }else if(x["Tissue"]=="MNC"){
    type <- "MNC"
  }else if(x["Tissue"]=="PB_gran"){
    type <- "PB_gran"
  }
  driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                            colnames(putative.drivers)[grep(type, colnames(putative.drivers))])]
  if(type=="CD34+"){
    driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                                colnames(driver.information)[grep( "_CD38", colnames(driver.information), invert = T)])]
  }
  driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(x["Sample"], "_"), colnames(driver.information))])]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = ncol(driver.information))
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = colnames(driver.information)[grep(paste0(x["Sample"], "_"), colnames(driver.information))])
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Sample"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  
  driver.information <- driver.information[which.max(driver.information$VAF),]
  vaf=as.numeric(driver.information["VAF"])
  DP=as.numeric(driver.information["Depth"])
  
  if(length(vaf)==0 || is.na(vaf)){
    vaf <- 0
    DP <- 0
  }
  lower=vaf+1.96*sqrt(vaf*(1-vaf)/DP) ## Wald approximation for 95% binomial CI
  return(lower)
})

pdf(paste0(analysis.directory, "Figures/Figure_5c"), width=4, height = 4)

ggplot(to.plot[to.plot$Driver.mean!=0 & to.plot$Sample %in% clearly.selected.samples,], aes(x=Median/2, xmin = lower/2, xmax = upper/2, color=CHIP.mutation,
                                                                                            y=Driver.mean, ymin = Driver.min, ymax=Driver.max)) + geom_point() +
  geom_errorbar() + geom_errorbarh() + scale_color_manual(values=CHIP.color) +
  geom_abline(slope = 1, intercept = 0, linetype=2) +
  facet_wrap(~Tissue, nrow=2) + scale_x_continuous("VAF population genetics model") +
  scale_y_continuous("VAF WGS") + theme(aspect.ratio = 1)

dev.off()


####################################################################################################################################################
## Fig. 5d, e/S7 b, c: Physiological parameters in selected cases

pdf(paste0(analysis.directory, "/Figures/Figure_5_d_e_S7_b_c.pdf"), width=5, height=3.5)

## Stem cell number

## selected clone
to.plot <- selected.parameters[selected.parameters$Parameter=="par_N" & selected.parameters$Sample %in% clearly.selected.samples,]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

## neutral as reference, take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.neutral <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                  neutral.parameters$Parameter=="par_N" &
                                                                  neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                              max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                  neutral.parameters$Parameter=="par_N"&
                                                                  neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                              x=c(0, sum( to.plot$Tissue=="CD34")))


## only HSCs

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median, ymin=lower, ymax=upper, col=CHIP.mutation)) + geom_pointrange() +
  geom_ribbon(data=to.plot.neutral, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="N (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~CHIP.mutation=="unknown driver", scales = "free")
   

## Fig. S7b, c, compare between all sorts

to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_N"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.selected.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$Tissue <- factor(to.plot$Tissue, levels=c("CD34", "MNC_minus_T", "MNC", "PB_gran"))
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed


ggplot(to.plot, 
       aes(x=Tissue, y=Median, ymin=lower, ymax=upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Stem cell number (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")


## N x tau

to.plot <- selected.parameters[selected.parameters$Parameter=="N_tau" & selected.parameters$Sample %in% clearly.selected.samples,]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown CHIP mutation", "healthy donor"))
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

## neutral as reference, take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.neutral <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="N_tau" &
                                                                      neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                              max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="N_tau"&
                                                                      neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                              x=c(0, sum( to.plot$Tissue=="CD34")))


## only HSCs 

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=paste( as.numeric(CHIP.mutation), ID), y=Median, ymin=lower, ymax=upper, col=CHIP.mutation)) + geom_pointrange() +
  geom_ribbon(data=to.plot.neutral, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="N x Tau") +
  theme(aspect.ratio = 1)+ expand_limits(y = 1) 

## Mutation rate

## selected clone
to.plot <- selected.parameters[selected.parameters$Parameter=="par_mu" & selected.parameters$Sample %in% clearly.selected.samples,]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

## neutral as reference, take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.neutral <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_mu" &
                                                                      neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                              max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_mu"&
                                                                      neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                              x=c(0, sum( to.plot$Tissue=="CD34")))


## only HSCs

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median, ymin=lower, ymax=upper, col=CHIP.mutation)) + geom_pointrange() +
  geom_ribbon(data=to.plot.neutral, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~CHIP.mutation=="unknown driver", scales = "free")

## Fig. S7b, c: compare between all sorts

to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_mu"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.selected.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed

ggplot(to.plot, 
       aes(x=Tissue, y=Median, ymin=lower, ymax=upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")


## Division rate

to.plot <- selected.parameters[selected.parameters$Parameter=="par_lambda_ss" & selected.parameters$Sample %in% clearly.selected.samples,]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

## neutral as reference, take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.neutral <- data.frame(min = rep(10^quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_lambda_ss" &
                                                                      neutral.parameters$Tissue=="CD34",]$lower, p = 0.025)*365, 2),
                              max = rep(10^quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_lambda_ss"&
                                                                      neutral.parameters$Tissue=="CD34",]$upper, p = 0.975)*365, 2),
                              x=c(0, sum( to.plot$Tissue=="CD34")))


## only HSCs

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x= ID, y=10^Median*365, ymin=10^lower*365, ymax=10^upper*365, col=CHIP.mutation)) + geom_pointrange() +
  geom_ribbon(data=to.plot.neutral, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~CHIP.mutation=="unknown driver", scales = "free")


## Fig. S7b, c: compare between all sorts

to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_lambda_ss"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.selected.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no driver"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed


ggplot(to.plot, 
       aes(x=Tissue, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")

dev.off()

####################################################################################################################################################
## Fig. S7_b_c: Compare the selective advantage by mutation

pdf(paste0(analysis.directory, "/Figures/Figure_S7_b_c_sel_adv.pdf"), width=5, height=3.5)

## compare between all sorts

to.plot <- melt(selected.parameters[selected.parameters$Parameter=="growth_per_year"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.selected.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$Tissue <- factor(to.plot$Tissue, levels=c("CD34", "MNC_minus_T", "MNC", "PB_gran"))
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed


ggplot(to.plot, 
       aes(x=Tissue, y=100*Median, ymin=100*lower, ymax=100*upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Growth per year (%)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")

dev.off()


####################################################################################################################################################
## Figure S7_b_c: Compare the age of the selected clone 

pdf(paste0(analysis.directory, "/Figures/Figure_S7_b_c_age.pdf"), width=5, height=3.5)

## compare between all sorts

to.plot <- melt(selected.parameters[selected.parameters$Parameter=="age_of_clone"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% clearly.selected.samples,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$Tissue <- factor(to.plot$Tissue, levels=c("CD34", "MNC_minus_T", "MNC", "PB_gran"))
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed

ggplot(to.plot, 
       aes(x=Tissue, y=Median, ymin=lower, ymax=upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of clone (years)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")

dev.off()

####################################################################################################################################################
## Fig. 5f-h: Compare physiological parameters between neutrally evolving cases and cases with selection

to.plot <- neutral.parameters[neutral.parameters$Sample %in% c(clearly.neutral.samples, neutral.driver),]
to.plot <- rbind(to.plot, selected.parameters[ selected.parameters$Sample %in% clearly.selected.samples,])
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID
to.plot$Age <- sample.info[to.plot$Sample,]$Age
to.plot$Type <- ifelse(to.plot$Sample %in% c(clearly.neutral.samples, neutral.driver), "Neutral", "Selection")

pdf(paste0(analysis.directory, "Figures/Figure_5_f_g_h.pdf"), width=5, height=2.5)

## Mutation rate
ggplot(to.plot[ to.plot$Tissue=="CD34" & to.plot$Parameter=="par_mu",], 
       aes(x=Type, y=Median, ymin=lower, ymax=upper, col = Type)) + geom_boxplot(outlier.shape = NA) + 
  geom_pointrange(position = position_dodge2(width=0.5)) + scale_color_manual(values=c(Neutral="darkgrey", Selection="red")) +
  scale_y_continuous(name ="Number of SSNVs per division") +
  expand_limits(y = 0)

## Stem cell number
ggplot(to.plot[ to.plot$Tissue=="CD34" & to.plot$Parameter=="par_N",], 
       aes(x=Type, y=Median, ymin=lower, ymax=upper, col = Type)) + geom_boxplot(outlier.shape = NA) + 
  geom_pointrange(position = position_dodge2(width=0.5)) + scale_color_manual(values=c(Neutral="darkgrey", Selection="red")) +
  scale_y_continuous(name ="Stem cell number") +
  expand_limits(y = 0)

## Division rate
ggplot(to.plot[ to.plot$Tissue=="CD34" & to.plot$Parameter=="par_lambda_ss",], 
       aes(x=Type, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col = Type)) + geom_boxplot(outlier.shape = NA) + 
  geom_pointrange(position = position_dodge2(width=0.5)) + scale_color_manual(values=c(Neutral="darkgrey", Selection="red")) +
  scale_y_continuous(name ="Division rate (1/y)") +
  expand_limits(y = 0)

ggplot(to.plot[ to.plot$Tissue=="CD34" & to.plot$Parameter=="par_lambda_ss",], 
       aes(x=Age, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col = Type)) +
  geom_point() + geom_errorbar() + scale_color_manual(values=c(Neutral="darkgrey", Selection="red")) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  expand_limits(x = 0, y = 0)

dev.off()

####################################################################################################################################################
## Figure 5i: compare the age of the selected clone among cases with known driver

to.plot <- selected.parameters[selected.parameters$Parameter=="age_of_clone"& selected.parameters$Sample %in% clearly.selected.samples,]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

pdf(paste0(analysis.directory, "Figures/Figure_5i.pdf"), width=3.5, height=2.5)

## only HSCs 

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median, ymin=lower, ymax=upper, col = CHIP.mutation)) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of clone (years)") +
  expand_limits(x = 0, y = 0)

dev.off()

####################################################################################################################################################
## Figure 5j: Compare the selective advantage between cases with and without known driver

to.plot <- selected.parameters[selected.parameters$Parameter=="growth_per_year" & selected.parameters$Sample %in% clearly.selected.samples,]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

pdf(paste0(analysis.directory, "Figures/Figure_5j.pdf"), width=5, height=2.5)

## only HSCs 

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median, ymin=lower, ymax=upper, col = CHIP.mutation)) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Growth per year (%)") +
  expand_limits(x = 0, y = 0)

dev.off()

####################################################################################################################################################
## Figure 6b, right panel: size of the selected clones

pdf(paste0(analysis.directory, "Figures/Figure_6b_right.pdf"), width=3.5, height=2.5)

to.plot <- selected.parameters[selected.parameters$Parameter=="size_of_clone" & selected.parameters$Sample %in% selection.no.driver,]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median/2, ymin=lower/2, ymax=upper/2)) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Size of clone (VAF)") +
  theme(aspect.ratio = 1) + expand_limits(y=0)

dev.off()

####################################################################################################################################################
## Fig. 6c/S7 d, e: Physiological parameters in selected cases

pdf(paste0(analysis.directory, "/Figures/Figure_6_c_S7_d_d.pdf"), width=3.5, height=3.5)

## Stem cell number

## selected clone
to.plot <- selected.parameters[selected.parameters$Parameter=="par_N" & selected.parameters$Sample %in% selection.no.driver,]
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

## neutral as reference, take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.neutral <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_N" &
                                                                      neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                              max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_N"&
                                                                      neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                              x=c(0, sum( to.plot$Tissue=="CD34")))


## only HSCs

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.neutral, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() +  
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="N (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) 


## Fig. S7d, e, compare between all sorts

to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_N"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selection.no.driver,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot$Tissue <- factor(to.plot$Tissue, levels=c("CD34", "MNC_minus_T", "MNC", "PB_gran"))
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed

ggplot(to.plot, 
       aes(x=Tissue, y=Median, ymin=lower, ymax=upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Stem cell number (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")


## Mutation rate

## selected clone
to.plot <- selected.parameters[selected.parameters$Parameter=="par_mu" & selected.parameters$Sample %in% selection.no.driver,]
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

## neutral as reference, take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.neutral <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_mu" &
                                                                      neutral.parameters$Tissue=="CD34",]$lower, p = 0.025), 2),
                              max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                      neutral.parameters$Parameter=="par_mu"&
                                                                      neutral.parameters$Tissue=="CD34",]$upper, p = 0.975), 2),
                              x=c(0, sum( to.plot$Tissue=="CD34")))

## only HSCs

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.neutral, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)

## Fig. S7d, e: compare between all sorts

to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_mu"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selection.no.driver,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed

ggplot(to.plot, 
       aes(x=Tissue, y=Median, ymin=lower, ymax=upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")


## Division rate

to.plot <- selected.parameters[selected.parameters$Parameter=="par_lambda_ss" & selected.parameters$Sample %in% selection.no.driver,]
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

## neutral as reference, take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.neutral <- data.frame(min = rep(10^quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                         neutral.parameters$Parameter=="par_lambda_ss" &
                                                                         neutral.parameters$Tissue=="CD34",]$lower, p = 0.025)*365, 2),
                              max = rep(10^quantile(neutral.parameters[neutral.parameters$Sample %in% clearly.neutral.samples & 
                                                                         neutral.parameters$Parameter=="par_lambda_ss"&
                                                                         neutral.parameters$Tissue=="CD34",]$upper, p = 0.975)*365, 2),
                              x=c(0, sum( to.plot$Tissue=="CD34")))


## only HSCs

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x= ID, y=10^Median*365, ymin=10^lower*365, ymax=10^upper*365)) + geom_pointrange() +
  geom_ribbon(data=to.plot.neutral, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) 


## Fig. S7d, e: compare between all sorts

to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_lambda_ss"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "Tissue"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selection.no.driver,]
to.plot$Patient.ID <- factor(sample.info[to.plot$Sample, ]$Paper_ID)
to.plot <- to.plot[to.plot$Patient.ID %in% names(table(to.plot$Patient.ID))[table(to.plot$Patient.ID)==4],] # only individuals with all 4 tissues analyzed

ggplot(to.plot, 
       aes(x=Tissue, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col=Tissue)) + geom_pointrange(fatten = 0.5) +
  geom_pointrange(fatten = 0.5) + scale_color_manual(values=sample.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + facet_wrap(~Patient.ID, nrow=1, scales = "free_y")

dev.off()

####################################################################################################################################################
## Figure 6d: age of the selected clone among cases with unknown driver

to.plot <- selected.parameters[selected.parameters$Parameter=="age_of_clone"& selected.parameters$Sample %in% selection.no.driver,]
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

pdf(paste0(analysis.directory, "Figures/Figure_6d.pdf"), width=5, height=2.5)

## only HSCs

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median, ymin=lower, ymax=upper)) +
  geom_pointrange() + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of clone") +
  expand_limits(x = 0, y = 0)

## compare to cases with known driver:

to.plot <- selected.parameters[selected.parameters$Parameter=="age_of_clone" & selected.parameters$Sample %in% c(selection.no.driver, clearly.selected.samples),]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown CHIP mutation"))

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=CHIP.mutation, y=Median, ymin=lower, ymax=upper, col = CHIP.mutation)) + geom_boxplot(outlier.shape = NA) + 
  geom_pointrange(position = position_dodge2(width=0.5)) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of clone") +
  expand_limits(x = 0, y = 0)


dev.off()

####################################################################################################################################################
## Figure 6e: Compare the selective advantage between cases with and without known driver

to.plot <- selected.parameters[selected.parameters$Parameter=="growth_per_year" & selected.parameters$Sample %in% selection.no.driver,]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

pdf(paste0(analysis.directory, "Figures/Figure_6e.pdf"), width=5, height=2.5)

## only HSCs 

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=ID, y=Median, ymin=lower, ymax=upper)) +
  geom_pointrange() +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Growth per year (%)") +
  expand_limits(x = 0, y = 0)

# compare to known drivers

to.plot <- selected.parameters[selected.parameters$Parameter=="growth_per_year" & selected.parameters$Sample %in% c(selection.no.driver, clearly.selected.samples),]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown CHIP mutation"))

ggplot(to.plot[ to.plot$Tissue=="CD34",], 
       aes(x=CHIP.mutation, y=Median, ymin=lower, ymax=upper, col = CHIP.mutation)) + geom_boxplot(outlier.shape = NA) + 
  geom_pointrange(position = position_dodge2(width=0.5)) +
  scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Growth per year (%)") +
  expand_limits(x = 0, y = 0)

dev.off()

####################################################################################################################################################
## Fig. 7: plot the incidence of CH driver acquisition

to.plot <- selected.parameters[selected.parameters$Parameter=="par_t_s_absolute"& 
                                 selected.parameters$Sample %in% c(selection.no.driver, clearly.selected.samples) &
                                 selected.parameters$Tissue=="CD34",]
to.plot$CHIP.mutation <- sample.info[to.plot$Sample,]$CHIP.mutation.associated.with.fit
to.plot$ID <- sample.info[to.plot$Sample,]$Paper_ID

## driver info:
to.plot.driver <- to.plot
to.plot.driver <- to.plot.driver[order(to.plot.driver$Median),]
to.plot.driver$y <- (1:nrow(to.plot.driver))/nrow(to.plot.driver)

## cohort summary
to.plot <- data.frame(Age = seq(0, 100), 
                      Incidence = sapply(seq(0,100), function(a){sum(to.plot$Median/365<=a)})/nrow(to.plot),
                      Lower = sapply(seq(0,100), function(a){sum(to.plot$upper/365<=a)})/nrow(to.plot),
                      Upper = sapply(seq(0,100), function(a){sum(to.plot$lower/365<=a)})/nrow(to.plot))


pdf(paste0(analysis.directory, "Figures/Figure_7.pdf"), width=5, height=2.5)

ggplot(to.plot, aes(x=Age, y=Incidence, ymin=Lower, ymax=Upper)) + geom_ribbon(fill="grey", alpha=0.5) + geom_line() +
  geom_point(data=to.plot.driver, aes(x=Median/365, y=y, col=CHIP.mutation), inherit.aes = F) + 
  geom_errorbarh(data=to.plot.driver, aes(xmin=lower/365, xmax=upper/365, col=CHIP.mutation, y=y), inherit.aes = F, height=0) +
  scale_color_manual(values=CHIP.color)

dev.off()
