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
load("RData/WGS_heme/SNVs.RData")

putative.drivers <- read.xlsx("MetaData/Supplementary Tables.xlsx", startRow = 8, sheet = 7)
sample.info <- read.xlsx("MetaData/Supplementary Tables.xlsx", sheet = 2, startRow = 6)
rownames(sample.info) <- sample.info$Paper_ID
patient.ids <- sample.info$Paper_ID


############################################################################################################################################
####### Set parameters as used for model fits
use.sensitivity <- F
sample.color["CD34+"] <- sample.color["CD34"]
sample.color["CD34+_deep"] <- sample.color["CD34"]
seq.type <- "bulk"

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Plot per patient and compare parameter across cohort

# collect highest density intervals of the parameters:
# .. overall ..
parameters <- data.frame()
# .. neutral fits only ..
neutral.parameters.1cm <- data.frame()
# .. selection fits only
selected.parameters.1cm <- data.frame()

# collect model support for selection (%)
model.support.selection <- matrix(NA, nrow=5, ncol=length(patient.ids), dimnames = list(c("MNC", "CD34", "CD34_deep", "MNC-T", "PB"),
                                                                                      patient.ids))

plotlist.model.vs.data <- list()

for(patient.id in patient.ids){
  print(patient.id)
  age <- sample.info[patient.id, ]$Age*365

  cell.sorts <- c("CD34", "MNC", "MNC-T", "PB", "CD34_deep")

  for(sort in cell.sorts){

    min.vaf <- 0.05
    min.clone.size <- 0.05
    min.prior.size <- 0.01
    
    # specify coverage and snvs of this sample
    if(sort == "CD34"){
      depth <- as.numeric(sample.info[patient.id,]$`Coverage.WGS.CD34+.1`)
      if(is.na(depth) | depth==0){next}
      snvs <- list(snvs.cd34[[patient.id]])
    }else if(sort == "MNC"){
      depth <- as.numeric(sample.info[patient.id,]$Coverage.WGS.BM.MNC)
      if(is.na(depth) | depth==0){next}
      snvs <- list(snvs.mnc[[patient.id]])
    }else if (sort =="MNC-T"){
      depth <- as.numeric(sample.info[patient.id,]$`Coverage.WGS.BM.MNC–T`)
      if(is.na(depth) | depth==0){next}
      snvs <- list(snvs.mnc_minus_t[[patient.id]])
    }else if (sort == "CD34_deep"){
      depth <- as.numeric(sample.info[patient.id,]$`Coverage.WGS.CD34+.1` + sample.info[patient.id,]$`Coverage.WGS.CD34+.2` )
      if(is.na(depth) | depth==0){next}
      snvs <- list(snvs.cd34.deep[[patient.id]])
      min.vaf <- 0.02
      min.clone.size <- 0.01
      min.prior.size <- 0.001
    }else{
      depth <- as.numeric(sample.info[patient.id,]$Coverage.WGS.PB.granulocytes)
      if(is.na(depth) | depth==0){next}
      snvs <- list(snvs.pb_gran[[patient.id]])
    }

    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
    ###### load observed data

    directory <- paste0(analysis.directory, "/Model_fits/WGS_heme/", paste(patient.id, sort, sep="_"))
    if(!file.exists(paste0(directory, "/Model_fit.csv"))){next}
    fits <- read.csv(paste0(directory, "/Model_fit.csv"))

    if(sort == "CD34"){
      sort.type = "CD34+"
    }else if(sort == "MNC-T"){
      sort.type = "MNC(-T)"
    }else if(sort=="MNC"){
      sort.type = "MNC"
    }else if(sort=="PB"){
      sort.type="PB_gran"
    }else if(sort=="CD34_deep"){
      sort.type="CD34+_deep"
    }

    ## get the drivers associated with this case
    driver.information <- putative.drivers[,c("FORMAT","CHROM", "POS", "REF", "ALT", "GENE", "AAchange",paste0(patient.id, "_", sort.type))]
    driver.information <- driver.information[!is.na(driver.information[,paste0(patient.id, "_", sort.type)]),]
    if(any(grepl(paste0(patient.id, "_"), colnames(driver.information)))){
      driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", type="snvs", mutationcaller="mpileup", sample.col.mpileup = paste0(patient.id, "_", sort.type))
      driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", type="snvs", mutationcaller="mpileup", sample.col.mpileup = paste0(patient.id, "_", sort.type))
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
    driver.information <- driver.information[driver.information$VAF!=1 & driver.information$GENE %in% c("DNMT3A", "TET2", "ASXL1") |
                                               (driver.information$GENE == "KMT2D" & patient.id == "T2"),]

    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
    #### Print parameter estimates

    pdf(paste0(directory, "/Parameter_estimates.pdf"), width=5, height=5)

    # the model converts the prior distribution for t_s and s into absolute values to match clones within the limits of min.prior.size and N. We here convert the parameter estimates accordingly.
    fits$par_t_s_absolute <- apply(fits, 1, function(x){
      min.t.s <- 0
      max.t.s <- age - log(min.prior.size*10^as.numeric(x["par_N"]))/10^as.numeric(x["par_lambda_ss"])
      
      t.s <- min.t.s +as.numeric(x["par_t_s"])*(max.t.s - min.t.s)
    })
    
    
    fits$par_s_absolute <- apply(fits, 1, function(x){
      ts <- as.numeric(x["par_t_s_absolute"])
      
      min.s <- (10^as.numeric(x["par_lambda_ss"])*(age-ts) - log(10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*
                                                                                                (age-ts))
      if(min.s < 0){
        min.s <- 0
      }
      max.s <- (10^as.numeric(x["par_lambda_ss"])*(age-ts) - log(min.prior.size*10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*(age-ts))
      
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
    hdinterval$Paper_ID <- patient.id
    hdinterval$Sort <- sort
    hdinterval$Depth <- depth
    
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
    fits.selection <- fits[fits$size_of_clone >= 2*min.vaf, ]
    fits.neutral <- fits[fits$size_of_clone < 2*min.vaf, ]

    if(nrow(fits.selection)>0){
      model.support.selection[sort, patient.id] <- nrow(fits.selection)/10
    }else{
      model.support.selection[sort,patient.id] <- 0
    }

    if(nrow(fits.selection)==0){

      parameters <- rbind(parameters, hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                   "mutations_per_year"),])

      neutral.parameters.1cm <- rbind(neutral.parameters.1cm,hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
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
      hdinterval.selection$Paper_ID <- patient.id
      hdinterval.selection$Sort <- sort
      hdinterval.selection$Depth <- depth
      
      p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free")

      print(p)

      parameters <- rbind(parameters, hdinterval)

      selected.parameters.1cm <- rbind(selected.parameters.1cm, hdinterval.selection)

      if(nrow(fits.neutral)>0){
        hdi.neutral <- as.data.frame(t(hdi(fits.neutral, credMass = 0.8)))
        hdi.neutral$Parameter <- rownames(hdi.neutral)
        hdi.neutral$Median <- apply(fits.neutral, 2, function(x){median(as.numeric(x))})
        hdi.neutral$Paper_ID <- patient.id
        hdi.neutral$Sort <- sort
        hdi.neutral$Depth <- depth
        
        neutral.parameters.1cm <- rbind(neutral.parameters.1cm, hdi.neutral[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                                      "mutations_per_year"),])

      }

      dev.off()

    }

    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
    ###### Plot fits for neutral and selected case

    source(paste0(custom.script.directory, "/Parameter_estimation/Bayesian_fit.R"))

    if( !file.exists(paste0(directory, "/Sim_trajectories.RData"))){

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

      data.vs.prediction <- data.frame(VAF=rep(vafs.of.interest), mean.data=mySumStatData$mutation.count[[1]],
                                       sd.data= mySumStatData$sampled.sd[[1]],
                                       min.model= min.pred, max.model=max.pred,
                                       Age=mySumStatData$age/365)

      save(sim, data.vs.prediction, file=paste0(directory, "/Sim_trajectories.RData"))

    }else{
      load(paste0(directory, "/Sim_trajectories.RData"))
    }

    to.plot <- data.vs.prediction

    # plot the model fits
    grDevices::pdf(paste0(directory, "/Model_fit.pdf"), width=3, height=2.5)

    max.y <- max(to.plot$max.model)

    p <- ggplot(data=to.plot, aes(x=VAF, y=mean.data, ymin=mean.data-sd.data, ymax=mean.data+sd.data)) +
      geom_ribbon(data=to.plot, aes(x=VAF, y=mean.data, ymin=min.model,ymax=max.model), alpha=1, fill=model.colors["selection"]) + ggtitle(sort.type)+
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + scale_y_continuous(name="Cumulative # of variants") +
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
                  fill=sample.color[sort.type]) +
      ggtitle(paste(sample.info[patient.id,]$Paper_ID, sort.type))+
      geom_pointrange(lwd=0.25, shape=1, fatten=1) +
      scale_y_continuous(name="Cumulative # of mutations") + theme(aspect.ratio = 1) +
      scale_x_continuous(limits=xlimits, breaks = c(5, 10, 20, 50), labels = c("0.2", "0.1", "0.05", "0.02"), name = "Variant allele frequency") +
      coord_cartesian(ylim=c(0, max.y))

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

    if(model.support.selection[sort, patient.id] > 15){

      p <- p + geom_ribbon(data = data.frame(x = unlist(hdinterval.selection[hdinterval.selection$Parameter=="size_of_clone",c("lower", "upper")]),
                                             ymin = c(0,0),
                                             ymax = c(max.y, max.y)), aes(x=2/x, ymin=ymin, ymax = ymax), inherit.aes = F, fill="grey", alpha = 0.5)
    }


    print(p)
    plotlist.model.vs.data[[paste(patient.id, sort.type)]] <- p

    dev.off()



  }

}

####################################################################################################################################################
## classify samples as selected or neutral according to CD34:

## CD34+ samples that are classified as neutral/selected at 90x WGS at a selection threshold of 15% posterior probability
neutral.samples.90 <- colnames(model.support.selection)[model.support.selection["CD34",]<15]
selected.samples.90 <- colnames(model.support.selection)[model.support.selection["CD34",]>=15]

## CD34+ samples that are classified as neutral/selected at 270x WGS at a selection threshold of 15%
all.270x.samples <- colnames(model.support.selection)[!is.na(model.support.selection["CD34_deep",])] ## we don't have 270x data for all samples
neutral.samples.270 <- colnames(model.support.selection)[!is.na(model.support.selection["CD34_deep",]) & model.support.selection["CD34_deep",]<15]
selected.samples.270 <- colnames(model.support.selection)[!is.na(model.support.selection["CD34_deep",]) & model.support.selection["CD34_deep",]>=15]

