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

putative.drivers <- read.xlsx("MetaData/Supplementary Tables.xlsx", startRow = 8, sheet = 7)
sample.info <- read.xlsx("MetaData/Supplementary Tables.xlsx", sheet = 2, startRow = 6)
rownames(sample.info) <- sample.info$Paper_ID
patient.ids <- sample.info$Paper_ID


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
# .. neutral fits only ..
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

  cell.sorts <- c("CD34", "MNC", "MNC_minus_T", "PB_gran", "CD34_deep")

  for(sort in cell.sorts){

    min.vaf <- 0.05
    min.clone.size <- 0.05
    min.prior.size <- 0.01
    
    # specify coverage and snvs of this sample
    if(sort == "CD34"){
      depth <- sample.info[patient.id,]$`Coverage.WGS.CD34+.1`
      if(depth==0){next}
      snvs <- list(snvs.cd34[[patient.id]])
    }else if(sort == "MNC"){
      depth <- sample.info[patient.id,]$Coverage.WGS.BM.MNC
      if(depth==0){next}
      snvs <- list(snvs.mnc[[patient.id]])
    }else if (sort =="MNC_minus_T"){
      depth <- sample.info[patient.id,]$`Coverage.WGS.BM.MNC–T`
      if(depth==0){next}
      snvs <- list(snvs.mnc_minus_t[[patient.id]])
    }else if (sort == "CD34_deep"){
      depth <- sample.info[patient.id,]$`Coverage.WGS.CD34+.1` + sample.info[patient.id,]$`Coverage.WGS.CD34+.2` 
      if(depth==0){next}
      snvs <- list(snvs.cd34.deep[[patient.id]])
      min.vaf <- 0.02
      min.clone.size <- 0.01
      min.prior.size <- 0.001
    }else{
      depth <- sample.info[patient.id,]$Coverage.WGS.PB.granulocytes
      if(depth==0){next}
      snvs <- list(snvs.pb_gran[[patient.id]])
    }

    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
    ###### load observed data

    directory <- paste0(analysis.directory, "/Model_fits/WGS/", paste(patient.id, sort, sep="_"))
    if(!file.exists(paste0(directory, "/Model_fit.csv"))){next}
    fits <- read.csv(paste0(directory, "/Model_fit.csv"))

    if(sort == "CD34"){
      sort.type = "CD34+"
    }else if(sort == "MNC_minus_T"){
      sort.type = "MNC_minus_T"
    }else if(sort=="MNC"){
      sort.type = "MNC"
    }else if(sort=="PB_gran"){
      sort.type="PB_gran"
    }else if(sort=="CD34_deep"){
      sort.type="CD34+_deep"
    }

    ## get the drivers associated with this case
    driver.information <- putative.drivers[putative.drivers$mutation_in_any_control<=5,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                                                                         colnames(putative.drivers)[grep(sort.type, colnames(putative.drivers))])]
    if(sort.type=="CD34+"){
      driver.information <- driver.information[,colnames(driver.information)[grep( "_CD38", colnames(driver.information), invert = T)]]
      driver.information <- driver.information[,colnames(driver.information)[grep( "deep", colnames(driver.information), invert = T)]]
    }else if(sort.type == "CD34+_deep"){
      driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                                colnames(putative.drivers)[grepl("CD34", colnames(putative.drivers)) &
                                                                             grepl("deep", colnames(putative.drivers))])]
      
    }
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
    }else{
      driver.information <- data.frame()
    }

    ## subset on drivers in "TET2", "DNMT3A", "ASXL1" and, in T2, KMT2D - the others are putative drivers with unknown consequences
    driver.information <- driver.information[driver.information$GENE %in% c("DNMT3A", "TET2", "ASXL1") |
                                               (driver.information$GENE == "KMT2D" & patient.id == "T2"),]

    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
    #### Print parameter estimates

    pdf(paste0(directory, "/Parameter_estimates.pdf"), width=5, height=5)

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
      hdinterval.selection$Paper_ID <- patient.id
      hdinterval.selection$Sort <- sort
      hdinterval.selection$Depth <- depth
      
      p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free")

      print(p)

      parameters <- rbind(parameters, hdinterval)

      selected.parameters <- rbind(selected.parameters, hdinterval.selection)

      if(nrow(fits.neutral)>0){
        hdi.neutral <- as.data.frame(t(hdi(fits.neutral, credMass = 0.8)))
        hdi.neutral$Parameter <- rownames(hdi.neutral)
        hdi.neutral$Median <- apply(fits.neutral, 2, function(x){median(as.numeric(x))})
        hdi.neutral$Paper_ID <- patient.id
        hdi.neutral$Sort <- sort
        hdi.neutral$Depth <- depth
        
        neutral.parameters <- rbind(neutral.parameters, hdi.neutral[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
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

      data.vs.prediction <- data.frame(VAF=rep(vafs.of.interest), mean=mySumStatData$mutation.count[[1]],
                                       sd= mySumStatData$sampled.sd[[1]],
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

    p <- ggplot(data=to.plot, aes(x=VAF, y=mean, ymin=mean-sd, ymax=mean+sd)) +
      geom_ribbon(data=to.plot, aes(x=VAF, y=mean.bm, ymin=min.model,ymax=max.model), alpha=1, fill=model.colors["selection"]) + ggtitle(sort.type)+
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
                  fill=sample.color[sort.type]) +
      ggtitle(paste(sample.info[patient.id,]$Paper_ID, sort.type))+
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
    plotlist.model.vs.data[[paste(patient.id, sort.type)]] <- p

    dev.off()



  }

}

####################################################################################################################################################
## classify samples as selected or neutral according to CD34:

## 90x
clearly.neutral.samples.90 <- intersect(normal.samples, colnames(model.support.selection)[model.support.selection["CD34",]<15])
clearly.selected.samples.90 <- setdiff(colnames(model.support.selection)[model.support.selection["CD34",]>=15],c(chip.samples.unknown.driver, normal.samples))
selection.no.driver.90 <-  intersect(c(chip.samples.unknown.driver, normal.samples), colnames(model.support.selection)[model.support.selection["CD34",]>=15])

## 270x
all.270x.samples <- colnames(model.support.selection)[!is.na(model.support.selection["CD34_deep",])]

clearly.neutral.samples.270 <- intersect(normal.samples, colnames(model.support.selection)[model.support.selection["CD34_deep",]<15])
clearly.selected.samples.270 <- setdiff(colnames(model.support.selection)[model.support.selection["CD34_deep",]>=15],c(chip.samples.unknown.driver, normal.samples))
selection.no.driver.270 <-  intersect(c(chip.samples.unknown.driver, normal.samples), colnames(model.support.selection)[model.support.selection["CD34_deep",]>=15])


####################################################################################################################################################
## Figures 4c/ 5a / 6c / S6, S7a/b: plot the model fits stratified by type

# no known CH driver, 90x:
pdf(paste0(analysis.directory, "/Figures/Figure_4c_S6a.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[paste(names(plotlist.model.vs.data)[
  grepl("N", names(plotlist.model.vs.data)) | grepl("U", names(plotlist.model.vs.data))], "CD34+")],
          nrow=6, ncol=6)
dev.off()

# no known CH driver, 2700x:
pdf(paste0(analysis.directory, "/Figures/Figure_6b.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[paste(names(plotlist.model.vs.data)[
  grepl("N", names(plotlist.model.vs.data)) | grepl("U", names(plotlist.model.vs.data))], "CD34+_deep")],
  nrow=6, ncol=6)
dev.off()


# known CH driver, 90x:
pdf(paste0(analysis.directory, "/Figures/Figure_5a_S7a.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[paste(names(plotlist.model.vs.data)[
  !grepl("N", names(plotlist.model.vs.data)) & !grepl("U", names(plotlist.model.vs.data))], "CD34+")],
  nrow=6, ncol=6)

dev.off()

# known CH driver, 270x
pdf(paste0(analysis.directory, "/Figures/Figure_6c.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[paste(names(plotlist.model.vs.data)[
  !grepl("N", names(plotlist.model.vs.data)) & !grepl("U", names(plotlist.model.vs.data))], "CD34+")],
  nrow=6, ncol=6)

dev.off()

## for sample U8, plot a zoom in:

pdf(paste0(analysis.directory, "/Figures/Figure_6c_U8_zoom.pdf"), width=3.5, height=3.5)

plotlist.model.vs.data$`U6 CD34` +
  scale_x_continuous( breaks = c(1/0.5, 1/0.4, 1/0.3, 1/0.2), labels = c("0.5", "0.4", "0.3", "0.2"), name="Variant allele frequency") +
  scale_y_continuous( breaks=seq(0,15), labels=c("0", "", "", "", "", "5", "", "", "", "", "10", "", "", "", "", "15"), name = "Number of SSNVs") +
  coord_cartesian(ylim=c(0, 15), xlim=c(1/0.5, 1/0.2))

dev.off()

## estimate time point of clone emergence: ~5-6 SSNVs

n.div <- 6/selected.parameters[selected.parameters$Paper_ID=="U6" & selected.parameters$Sort=="CD34" &
                                 selected.parameters$Parameter=="par_mu" ,c("lower", "Median", "upper")]


####################################################################################################################################################
## Fig. 4b,d 5b, 6a, Supplementary Figures 3: plot the posterior probability for the neutral and the selection model for each sample

to.plot <- melt(t(model.support.selection), value.name = "P_selection")
colnames(to.plot)[c(1,2)] <- c("Patient", "Sort")
to.plot$P_neutral <- 100 - to.plot$P_selection
to.plot <- melt(to.plot, value.name = "Posterior probability")
to.plot$Sort <- as.character(to.plot$Sort)
to.plot$Clone_size <- apply(to.plot, 1, function(x){
  res <- selected.parameters[selected.parameters$Sort==x[1] & selected.parameters$Parameter=="size_of_clone" &
                        selected.parameters$Sort==x[2] & x[3]=="P_selection",]$Median
  if(length(res)==0){
    return(0)}else{
      return(res)
    }
})

to.plot$Sort <- factor(to.plot$Sort, levels = c("CD34", "CD34_deep", "MNC_minus_T", "MNC", "PB_gran"))
to.plot$variable <- factor(to.plot$variable, levels=c("P_neutral", "P_selection"))
to.plot$Clone_size[to.plot$Clone_size==0] <- NA
to.plot$ID <- sample.info[as.character(to.plot$Patient),]$Paper_ID
to.plot$CHIP_driver <- sample.info[as.character(to.plot$Patient),]$CHIP.mutation.associated.with.fit

## samples w/o CH driver

pdf(paste0(analysis.directory, "/Figures/Figures_4b_d.pdf"), width=6, height=6)

## 90x
to.plot.noCH <- to.plot[to.plot$Sort=="CD34" & !is.na(to.plot$`Posterior probability`) & 
                          to.plot$ID %in% sample.info[sample.info$CH.driver.found == "no",]$Paper_ID,]
to.plot.noCH$ID <- factor(to.plot.noCH$ID, levels=unique(to.plot.noCH$ID)[order(sapply(unique(to.plot.noCH$ID), function(x){
  sample.info[sample.info$Paper_ID==x,]$Age}))])

ggplot(to.plot.noCH,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34 no CH driver, 90x")+ theme( strip.background = element_blank() )+ geom_hline(yintercept = 15, linetype=2)

## 270x
to.plot.noCH <- to.plot[to.plot$Sort=="CD34_deep" & !is.na(to.plot$`Posterior probability`)  & 
                          to.plot$ID %in% sample.info[sample.info$CH.driver.found == "no",]$Paper_ID,]
to.plot.noCH$ID <- factor(to.plot.noCH$ID, levels=unique(to.plot.noCH$ID)[order(sapply(unique(to.plot.noCH$ID), function(x){
  sample.info[sample.info$Paper_ID==x,]$Age}))])

ggplot(to.plot.noCH,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34 no CH driver, 270x")+ theme( strip.background = element_blank() )+ geom_hline(yintercept = 15, linetype=2)

dev.off()

## samples with known CH driver

pdf(paste0(analysis.directory, "/Figures/Figures_5b.pdf"), width=6, height=6)

# 90x
to.plot.selected <- to.plot[to.plot$Sort=="CD34" & !is.na(to.plot$`Posterior probability`)  & 
                              to.plot$ID %in% sample.info[!sample.info$CH.driver.found == "yes",]$Paper_ID,]
ggplot(to.plot.selected,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, known CH driver, 90x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

# 270x
to.plot.selected <- to.plot[to.plot$Sort=="CD34_deep" & !is.na(to.plot$`Posterior probability`)  & 
                              to.plot$ID %in% sample.info[!sample.info$CH.driver.found == "yes",]$Paper_ID,]

ggplot(to.plot.selected,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, known CH driver, 270x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

dev.off()

## samples with unknown CH driver

pdf(paste0(analysis.directory, "/Figures/Figures_6bd.pdf"), width=6, height=6)

# 90x
to.plot.selected <- to.plot[to.plot$Sort=="CD34" & !is.na(to.plot$`Posterior probability`)  & 
                              grepl("U", to.plot$ID),]
ggplot(to.plot.selected,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, unknown CH driver, 90x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

# 270x
to.plot.selected <- to.plot[to.plot$Sort=="CD34_deep" & !is.na(to.plot$`Posterior probability`)  & 
                              grepl("U", to.plot$ID),]

ggplot(to.plot.selected,
       aes(x=ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, unknown CH driver, 270x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

dev.off()


## compare with other sorts
pdf(paste0(analysis.directory, "/Figures/Supplementary_Figure_3.pdf"), width=6, height=6)

to.plot <- to.plot[!is.na(to.plot$`Posterior probability`),]
to.plot <- to.plot[to.plot$ID %in% names(table(to.plot$ID))[table(to.plot$ID)==2*4],] # select samples with all 4 sorts sequenced

ggplot(to.plot, aes(x=Sort, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col=NA) +
  facet_rep_wrap(~ID) +   geom_hline(yintercept = 15, linetype=2)+
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  theme( strip.background = element_blank() )

dev.off()


####################################################################################################################################################
## Figure 4c: plot raw data w/o fits for normal blood

to.plot <- data.frame(VAF= c(), MutationCount=c(), Paper_ID=c(), Age=c(), Depth=c())

## 90x
vafs.of.interest <- seq(0.05, 1, 0.01)

for(i in sample.info[sample.info$CH.driver.found=="no",]$Paper_ID){
  tmp <- snvs.cd34[[i]][snvs.cd34[[i]]$varCounts >=3 & snvs.cd34[[i]]$Depth >=10 &
                          snvs.cd34[[i]]$Depth <=300,]

  M <- sapply(vafs.of.interest, function(x){
    sum(tmp$VAF >=x)
  })

  M <- M - min(M)
  age <- sample.info[sample.info$Paper_ID==i,]$Age

  to.plot <- rbind(to.plot, data.frame(VAF=vafs.of.interest,
                                       MutationCount=M,
                                       Paper_ID=i,
                                       Age=age,
                                       Depth=90))
}

## 270x
vafs.of.interest <- seq(0.02, 1, 0.01)

for(i in sample.info$Paper_ID[sample.info$CH.driver.found=="no" &
                            sample.info$`Coverage.WGS.CD34+.2`>0]){
  tmp <- snvs.cd34.deep[[i]]
  
  M <- sapply(vafs.of.interest, function(x){
    sum(tmp$VAF >=x)
  })
  
  M <- M - min(M)
  age <- sample.info[sample.info$Paper_ID==i,]$Age
  
  to.plot <- rbind(to.plot, data.frame(VAF=vafs.of.interest,
                                       MutationCount=M,
                                       Paper_ID=i,
                                       Age=age,
                                       Depth = 270))
}



pdf(paste0(analysis.directory, "/Figures/Figure_4a.pdf"), width=3.5, height=3.5)

ggplot(to.plot[to.plot$Depth==90,], aes(x=1/VAF, y=MutationCount, col=Age, group=Paper_ID)) + geom_point() + geom_line() +
  scale_color_gradientn(colors=hcl.colors(n=7, palette="Zissou 1"), limits=c(25, 80)) +
  scale_x_continuous(breaks = c(5, 10, 20), labels = c("0.2", "0.1", "0.05"), name="Variant allele frequency") +
  theme(aspect.ratio = 1)+ ggtitle("CD34, 90x")

ggplot(to.plot, aes(x=1/VAF, y=MutationCount, col=Age, group=paste(Paper_ID, Depth), 
                    linetype = factor(Depth))) + geom_point() + geom_line() +
  scale_color_gradientn(colors=hcl.colors(n=7, palette="Zissou 1"), limits=c(25, 80)) +
  scale_x_continuous(breaks = c(1/0.2, 1/0.1, 1/0.05, 1/0.02), labels = c("0.2", "0.1", "0.05", "0.02"), name="Variant allele frequency") +
  theme(aspect.ratio = 0.02/0.05)+ ggtitle("CD34, 90x and 270x")


dev.off()

####################################################################################################################################################
## Figure 4e: compare the physiological parameters in neutrally evolving samples estimates across the cohort

pdf(paste0(analysis.directory, "/Figures/Figure_4e.pdf"), width=5, height=w)

molten.par <- reshape2::melt(neutral.parameters[((grepl("N", neutral.parameters$Paper_ID) & neutral.parameters$Depth %in% c(90, 120)) |
                                                  (grepl("N", neutral.parameters$Paper_ID) & neutral.parameters$Depth %in% c(270, 300))) &
                                                  grepl("CD34", neutral.parameters$Sort),],
                             id.vars = c("Paper_ID", "lower", "upper", "Parameter", "Sort", "Depth"), value.name = "Median")

molten.par$Sort <- replace(molten.par$Sort, molten.par$Sort == "CD34_deep", "CD34")
molten.par$Paper_ID <- factor(molten.par$Paper_ID,
                                levels=unique(sample.info[molten.par$Paper_ID, ]$Paper_ID[order(sample.info[molten.par$Paper_ID, ]$Age)]))


## stem cell number
to.plot <- molten.par[molten.par$Parameter == "par_N",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # merge 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # merge 270x and 300x
## take 95% CI of lower and upper quantiles, measured with 90x, as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(molten.par[molten.par$Parameter=="par_N" & molten.par$Depth %in% c(90, 120),]$lower, p = 0.025), 2),
                                max = rep(quantile(molten.par[molten.par$Parameter=="par_N" & molten.par$Depth %in% c(90, 120),]$upper, p = 0.975), 2),
                                x=c(0, sum(molten.par$Parameter=="par_N"  & molten.par$Depth %in% c(90, 120))))

to.plot$x <- as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth))) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of stem cells (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "")


## N x tau
to.plot <- molten.par[molten.par$Parameter=="N_tau",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # merge 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # merge 270x and 300x
## take 95% CI of lower and upper quantiles, measured with 90x, as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(molten.par[molten.par$Parameter=="N_tau" & molten.par$Depth %in% c(90, 120),]$lower, p = 0.025), 2),
                                max = rep(quantile(molten.par[molten.par$Parameter=="N_tau" & molten.par$Depth %in% c(90, 120),]$upper, p = 0.975), 2),
                                x=c(0, sum( molten.par$Parameter=="N_tau" & molten.par$Depth %in% c(90, 120))))

to.plot$x <- as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth))) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="N x Tau", limits = c(1, 10^7)) +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "")

### Division rate
to.plot <- molten.par[molten.par$Parameter == "par_lambda_ss",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # merge 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # merge 270x and 300x
## take 95% CI of lower and upper quantiles, measured with 90x, as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(molten.par[molten.par$Parameter=="par_lambda_ss" & molten.par$Depth %in% c(90, 120),]$lower, p = 0.025), 2),
                                max = rep(quantile(molten.par[molten.par$Parameter=="par_lambda_ss" & molten.par$Depth %in% c(90, 120),]$upper, p = 0.975), 2),
                                x=c(0, sum(molten.par$Parameter=="par_lambda_ss" & molten.par$Depth %in% c(90, 120))))

to.plot$x <- as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, linetype=factor(Depth))) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = 365*10^min, ymax = 365*10^max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ 
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "")


### Mutation rate
to.plot <- molten.par[molten.par$Parameter == "par_mu",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # merge 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # merge 270x and 300x
## take 95% CI of lower and upper quantiles, measured with 90x, as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(molten.par[molten.par$Parameter=="par_mu" & molten.par$Depth %in% c(90, 120),]$lower, p = 0.025), 2),
                                max = rep(quantile(molten.par[molten.par$Parameter=="par_mu" & molten.par$Depth %in% c(90, 120),]$upper, p = 0.975), 2),
                                x=c(0, sum( molten.par$Parameter=="par_mu" & molten.par$Depth %in% c(90, 120))))

to.plot$x <- as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth))) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) + 
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "")

dev.off()


####################################################################################################################################################
## Figure 5c: compare the estimated clone size with the driver VAF

to.plot <- selected.parameters[selected.parameters$Parameter=="size_of_clone" &
                                 selected.parameters$Sort == "CD34" & 
                                 selected.paramters$Depth %in% c(90, 120),]
to.plot$CHIP.mutation <- apply(to.plot, 1, function(x){
  as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"])
})
to.plot <- to.plot[to.plot$CHIP.mutation!="no selected clone" & !is.na(to.plot$CHIP.mutation),]

# get ML estimate of the driver VAF and compute 95% CI based on Binomial distr.
to.plot$Driver.mean <- apply(to.plot, 1, function(x){
  
  type <- "CD34+"
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                            colnames(putative.drivers)[grep(type, colnames(putative.drivers))])]
  driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(x["Paper_ID"], "_"), colnames(driver.information))])]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = ncol(driver.information))
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  
  driver.information <- driver.information[driver.information$VAF < 0.9,]
  driver.information <- driver.information[which.max(driver.information$VAF),]
  vaf=as.numeric(driver.information["VAF"])
  if(is.na(vaf)){
    vaf <- 0
  }
  return(vaf)
})
to.plot$Driver.min <- apply(to.plot, 1, function(x){
  
  type <- "CD34+"
  sample.col <- colnames(putative.drivers)[grepl(x["Paper_ID"], colnames(putative.drivers)) & 
                                             grepl("CD34", colnames(putative.drivers)) &
                                             !grepl("deep", colnames(putative.drivers))]
  

  driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                            colnames(putative.drivers)[grep(type, colnames(putative.drivers))])]
  driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(x["Paper_ID"], "_"), colnames(driver.information))])]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = sample.col)
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = sample.col)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  
  driver.information <- driver.information[driver.information$VAF < 0.9,]
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

  type <- "CD34+"
  sample.col <- colnames(putative.drivers)[grepl(x["Paper_ID"], colnames(putative.drivers)) & 
                                             grepl("CD34", colnames(putative.drivers)) &
                                             !grepl("deep", colnames(putative.drivers)) ]

  driver.information <- putative.drivers[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                            colnames(putative.drivers)[grep(type, colnames(putative.drivers))])]
  driver.information <- driver.information[,c("CHROM", "POS", "REF", "ALT", "GENE", "AAchange",
                                              colnames(driver.information)[grep(paste0(x["Paper_ID"], "_"), colnames(driver.information))])]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = sample.col)
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = sample.col)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  
  driver.information <- driver.information[driver.information$VAF < 0.9,]
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

ggplot(to.plot[to.plot$Driver.mean!=0 & !grepl("N", to.plot$Paper_ID) &
                 !grepl("U", to.plot$Paper_ID),], aes(x=Median/2, xmin = lower/2, xmax = upper/2, color=CHIP.mutation,
                                                                                            y=Driver.mean, ymin = Driver.min, ymax=Driver.max)) + geom_point() +
  geom_errorbar() + geom_errorbarh() + scale_color_manual(values=CHIP.color) +
  geom_abline(slope = 1, intercept = 0, linetype=2) + 
  scale_x_continuous("VAF population genetics model") +
  scale_y_continuous("VAF WGS") + theme(aspect.ratio = 1)

dev.off()


####################################################################################################################################################
## Fig. 5d, e: Physiological parameters in selected cases

pdf(paste0(analysis.directory, "/Figures/Figure_5_d_e.pdf"), width=5, height=3.5)

molten.par <- reshape2::melt(selected.parameters[((selected.parameters$Paper_ID %in% clearly.selected.samples.90 & selected.parameters$Depth %in% c(90, 120)) |
                                                   (selected.parameters$Paper_ID %in% clearly.selected.samples.270 & selected.parameters$Depth %in% c(270, 300))) &
                                                   grepl("CD34", selected.parameters$Sort),],
                             id.vars = c("Paper_ID", "lower", "upper", "Parameter", "Sort", "Depth"), value.name = "Median")
molten.par$CHIP.mutation <- sample.info[molten.par$Paper_ID,]$CHIP.mutation.associated.with.fit
molten.par$CHIP.mutation <- factor(molten.par$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))
molten.par$Sort <- replace(molten.par$Sort, molten.par$Sort == "CD34_deep", "CD34")
molten.par$Paper_ID <- factor(sample.info[molten.par$Paper_ID, ]$Paper_ID,
                                levels=unique(sample.info[molten.par$Paper_ID, ]$Paper_ID[order(sample.info[molten.par$Paper_ID, ]$Age)]))

## stem cell number
to.plot <- molten.par[molten.par$Parameter == "par_N",]

to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # combine 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # combine 270x and 300x

## take 95% CI of lower and upper quantiles of the neutral cases as across the cohort for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_N" & neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", neutral.parameters$Paper_ID),]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_N" & neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", neutral.parameters$Paper_ID),]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Sort=="CD34" & to.plot$Parameter=="par_N"  & to.plot$Depth %in% c(90, 120))))

to.plot$x <- as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth), col = CHIP.mutation)) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of stem cells (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "")


## N x tau

to.plot <- molten.par[molten.par$Parameter == "N_tau",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # combine 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # combine 270x and 300x

## take 95% CI of lower and upper quantiles of the neutral cases as across the cohort for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="N_tau" & neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", neutral.parameters$Paper_ID),]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="N_tau" & neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", neutral.parameters$Paper_ID),]$upper, p = 0.975), 2),
                                x=c(0, sum(to.plot$Parameter=="N_tau"  & to.plot$Depth %in% c(90, 120))))

to.plot$x <- as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth), col = CHIP.mutation)) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="N x Tau", limits = c(1, 10^8)) +
  theme(aspect.ratio = 1)+ 
  scale_x_continuous(breaks = unique(as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "")

# Division rate
to.plot <- molten.par[molten.par$Parameter == "par_lambda_ss",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # combine 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # combine 270x and 300x

## take 95% CI of lower and upper quantiles of the neutral cases as across the cohort for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss" & neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", neutral.parameters$Paper_ID),]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss" & neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", neutral.parameters$Paper_ID),]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Parameter=="par_lambda_ss" & to.plot$Depth %in% c(90, 120))))

to.plot$x <- as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, linetype=factor(Depth), col = CHIP.mutation)) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = 365*10^min, ymax = 365*10^max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ 
  scale_x_continuous(breaks = unique(as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "")


### Mutation rate
to.plot <- molten.par[molten.par$Parameter == "par_mu",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # combine 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # combine 270x and 300x

## take 95% CI of lower and upper quantiles of the neutral cases as across the cohort for comparison

to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_mu" & neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", neutral.parameters$Paper_ID),]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_mu" & neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", neutral.parameters$Paper_ID),]$upper, p = 0.975), 2),
                                x=c(0, sum(to.plot$Parameter=="par_mu" & to.plot$Depth %in% c(90, 120) )))

to.plot$x <- as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth), col = CHIP.mutation)) + geom_pointrange(fatten = 0.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x= 0, y = 0)+
  scale_x_continuous(breaks = unique(as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "")

dev.off()


####################################################################################################################################################
## Fig. 5f-h: Compare physiological parameters between neutrally evolving cases and cases with selection

# where available, use 270x estimates and else the 90x ones
to.plot <- rbind(neutral.parameters[grepl("N", neutral.parameters$Paper_ID) & neutral.parameters$Depth %in% c(90, 120) &
                                      sample.info[neutral.parameters$Paper_ID,]$`Coverage.WGS.CD34+.2`==0,],
                 neutral.parameters[grepl("N", neutral.parameters$Paper_ID) & neutral.parameters$Depth %in% c(270, 300),])
to.plot <- rbind(to.plot, selected.parameters[ !grepl("N", selected.parameters$Paper_ID) & !grepl("U", selected.parameters$Paper_ID) & selected.parameters$Depth %in% c(90, 120) &
                                                 sample.info[neutral.parameters$Paper_ID,]$`Coverage.WGS.CD34+.2`==0,],
                 selected.parameters[ !grepl("N", selected.parameters$Paper_ID) & !grepl("U", selected.parameters$Paper_ID) & selected.parameters$Depth %in% c(270, 300),])

# subset on CD34
to.plot <- to.plot[grepl("CD34", to.plot$Sort),]

to.plot$CHIP.mutation <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit

# annotate age of the individuals
to.plot$Age <- sample.info[to.plot$Paper_ID,]$Age
## make ages unique to facilitate visualization
to.plot <- to.plot[order(to.plot$Age),]
to.plot$Age <- unlist(sapply(unique(to.plot$Age), function(x){
  tmp <- to.plot[to.plot$Age==x,]
  for(i in unique(tmp$Paper_ID)){
    tmp[tmp$Paper_ID==i,]$Age <- tmp[tmp$Paper_ID==i,]$Age + (which(unique(tmp$Paper_ID)==i)-1)/2
  }
  tmp$Age
}))

to.plot$Type <- ifelse(grepl("N", to.plot$Paper_ID), "N", "Selection")
to.plot$Type <- factor(to.plot$Type, levels=c("Neutral", "Selection"))


pdf(paste0(analysis.directory, "Figures/Figure_5_f_g_h.pdf"), width=5, height=2.5)

## Mutation rate
ggplot(to.plot[to.plot$Parameter == "par_mu",],
       aes(x=Type, y=Median, ymin=lower, ymax=upper, col = Type)) + geom_boxplot(outlier.shape = NA) +
  geom_pointrange(position = position_dodge2(width=0.5)) + scale_color_manual(values=c(Neutral="darkgrey", Selection="red"))  +
  scale_y_continuous(name ="Number of SSNVs per division") +
  expand_limits(y = 0) 

## Stem cell number
ggplot(to.plot[to.plot$Parameter == "par_N",],
       aes(x=Type, y=Median, ymin=lower, ymax=upper, col = Type)) + geom_boxplot(outlier.shape = NA) +
  geom_pointrange(position = position_dodge2(width=0.5)) + scale_color_manual(values=c(Neutral="darkgrey", Selection="red"))  +
  scale_y_continuous(name ="Stem cell number") +
  expand_limits(y = 0)


## Division rate
ggplot(to.plot[to.plot$Parameter=="par_lambda_ss",],
       aes(x=Age, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col = Type)) +
  geom_point() + geom_errorbar() + scale_color_manual(values=c(Neutral="darkgrey", Selection="red"))  +
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Division rate (1/y)") +
  expand_limits(x = 0) 

ggplot(to.plot[to.plot$Parameter == "par_lambda_ss"],
       aes(x=Type, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col = Type)) + geom_boxplot(outlier.shape = NA) +
  geom_pointrange(position = position_dodge2(width=0.5)) + scale_color_manual(values=c(Neutral="darkgrey", Selection="red"))  +
  scale_y_log10(name ="Division rate (1/y)") 

dev.off()


####################################################################################################################################################
## Figure 6b, size of the selected clones in cases with unknown drivers

molten.par <- reshape2::melt(selected.parameters[grepl("CD34", selected.parameters$Sort) & grepl("U", selected.parameters$Paper_ID) &
                                                   ((selected.parameters$Paper_ID %in% selection.no.driver.90 & selected.parameters$Depth %in% c(90, 120)) |
                                                   (selected.parameters$Paper_ID %in% selection.no.driver.270 & selected.parameters$Depth %in% c(270, 300))),],
                             id.vars = c("Paper_ID", "lower", "upper", "Parameter", "Sort", "Depth"), value.name = "Median")
molten.par$Paper_ID <- factor(molten.par$Paper_ID,
                                levels=unique(molten.par$Paper_ID[order(sample.info[molten.par$Paper_ID, ]$Age)]))


pdf(paste0(analysis.directory, "Figures/Figure_6b.pdf"), width=3.5, height=2.5)

## Subset on the size of the clone
to.plot <- molten.par[molten.par$Parameter == "size_of_clone",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # combine 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # combine 270x and 300x

to.plot$x <- as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth), col = factor(Depth))) + geom_pointrange(fatten = 1.5) +
  geom_pointrange(fatten = 1.5) + scale_color_manual(values=c("90" = "black", "270" = "grey")) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Size of clone") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y =0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.3 - 1),
                     labels = (unique(as.character(to.plot$Paper_ID))), name = "")

dev.off()

####################################################################################################################################################
## Fig. 6d: Physiological parameters in selected cases

molten.par <- reshape2::melt(rbind(selected.parameters[(grepl("CD34", selected.parameters$Sort) & selected.parameters$Paper_ID %in% selection.no.driver.90 & selected.parameters$Depth %in% c(90, 120)) |
                                                         (grepl("CD34", selected.parameters$Sort) & selected.parameters$Paper_ID %in% selection.no.driver.270 & selected.parameters$Depth %in% c(270, 300)),],
                                   neutral.parameters[grepl("CD34", neutral.parameters$Sort) & neutral.parameters$Paper_ID %in% setdiff(selection.no.driver.270, selection.no.driver.90) & neutral.parameters$Depth %in% c(90, 120),]),
                             id.vars = c("Paper_ID", "lower", "upper", "Parameter", "Sort", "Depth"), value.name = "Median")
molten.par$CHIP.mutation <- sample.info[molten.par$Paper_ID,]$CHIP.mutation.associated.with.fit

molten.par$Sort <- replace(molten.par$Sort, molten.par$Sort == "CD34_deep", "CD34")
molten.par$Paper_ID <- factor(molten.par$Paper_ID,
                                levels=unique(molten.par$Paper_ID[order(sample.info[molten.par$Paper_ID, ]$Age)]))


pdf(paste0(analysis.directory, "/Figures/Figure_6_d.pdf"), width=3.5, height=3.5)

## stem cell number
to.plot <- molten.par[molten.par$Parameter == "par_N",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # combine 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # combine 270x and 300x
## take 95% CI of lower and upper quantiles of the neutrally evolving cases for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_N"& neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", sample.info$Paper_ID),]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_N"& neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) & grepl("N", sample.info$Paper_ID),]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Parameter=="par_N" & to.plot$Depth %in% c(90, 120) & grepl("U", to.plot$Paper_ID))))

to.plot$x <- as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth))) + geom_pointrange(fatten = 1.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of stem cells (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.3 - 1),
                     labels = (unique(as.character(to.plot$Paper_ID))), name = "")


### Division rate
to.plot <- molten.par[molten.par$Parameter == "par_lambda_ss",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # combine 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # combine 270x and 300x
## take 95% CI of lower and upper quantiles of the neutrally evolving cases for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss"& neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) &  grepl("N", neutral.parameters$Paper_ID),]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss"& neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) &  grepl("N", neutral.parameters$Paper_ID),]$upper, p = 0.975), 2),
                                x=c(0, sum(  to.plot$Parameter=="par_lambda_ss" & to.plot$Depth %in% c(90, 120) & grepl("U", to.plot$Paper_ID))))

to.plot$x <- as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, linetype=factor(Depth))) + geom_pointrange(fatten = 1.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = 365*10^min, ymax = 365*10^max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ 
  scale_x_continuous(breaks = unique(as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.3 - 1),
                     labels = (unique(as.character(to.plot$Paper_ID))), name = "")

### Mutation rate
to.plot <- molten.par[molten.par$Parameter == "par_mu",]
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==120, 90) # combine 90x and 120x
to.plot$Depth <- replace(to.plot$Depth, to.plot$Depth==300, 270) # combine 270x and 300x
## take 95% CI of lower and upper quantiles of the neutrally evolving cases for comparison
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_mu"& neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) &  grepl("N", neutral.parameters$Paper_ID),]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Parameter=="par_mu"& neutral.parameters$Sort=="CD34" & neutral.parameters$Depth %in% c(90, 120) &  grepl("N", neutral.parameters$Paper_ID),]$upper, p = 0.975), 2),
                                x=c(0, sum( to.plot$Sort=="CD34" & to.plot$Parameter=="par_mu" & to.plot$Depth %in% c(90, 120) & grepl("U", to.plot$Paper_ID) )))

to.plot$x <- as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.2* 
  as.numeric(factor(to.plot$Depth, levels = c(90, 270))) -1

ggplot(to.plot,
       aes(x=x, y=Median, ymin=lower, ymax=upper, linetype=factor(Depth))) + geom_pointrange(fatten = 1.5) +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x-0.2), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange(fatten = 1.5) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(as.character(to.plot$Paper_ID), levels = sort(unique(as.character(to.plot$Paper_ID)))))+0.3 - 1),
                     labels = (unique(as.character(to.plot$Paper_ID))), name = "")


dev.off()

####################################################################################################################################################
## In Figure 7, Extended Data Fig. 8g and Extended Data Fig.9, we show parameter estimates obtained with both the one-clone model (see code above)
## and the 2-clone model. To load the parameters obtained with the 2-clone model, we source the file "Plot_fits_WGS_data_2_clone_model.R"

source(paste0(custom.script.directory, "/Analysis_and_figures/Plot_fits_WGS_data_2_clone_model.R"))


## cases with evidence for a second selected clone, parameters according to a two-clone model
two.clone.model.parameters <- selected.parameters.2cm.2clones[paste(selected.parameters.2cm.2clones$Paper_ID, selected.parameters.2cm.2clones$Mode, sep="_") %in%
                                                      colnames(model.support.selection.2_clones[,model.support.selection.2_clones["Clone 2",] > 15]),]
two.clone.model.parameters$Subclones <- "2 clones"


## cases with evidence for selection, parameters according to a one-clone model
one.clone.model.parameters <- selected.parameters[!grepl("N", selected.parameters$Paper_ID) &
                                              selected.parameters$Sort %in% c("CD34", "CD34_deep"),]
one.clone.model.parameters$Subclones <- "1 clone"
one.clone.model.parameters$Mode <- ""


####################################################################################################################################################
## Extended Data Fig. 8g, compare the physiological parameter estimates obtained with the 1-clone and the 2-clone model

## Collect the physiological parameters for all cases in one data frame; for the 2 clone model take the average between linear and branched evolution, whenever both topologies suggested selection
                                                                                                                                                                                   
parameters.2clones <- two.clone.model.parameters[two.clone.model.parameters$Parameter %in% c("par_N", "par_lambda_ss", "par_mu"),intersect(colnames(two.clone.model.parameters),
                                                                                                                                 colnames(one.clone.model.parameters))]
parameters.1clone <- one.clone.model.parameters[one.clone.model.parameters$Paper_ID %in% two.clone.model.parameters$Paper_ID &
                                                  one.clone.model.parameters$Parameter %in% c("par_N", "par_lambda_ss", "par_mu") &
                                                  one.clone.model.parameters$Sort == "CD34_deep",intersect(colnames(two.clone.model.parameters),
                                                                                                                                          colnames(one.clone.model.parameters))]

## for plotting, average the two 2-clone topologies, in case both topologies fit the data

to.plot <- data.frame()

for(i in unique(parameters.2clones$Paper_ID)){
  for(j in c("par_mu", "par_N", "par_lambda_ss")){
    tmp.2clones <- parameters.2clones[parameters.2clones$Paper_ID == i &
                                        grepl(j, parameters.2clones$Parameter),]
    
    tmp.1clone <- parameters.1clone[parameters.1clone$Paper_ID == i &
                                      grepl(j, parameters.1clone$Parameter),]
    
    to.plot <- rbind(to.plot,
                     data.frame(lower.1clone = mean(tmp.1clone$lower),
                                upper.1clone = mean(tmp.1clone$upper),
                                Parameter = j,
                                Median.1clone = mean(tmp.1clone$Median),
                                Paper_ID = i,
                                Subclones = tmp.2clones$Subclones[1],
                                lower.2clones = mean(tmp.2clones$lower),
                                upper.2clones = mean(tmp.2clones$upper),
                                Median.2clones = mean(tmp.2clones$Median)))
    
  }
}

rm(parameters.1clone)
rm(parameters.2clones)

pdf(paste0(analysis.directory, "/Figures/Figure_S8g.pdf"), width=3.5, height=2.5)

## stem cell number
ggplot(to.plot[to.plot$Parameter == "par_N",], 
       aes(x = Median.1clone, xmin = lower.1clone, xmax = upper.1clone,
           y = Median.2clones, ymin = lower.2clones, ymax = upper.2clones)) + geom_point() +
  geom_errorbarh(height = 0) + geom_errorbar(width = 0) + theme(aspect.ratio = 1) +
  scale_x_continuous("Stem cell number (1 clone; log10)") + scale_y_continuous("Stem cell number (2 clones; log10)") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + expand_limits(x = 0, y = 0) 

## division rate
ggplot(to.plot[to.plot$Parameter == "par_lambda_ss",], 
       aes(x = 365*10^Median.1clone, xmin = 365*10^lower.1clone, xmax = 365*10^upper.1clone,
           y = 365*10^Median.2clones, ymin = 365*10^lower.2clones, ymax = 365*10^upper.2clones)) + geom_point() +
  geom_errorbarh(height = 0) + geom_errorbar(width = 0) + theme(aspect.ratio = 1) +
  scale_x_log10("Stem cell divisions per year (1 clone)") + scale_y_log10("Stem cell divisions per year (2 clones)") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + expand_limits(x = 0, y = 0) 

## mutation rate
ggplot(to.plot[to.plot$Parameter == "par_mu",], 
       aes(x = Median.1clone, xmin = lower.1clone, xmax = upper.1clone,
           y = Median.2clones, ymin = lower.2clones, ymax = upper.2clones)) + geom_point() +
  geom_errorbarh(height = 0) + geom_errorbar(width = 0) + theme(aspect.ratio = 1) +
  scale_x_continuous("Number of SSNVs per division (1 clone)") + scale_y_continuous("Number of SSNVs per division (2 clones)") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + expand_limits(x = 0, y = 0) 

dev.off()

####################################################################################################################################################
## Extended Data Fig. 9a/c: compare age and growth rate of the leading (largest) selected clone between the one-clone model (90x/270x) and the 2-clone model (270x)

parameters.2clones <- two.clone.model.parameters[two.clone.model.parameters$Parameter %in% c("age_of_clone_1", "growth_per_year1"),]

## for plotting, average the two 2-clone topologies, in case both topologies fit the data

to.plot <- data.frame()

for(i in unique(parameters.2clones$Paper_ID)){
  for(j in c("age_of_clone", "growth_per_year")){
    tmp <- parameters.2clones[parameters.2clones$Paper_ID == i &
                                grepl(j, parameters.2clones$Parameter),]
    
    if(nrow(tmp)==0){next}
    to.plot <- rbind(to.plot,
                     data.frame(lower = mean(tmp$lower),
                                upper = mean(tmp$upper),
                                Parameter = j,
                                Median = mean(tmp$Median),
                                Paper_ID = i,
                                Subclones = tmp$Subclones[1]))
    
  }
}

to.plot$Depth <- 270

rm(parameters.2clones)

parameters.1clone <- one.clone.model.parameters[one.clone.model.parameters$Parameter %in% c("age_of_clone", "growth_per_year"),]

parameters.1clone$Subclones <- "1 clone"

to.plot <- rbind(to.plot, parameters.1clone[,colnames(to.plot)])

rm(parameters.1clone)


pdf(paste0(analysis.directory, "/Figures/Figure_S9ab.pdf"), width=5, height=2.5)

to.plot$CHIP.mutation <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))
to.plot$CHIP.mutation.2 <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit.2
to.plot$CHIP.mutation.2 <- factor(to.plot$CHIP.mutation.2, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))

## Age of leading clone 
to.plot. <- to.plot[to.plot$Parameter == "age_of_clone",]
to.plot.$Depth <- replace(to.plot.$Depth, to.plot.$Depth==120, 90) ## merge 90x and 120x
to.plot.$Depth <- replace(to.plot.$Depth, to.plot.$Depth==300, 270) ## merge 270x and 300x
to.plot.$x <- as.numeric(factor(to.plot.$Paper_ID, levels = sort(unique(to.plot.$Paper_ID))))+0.2* 
  as.numeric(factor(to.plot.$Depth, levels = c(90, 270))) + 0.2*
  as.numeric(factor(to.plot.$Subclones, levels = c("1 clone", "2 clones"))) -1 - 0.2
to.plot.$Linetype <- factor(paste(to.plot.$Depth, to.plot.$Subclones), levels = c("90 1 clone", "270 1 clone", "270 2 clones"))

ggplot(to.plot.,
       aes(x=x, y=Median, ymin=lower, ymax=upper, col = CHIP.mutation, linetype = Linetype, shape = Linetype)) + 
  geom_pointrange(fatten = 2) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of leading clone (years)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot.$Paper_ID, levels = sort(unique(to.plot.$Paper_ID))))+0.4 - 1),
                     labels = (unique(to.plot.$Paper_ID)), name = "")


## Clonal growth (leading clone)
to.plot. <- to.plot[to.plot$Parameter == "growth_per_year",]
to.plot.$Depth <- replace(to.plot.$Depth, to.plot.$Depth==120, 90) ## merge 90x and 120x
to.plot.$Depth <- replace(to.plot.$Depth, to.plot.$Depth==300, 270) ## merge 270x and 300x
to.plot.$x <- as.numeric(factor(to.plot.$Paper_ID, levels = sort(unique(to.plot.$Paper_ID))))+0.2* 
  as.numeric(factor(to.plot.$Depth, levels = c(90, 270))) + 0.2*
  as.numeric(factor(to.plot.$Subclones, levels = c("1 clone", "2 clones"))) -1 - 0.2
to.plot.$Linetype <- factor(paste(to.plot.$Depth, to.plot.$Subclones), levels = c("90 1 clone", "270 1 clone", "270 2 clones"))

ggplot(to.plot.,
       aes(x=x, y=100*Median, ymin=100*lower, ymax=100*upper, col = CHIP.mutation,
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Growth per year leading clone (%)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot.$Paper_ID, levels = sort(unique(to.plot.$Paper_ID))))+0.4 - 1),
                     labels = (unique(to.plot.$Paper_ID)), name = "")

dev.off()


####################################################################################################################################################
## Extended Data Fig. 9b/d: compare age and growth rate of the minor selected clone across the cohort

selection.parameters.second.clone <- two.clone.model.parameters[two.clone.model.parameters$Parameter %in% c("growth_per_year2", "age_of_clone_2",),]

## for plotting, average the two 2-clone topologies, in case both topologies fit the data

to.plot <- data.frame()

for(i in unique(selection.parameters.second.clone$Paper_ID)){
  for(j in unique(selection.parameters.second.clone$Parameter)){
    tmp <- selection.parameters.second.clone[selection.parameters.second.clone$Paper_ID == i &
                                               selection.parameters.second.clone$Parameter == j,]
    to.plot <- rbind(to.plot,
                     data.frame(lower = mean(tmp$lower),
                                upper = mean(tmp$upper),
                                Parameter = j,
                                Median = mean(tmp$Median),
                                Paper_ID = i,
                                Subclones = tmp$Subclones[1]))
  }
}


pdf(paste0(analysis.directory, "/Figures/Figure_S9bd.pdf"), width=5, height=2)

to.plot$CHIP.mutation <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))
to.plot$CHIP.mutation.2 <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit.2
to.plot$CHIP.mutation.2 <- factor(to.plot$CHIP.mutation.2, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))

## Age of clone 2
to.plot. <- to.plot[to.plot$Parameter %in% c("age_of_clone_2") & !is.na(to.plot$Median) ,]
to.plot.$Paper_ID <- factor(to.plot.$Paper_ID, levels = sort(unique(as.character(to.plot.$Paper_ID))))

ggplot(to.plot.,
       aes(x=Paper_ID, y=Median, ymin=lower, ymax=upper, col = CHIP.mutation.2)) + geom_pointrange() +
  geom_pointrange() + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of second clone (years)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) 


## Clonal growth (clone 2)
to.plot. <- to.plot[to.plot$Parameter %in% c("growth_per_year2") ,]
to.plot.$Paper_ID <- factor(to.plot.$Paper_ID, levels = sort(unique(as.character(to.plot.$Paper_ID))))

ggplot(to.plot.[!is.na(to.plot.$Median),],
       aes(x=Paper_ID, y=100*Median, ymin=100*lower, ymax=100*upper, col = CHIP.mutation.2)) + geom_pointrange() +
  geom_pointrange() + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Clonal growth per year (2nd clone) (%)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) 

dev.off()


####################################################################################################################################################
## Figure 7a-c compare age and growth across the cohort; combining estimates from one- and two-clone model:
## where there is evidence for a second selected clone, use the estimates from the 2-clone model
## where there is evidence for only 1 selected clone but 270x WGS data available, use the estimates from the 1-clone model on 270x data
## where there is evidence for only 1 selected clone and no 270x WGS data available, use the estimates from the 1-clone model on 90x data

selection.parameters.combined <- rbind(two.clone.model.parameters[two.clone.model.parameters$Parameter %in% c("growth_per_year1", "growth_per_year2",
                                                                                                     "age_of_clone_1", "age_of_clone_2"),intersect(colnames(two.clone.model.parameters),
                                                                                                                                                              colnames(one.clone.model.parameters))],
                                          one.clone.model.parameters[!one.clone.model.parameters$Paper_ID %in% two.clone.model.parameters &
                                                                       one.clone.model.parameters$Sort == "CD34_deep" & 
                                            one.clone.model.parameters$Parameter %in% c("growth_per_year", "age_of_clone"),intersect(colnames(two.clone.model.parameters),colnames(one.clone.model.parameters))],
                                          one.clone.model.parameters[!one.clone.model.parameters$Paper_ID %in% two.clone.model.parameters &
                                                                       sample.info[one.clone.model.parameters$Paper_ID,]$`Coverage.WGS.CD34+.2`==0 &
                                                                       one.clone.model.parameters$Sort == "CD34" & 
                                                                       one.clone.model.parameters$Parameter %in% c("growth_per_year", "age_of_clone"),intersect(colnames(two.clone.model.parameters),colnames(one.clone.model.parameters))])

## for plotting, average the two 2-clone topologies, in case both topologies fit the data

to.plot <- data.frame()

for(i in unique(selection.parameters.combined$Paper_ID)){
  for(j in unique(selection.parameters.combined$Parameter)){
    tmp <- selection.parameters.combined[selection.parameters.combined$Paper_ID == i &
                                           selection.parameters.combined$Parameter == j,]
    to.plot <- rbind(to.plot,
                     data.frame(lower = mean(tmp$lower),
                                upper = mean(tmp$upper),
                                Parameter = j,
                                Median = mean(tmp$Median),
                                Paper_ID = i,
                                Subclones = tmp$Subclones[1]))
  }
}

pdf(paste0(analysis.directory, "/Figures/Figure_7abc.pdf"), width=5, height=2)

to.plot$CHIP.mutation <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))
to.plot$CHIP.mutation.2 <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit.2
to.plot$CHIP.mutation.2 <- factor(to.plot$CHIP.mutation.2, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))

## Age of leading clone 
to.plot. <- to.plot[to.plot$Parameter %in% c("age_of_clone_1", "age_of_clone") ,]

## summarize by driver
ggplot(to.plot.,
       aes(x=CHIP.mutation, y=Median, ymin=lower, ymax=upper, col = CHIP.mutation)) + geom_boxplot(outliers = F) +
  geom_pointrange(position = position_jitter(width = 0.25)) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=30)) + scale_y_continuous(name ="Age of first clone (years)") +
  expand_limits(x = 0, y = 0) 

## Age of minor clone
to.plot. <- to.plot[to.plot$Parameter %in% c("age_of_clone_2") & !is.na(to.plot$Median) ,]

## summarize by driver
ggplot(to.plot.,
       aes(x=CHIP.mutation.2, y=Median, ymin=lower, ymax=upper, col = CHIP.mutation.2)) + geom_boxplot(outliers = F) +
  geom_pointrange(position = position_jitter(width = 0.25)) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=30)) + scale_y_continuous(name ="Age of second clone (years)") +
  expand_limits(x = 0, y = 0) 

## Clonal growth (both clones)
to.plot. <- to.plot[to.plot$Parameter %in% c("growth_per_year1", "growth_per_year", "growth_per_year2") ,]
to.plot.$Paper_ID <- factor(to.plot.$Paper_ID, levels = sort(unique(as.character(to.plot.$Paper_ID))))
to.plot.$CHIP.mutation.this.clone <- ifelse(to.plot.$Parameter %in% c("growth_per_year", "growth_per_year1"), 
                                            as.character(to.plot.$CHIP.mutation), as.character(to.plot.$CHIP.mutation.2))

## summarize by driver
ggplot(to.plot.[!is.na(to.plot.$Median),],
       aes(x=CHIP.mutation.this.clone, y=100*Median, ymin=100*lower, ymax=100*upper, col = CHIP.mutation.this.clone)) + geom_boxplot(outliers = F) +
  geom_pointrange(position = position_jitter(width = 0.25)) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=30)) + scale_y_log10(name ="Clonal growth per year (both clones) (%)") +
  expand_limits(x = 0, y = 0) 

dev.off()



####################################################################################################################################################
## Fig. 7d: plot the incidence of CH driver acquisition
## where there is evidence for a second selected clone, use the estimates from the 2-clone model
## where there is evidence for only 1 selected clone but 270x WGS data available, use the estimates from the 1-clone model on 270x data
## where there is evidence for only 1 selected clone and no 270x WGS data available, use the estimates from the 1-clone model on 90x data


to.plot. <- rbind(two.clone.model.parameters[two.clone.model.parameters$Parameter %in% c("par_ts1", "par_ts2"),intersect(colnames(two.clone.parameters),
                                                                                                                         colnames(one.clone.model.parameters))],
                  one.clone.model.parameters[one.clone.model.parameters$Parameter == "par_t_s_absolute" &
                                               !one.clone.model.parameters$Paper_ID %in% two.clone.model.parameters$Paper_ID &
                                               one.clone.model.parameters$Sort == "CD34_deep", intersect(colnames(two.clone.parameters),
                                                                                                         colnames(one.clone.model.parameters))],
                  one.clone.model.parameters[one.clone.model.parameters$Parameter == "par_t_s_absolute" &
                                               !one.clone.model.parameters$Paper_ID %in% two.clone.model.parameters$Paper_ID &
                                               sample.info[one.clone.model.parameters$Paper_ID,]$`Coverage.WGS.CD34+.2`==0 &
                                               one.clone.model.parameters$Sort == "CD34", intersect(colnames(two.clone.parameters), colnames(one.clone.model.parameters.neutral))])

## for plotting, average the two 2-clone topologies, in case both topologies fit the data

to.plot.driver <- data.frame()

for(i in unique(to.plot.$Paper_ID)){
  for(j in unique(to.plot.$Parameter)){
    tmp <- to.plot.[to.plot.$Paper_ID == i &
                      to.plot.$Parameter == j,]
    if(nrow(tmp)==0){next}
    to.plot.driver <- rbind(to.plot.driver,
                            data.frame(lower = mean(tmp$lower),
                                       upper = mean(tmp$upper),
                                       Parameter = j,
                                       Median = mean(tmp$Median),
                                       Paper_ID = i,
                                       Subclones = tmp$Subclones[1]))
  }
}

## driver info:
to.plot.driver <- to.plot.driver[order(to.plot.driver$Median),]
to.plot.driver$y <- (1:nrow(to.plot.driver))/nrow(to.plot.driver)
to.plot.driver$CHIP.mutation <- apply(to.plot.driver, 1, function(x){
  if(x["Parameter"] != "par_ts2"){
    as.character(sample.info[x["Paper_ID"],]$CHIP.mutation.associated.with.fit)
  }else{
    as.character(sample.info[x["Paper_ID"],]$CHIP.mutation.associated.with.fit.2)
  }
})
to.plot.driver$CHIP.mutation <- factor(to.plot.driver$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))

## cohort summary
to.plot. <- data.frame(Age = seq(0, ceiling(max(to.plot.driver$upper)/365/25)*25),
                       Incidence = sapply(seq(0,ceiling(max(to.plot.driver$upper)/365/25)*25), function(a){sum(to.plot.driver$Median/365<=a)})/nrow(to.plot.driver),
                       Lower = sapply(seq(0, ceiling(max(to.plot.driver$upper)/365/25)*25), function(a){sum(to.plot.driver$upper/365<=a)})/nrow(to.plot.driver),
                       Upper = sapply(seq(0,ceiling(max(to.plot.driver$upper)/365/25)*25), function(a){sum(to.plot.driver$lower/365<=a)})/nrow(to.plot.driver))


pdf(paste0(analysis.directory, "/Figures/Figure_7d.pdf"), width=3.5, height=2.5)

ggplot(to.plot., aes(x=Age, y=Incidence, ymin=Lower, ymax=Upper)) + geom_ribbon(fill="grey", alpha=0.5) + geom_line() +
  geom_point(data=to.plot.driver, aes(x=Median/365, y=y, col=CHIP.mutation), inherit.aes = F) +
  geom_errorbarh(data=to.plot.driver, aes(xmin=lower/365, xmax=upper/365, col=CHIP.mutation, y=y), inherit.aes = F, height=0) +
  scale_color_manual(values=CHIP.color) + scale_x_continuous(breaks = seq(0, 75, 25)) + theme(aspect.ratio = 1)

dev.off()


####################################################################################################################################################
## In Extended Data Fig. 7e we compare parameter estimates obtained with the one-clone model assuming size compensation for the selected clone and the one-clone model without size-compensation
## To load the parameters obtained without size compensation, we source the file "Plot_fits_WGS_data_no_size_compensation.R"

source(paste0(custom.script.directory, "/Analysis_and_figures/Plot_fits_WGS_data_no_size_compensation.R"))

## parameters obtained without size compensation
tmp.1 <- selected.parameters.nsc
tmp.1$Size_compensation <- "no"

## parameters obtained with size compensation
tmp.2 <- selected.parameters
tmp.2$Size_compensation <- "yes"

to.plot <- rbind(tmp.2[tmp.2$Paper_ID %in% tmp.1$Paper_ID & tmp.2$Sort == "CD34" & tmp.2$Depth %in% c(90, 120),],
                 tmp.1)

rm(tmp.1, tmp.2)


pdf(paste0(analysis.directory, "/Figures/Figure_S7e.pdf"), width=3.5, height=2.5)


to.plot$x <- as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.2* 
  as.numeric(factor(to.plot$Size_compensation, levels = c("yes", "no"))) -1

## Stem cell number

ggplot(to.plot[to.plot$Parameter=="par_N",],
       aes(x=x, y=Median, ymin=lower, ymax=upper, col=factor(Size_compensation, levels = c("yes", "no")))) + 
  geom_pointrange(fatten = 1.5) + scale_color_manual(values = c(yes = "black", no = "orange")) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of stem cells (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "") +
  guides(color=guide_legend(title="Size compensation"))


## Division rate

ggplot(to.plot[to.plot$Parameter=="par_lambda_ss",],
       aes(x=x, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col=factor(Size_compensation, levels = c("yes", "no")))) + 
  geom_pointrange(fatten = 1.5) + scale_color_manual(values = c(yes = "black", no = "orange")) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "") +
  guides(color=guide_legend(title="Size compensation"))


## Mutation rate

ggplot(to.plot[to.plot$Parameter=="par_mu",],
       aes(x=x, y=Median, ymin=lower, ymax=upper, col=factor(Size_compensation, levels = c("yes", "no")))) + 
  geom_pointrange(fatten = 1.5) + scale_color_manual(values = c(yes = "black", no = "orange")) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "") +
  guides(color=guide_legend(title="Size compensation"))

## Age of selected clone

ggplot(to.plot[to.plot$Parameter=="age_of_clone",],
       aes(x=x, y=Median, ymin=lower, ymax=upper, col=factor(Size_compensation, levels = c("yes", "no")))) + 
  geom_pointrange(fatten = 1.5) + scale_color_manual(values = c(yes = "black", no = "orange")) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of clone (years)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "") +
  guides(color=guide_legend(title="Size compensation"))


## Growth of selected clone

ggplot(to.plot[to.plot$Parameter=="growth_per_year",],
       aes(x=x, y=100*Median, ymin=100*lower, ymax=100*upper, col=factor(Size_compensation, levels = c("yes", "no")))) + 
  geom_pointrange(fatten = 1.5) + scale_color_manual(values = c(yes = "black", no = "orange")) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Clonal growth per year (%)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot$Paper_ID, levels = sort(unique(to.plot$Paper_ID))))+0.3 - 1),
                     labels = sort(unique(to.plot$Paper_ID)), name = "") +
  guides(color=guide_legend(title="Size compensation"))


dev.off()



