###### This scripts plots the output from the population genetics model for previously published data. It produces
###### - For each patient the model fit vs data and the posterior probabilities of the parameters
###### - A summary of the 80%HDI estimates across the population
############################################################################################################################################
###### libraries and functions
library(cdata)
library(scales)
library(openxlsx)
library(FLORENCE)

source("~/Settings.R")

patient.ids <-sample.info.published.data$SAMPLE

seq.type="sc"

############################################################################################################################################
## Summary parameters and individual model fits (as shown in Fig. 3/S2)

parameters <- data.frame()
neutral.parameters <- data.frame()
selected.parameters <- data.frame()
evidence.for.subclone <- rep(NA, nrow(sample.info.published.data))
names(evidence.for.subclone) <- sample.info.published.data$ID

plotlist.model.vs.data <- list()

for(patient in 1:nrow(sample.info.published.data)){
 
  print(patient)
  
  patient.id <- sample.info.published.data[patient,]$ID
  age <- sample.info.published.data[patient, ]$Age*365
  ncells <- sample.info.published.data[patient,]$ncells
  
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    ###### load observed data
    
    min.vaf <- 0.01
    min.clone.size = 0.01
    min.prior.size = 0.001
   
    directory <- paste0(analysis.directory, "Model_fits/Published_data/", sample.info.published.data[patient,]$Path)
    
    if(!file.exists(paste0(directory, "/Model_fit.csv"))){next}
    fits <- read.csv(paste0(directory, "/Model_fit.csv"))
    
    load(paste0(analysis.directory, "/RData/", sample.info.published.data[patient,]$SNVs_file_name))
    
    snvs.bm <- lapply(snvs, function(x){
      data.frame(VAF=x, Depth=100, varCounts=x*100)
    })
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    #### Print parameter estimates
    
    pdf(paste0(directory, "/Parameter_estimates.pdf"), width=5, height=5)
    
    # the model converts the prior distribution for t_s and s into absolute values to match clones within the limits of min.prior.size and N. We here convert the parameter estimates accordingly.
    fits$par_t_s_absolute <- apply(fits, 1, function(x){
      min.t.s <- 0
      max.t.s <- age - log(min.prior.size*10^as.numeric(x["par_N"]))/10^as.numeric(x["par_lambda_ss"])
      t.s <- min.t.s + as.numeric(x["par_t_s"])*(max.t.s - min.t.s)
    })
    
    fits$par_s_absolute <- apply(fits, 1, function(x){
      ts <- as.numeric(x["par_t_s_absolute"])
      
      min.s <- (10^as.numeric(x["par_lambda_ss"])*(age-as.numeric(ts)) - log(10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*
                                                                                                            (age-as.numeric(ts)))
      if(min.s < 0){
        min.s <- 0
      }
      max.s <- (10^as.numeric(x["par_lambda_ss"])*(age-as.numeric(ts)) - log(min.prior.size*10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*(age-as.numeric(ts)))
      
      s <- min.s + as.numeric(x["par_s"])*(max.s-min.s)
      
    })
    
    fits$size_of_clone <- exp(10^fits$par_lambda_ss*(1-fits$par_s_absolute)*(age - fits$par_t_s_absolute))/10^fits$par_N
    
    fits$size_of_clone[fits$size_of_clone>1] <- 1
    
    fits$age_of_clone <- (age - fits$par_t_s_absolute)/365
    
    ## how long did the clone grow until reaching 2% of the stem cell compartment?
    
    fits$age_of_clone_at_4pct <- log(0.04*10^fits$par_N)/(10^fits$par_lambda_ss*(1-fits$par_s_absolute))/365
    
    fits$growth_per_year <- exp(10^fits$par_lambda_ss*(1-fits$par_s_absolute)*365)-1
    
    fits$mutations_per_year <- fits$par_mu*10^fits$par_lambda_ss*365
    
    to.plot <- fits[,c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "par_t_s_absolute",
                       "par_s_absolute", "size_of_clone", "growth_per_year",
                       "mutations_per_year", "age_of_clone")]
    
    to.plot$par_t_s_absolute <- to.plot$par_t_s_absolute/365
    
    to.plot <- melt(to.plot)
    
    p <- ggplot(to.plot, aes(x=value)) + geom_density(fill="grey") + facet_wrap("variable", scales="free")
    
    print(p)
    
    ## print selected parameters with own prior as x-axis 
    p <- list()
    
    p[[1]] <- ggplot(to.plot[to.plot$variable=="par_N",], aes(x=10^value,  y=0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_log10(limits=10^c(2, 8.5), name="N") + 
      scale_color_manual(values=sample.color)+ scale_fill_manual(values=sample.color)+
      theme(aspect.ratio = 1)
    
    p[[2]] <- ggplot(to.plot[to.plot$variable=="par_mu",], aes(x=value,  y=0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_continuous(limits=c(0, 10), name="Number of SSNVs per division") + 
      scale_color_manual(values=sample.color)+ scale_fill_manual(values=sample.color)+
      theme(aspect.ratio = 1)
    
    p[[3]] <- ggplot(to.plot[to.plot$variable=="par_lambda_ss",], aes(x=value, y = 0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_continuous(limits=c(-3, -1), name="Division rate (1/d)") + 
      scale_color_manual(values=sample.color)+ scale_fill_manual(values=sample.color)+
      theme(aspect.ratio = 1)
    
    p[[4]] <- ggplot(to.plot[to.plot$variable=="size_of_clone",], aes(x=value, y = 0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_continuous(limits=c(0, 1), name="Size of selected clone") + 
      scale_color_manual(values=sample.color)+ scale_fill_manual(values=sample.color)+
      theme(aspect.ratio = 1)
    
    p[[5]] <- ggplot(to.plot[to.plot$variable=="age_of_clone",], aes(x=value,  y=0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_continuous(limits=c(0, age/365), name="Age of selected clone (years)") + 
      scale_color_manual(values=sample.color)+ scale_fill_manual(values=sample.color)+
      theme(aspect.ratio = 1)
    
    p[[6]] <- ggplot(to.plot[to.plot$variable=="growth_per_year",], aes(x=value*100,  y = 0)) +
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_log10( name="Growth of selected clone (% per year)") + 
      scale_color_manual(values=sample.color)+ scale_fill_manual(values=sample.color) +
      theme(aspect.ratio = 1)
    
    
    print(ggarrange(plotlist=p, nrow=2, ncol=3, common.legend = T))
    
    hdinterval <- as.data.frame(t(hdi(fits, credMass=0.8)))
    hdinterval$Parameter <- rownames(hdinterval) 
    hdinterval$Median <- apply(fits, 2, function(x){median(as.numeric(x))})
    hdinterval$Sample <- patient.id

    ### 2D-correlations
    
    ## parameters to plot
    parameter.names <- c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_t_s_absolute",
                         "par_s_absolute", "size_of_clone", "growth_per_year")
    
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
    
    for(i in 1:(sqrt(length(unique(to.plot$par_key))))){
      p[[length(p)+1]] <- ggplot(data.frame()) + geom_point()+
        theme_bw() + theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"),
                           panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
      
    }
    
    for(i in unique(to.plot$par_key)){
      
      
      tmp <- to.plot[to.plot$par_key==i,]
      
      if(tmp$xv[1]=="psurv" | tmp$yv[1]=="psurv"){next}
      
      if(tmp$xv[1]==tmp$yv[1]){
        p[[length(p)+1]] <- ggplot(tmp, aes(x=x)) + 
          geom_histogram() + scale_x_continuous(name=tmp$xv[1]) + scale_y_continuous(name=tmp$yv[1])+
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
          theme(legend.position = "none")
        
      }else{
        tmp[tmp$x==1,]$x <- 1 - abs(rnorm(sum(tmp$x==1), mean = 0, sd = 0.01))
        tmp[tmp$y==1,]$y <- 1 - abs(rnorm(sum(tmp$y==1), mean = 0, sd = 0.01))
        
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
      
      if(tmp$rightm[1]){
        p[[length(p)+1]] <- ggplot(data.frame()) + geom_point()+
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
          theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"),
                legend.position = "none")
        
      }
      
    }
    
    #ggarrange(plotlist=p, nrow=10, ncol=10, align="hv")
    
    
    ## for selection associated values, store only the parameters associated with a clone >= 0.1
    
    fits.selection <- fits[fits$size_of_clone >= 2*0.05, ]
    fits.neutral <- fits[fits$size_of_clone < 2*0.05, ]
    
    if(nrow(fits.selection)>0){
      size.of.selected.clone <- median(fits$size_of_clone[fits$size_of_clone>=0.05])
      evidence.for.subclone[patient.id] <- nrow(fits.selection)/10
    }else{
      evidence.for.subclone[patient.id] <- 0
    }
    
    
    if(nrow(fits.selection)==0){
      
      parameters <- rbind(parameters, hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset",
                                                   "mutations_per_year"),])
      neutral.parameters <- rbind(neutral.parameters,hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset",
                                                                  "mutations_per_year"),])
      
      p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free") 
      
      print(p)
      dev.off()
    }else{
      hdinterval.selection <- as.data.frame(t(hdi(fits.selection, credMass=0.8)))
      hdinterval.selection <- as.data.frame(t(hdi(fits.selection, credMass=0.8)))
      hdinterval.selection$Parameter <- rownames(hdinterval.selection) 
      hdinterval.selection$Median <- apply(fits.selection, 2, function(x){median(as.numeric(x))})
      hdinterval.selection$Sample <- patient.id

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

        neutral.parameters <- rbind(neutral.parameters, hdi.neutral[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset",
                                                                      "mutations_per_year"),])
        
      }
      dev.off()
    }
    
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    ###### Plot fits for neutral and selected case
    
    ## no subsampling due to single-cell info
    depth <- 10000
    use.sensitivity <- F
    
    snvs <- list(snvs.bm[[ sample.info.published.data[patient,]$SAMPLE]])
    
    source("./Parameter_estimation/Bayesian_fit.R")
    
    if( !file.exists(paste0(directory, "/Sim_trajectories.RData"))  ){
      sim.BM <- matrix(0, nrow=100, ncol=length(mySumStatData$mutation.count[[1]]))
      for(j in 1:100){
        print(j)
        parms <- list(mu=fits$par_mu[j], N=fits$par_N[j], delta_exp = fits$par_delta_exp[j], lambda_ss=fits$par_lambda_ss[j],
                      offset=fits$par_offset[j], t_s=fits$par_t_s[j], s=fits$par_s[j])
        
        model <- myModel(parms)
        sim.BM[j,] <- model$modelResult[[1]]
        
      }
      
      max.pred.BM <- apply(sim.BM, 2, quantile, p=0.975)
      min.pred.BM <- apply(sim.BM, 2, quantile, p=0.025)
      
      best.fit <- which.min(fits$distance)
      best.parms <- list(mu=fits$par_mu[best.fit], N=fits$par_N[best.fit], delta_exp = fits$par_delta_exp[best.fit], lambda_ss=fits$par_lambda_ss[best.fit],
                         offset=fits$par_offset[best.fit], t_s=fits$par_t_s[best.fit], s=fits$par_s[best.fit])
      
      best.fit <- myModel(best.parms)
      best.fit <- best.fit$modelResult
      
      data.vs.prediction <- data.frame(VAF=rep(vafs.of.interest), mean.BM=mySumStatData$mutation.count[[1]],
                                       sd.BM = mySumStatData$sampled.sd[[1]],
                                       min.model.BM = min.pred.BM, max.model.BM=max.pred.BM, 
                                       best.fit.model.BM = best.fit[[1]], 
                                       Age=mySumStatData$age/365)
      
      save(sim.BM, data.vs.prediction, file=paste0(directory, "/Sim_trajectories.RData"))
      
    }else{
      load(paste0(directory, "/Sim_trajectories.RData"))
    }
    
    to.plot <- data.vs.prediction
    
    grDevices::pdf(paste0(directory, "/Model_fit.pdf"), width=3, height=2.5)
    
    max.y <-  max(to.plot$max.model.BM, to.plot$mean.BM + to.plot$sd.BM)
    
    p <- ggplot(data=to.plot, aes(x=VAF, y=mean.BM, ymin=mean.BM-sd.BM, ymax=mean.BM+sd.BM)) +
      geom_ribbon(data=to.plot, aes(x=VAF, y=mean.BM, ymin=min.model.BM,ymax=max.model.BM), alpha=1, fill=model.colors["selection"]) + ggtitle("Bone marrow")+
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + scale_y_continuous(name="Cumulative # of mutations") + 
      scale_x_continuous(limits=c(0, 0.6)) + geom_line(aes(x=VAF, y=best.fit.model.BM), inherit.aes = F, col="black")
    
    print(p)
    
    xlimits <- c(0,1/min.clone.size)
    
    p <- ggplot(data=to.plot, aes(x=1/VAF, y=mean.BM, ymin=mean.BM-sd.BM, ymax=mean.BM+sd.BM)) +
      geom_ribbon(data=to.plot, aes(x=1/VAF, y=mean.BM, ymin=min.model.BM,ymax=max.model.BM), alpha=1, fill=model.colors["selection"]) + ggtitle("Bone marrow")+
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + 
      scale_y_continuous(name="Cumulative # of mutations") + 
      geom_line(aes(x=1/VAF, y=best.fit.model.BM), inherit.aes = F, col="black")+ theme(aspect.ratio = 1) +
      scale_x_continuous(limits=xlimits) + coord_cartesian(ylim=c(0, max.y))
    
    ## if the selection model fits the data, illustrate the position of the selected clone
    if(evidence.for.subclone[patient.id]>=15){
      p <- p + geom_ribbon(data= data.frame(x=2/(unlist(hdinterval.selection[hdinterval.selection$Parameter=="size_of_clone",c("lower", "upper")])),
                                            ymin=c(0,0), ymax=c(max.y, max.y)),
                           aes(x = x, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = F)
    }
    
    print(p)
    plotlist.model.vs.data[[patient.id]] <- p
    
    dev.off()
    
  

}

save(parameters, neutral.parameters, selected.parameters, evidence.for.subclone, plotlist.model.vs.data, file= paste0(rdata.directory, "Cohort_parameters_published_data.RData"))

############################################################################################################################################
## Classify samples

## normal samples (no driver, no evidence for selection)
normal.samples.published.data <- sample.info.published.data[sample.info.published.data$CHIP.mutation %in% c("healthy donor", "unknown driver"),]$ID
normal.samples.published.data <- intersect(normal.samples.published.data, names(evidence.for.subclone)[evidence.for.subclone<15]) ## require < 15 for clear-cut normals

## selected samples (driver and evidence for selection)
selected.samples.published.data <- sample.info.published.data[!sample.info.published.data$CHIP.mutation %in% c("healthy donor", "unknown driver", "multiple drivers"),]$ID
selected.samples.published.data <- intersect(selected.samples.published.data, names(evidence.for.subclone)[evidence.for.subclone>=15]) ## require < 15 for clear-cut normals

## selected samples (no driver, but evidence for selection)
selected.no.driver.published.data <- sample.info.published.data[sample.info.published.data$CHIP.mutation %in% c("healthy donor", "unknown driver", "multiple drivers"),]$ID
selected.no.driver.published.data <- intersect(selected.no.driver.published.data, names(evidence.for.subclone)[evidence.for.subclone>=15]) ## require < 15 for clear-cut normals

############################################################################################################################################
## Fig. 3c, e, m, Plot evidence for selection per sample

to.plot <- data.frame(Patient = names(evidence.for.subclone),
                      P_selection = evidence.for.subclone)
to.plot$P_neutral <- 100 - to.plot$P_selection
to.plot <- melt(to.plot, value.name = "Posterior probability")

to.plot$Clone_size <- apply(to.plot, 1, function(x){
  res <- selected.parameters[selected.parameters$Sample==x[1] & selected.parameters$Parameter=="size_of_clone" 
                             & x[2]=="P_selection",]$Median
  if(length(res)==0){
    return(0)}else{
      return(res)
    }
})
to.plot$variable <- factor(to.plot$variable, levels=c("P_neutral", "P_selection"))
to.plot$Clone_size[to.plot$Clone_size==0] <- NA

pdf("./Figures/Figure_3c_e_m.pdf", width=6, height=6)

## normals only
to.plot.normals <- to.plot[to.plot$Patient %in% normal.samples.published.data,]
to.plot.normals$Patient <- factor(to.plot.normals$Patient, levels=c("Lee_Six_Caveman", "Lee_Six_Mutect_Strelka", "AX001", "KX001", "KX002"))

ggplot(to.plot.normals, aes(x=Patient, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
   scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  theme( strip.background = element_blank() )

## selected only

to.plot.selection <- to.plot[to.plot$Patient %in% selected.samples.published.data,]

ggplot(to.plot.selection, aes(x=Patient, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

## unknown drivers

to.plot.normals.w.selection <- to.plot[to.plot$Patient %in% selected.no.driver.published.data,]
to.plot.normals.w.selection$Patient <- factor(to.plot.normals.w.selection$Patient, 
                                              levels=unique(to.plot.normals.w.selection$Patient)[order(to.plot.normals.w.selection$`Posterior probability`[to.plot.normals.w.selection$variable=="P_selection"],
                                                                                                       decreasing = T)])

ggplot(to.plot.normals.w.selection, aes(x=Patient, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

dev.off()


############################################################################################################################################
## Figure 3d: plot parameters across cohort: normal samples -- no driver, no evidence for selection!

load("RData/Cohort_parameters_published_data.RData")

pdf("Figures/Figure_3d.pdf", width=3.5, height=2)

## Stem cell number
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_N",], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% normal.samples.published.data,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
# for comparison, add the estimate from Lee-Six et al.
to.plot <- rbind(to.plot, data.frame(Sample="Orig. Study", lower = log10(50000), upper = log10(200000),
                                     Parameter="par_N", variable = "Median", Median = log10(100000)))
## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% normal.samples.published.data &
                                                                        neutral.parameters$Parameter=="par_N",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% normal.samples.published.data &
                                                                        neutral.parameters$Parameter=="par_N",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=unique(c("Orig. Study", "Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six"))))

ggplot(to.plot, 
       aes(x=Sample, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Stem cell number (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  



## lambda_ss
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss",], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% normal.samples.published.data,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
# for comparison, add the estimate from Lee-Six et al.
to.plot <- rbind(to.plot, data.frame(Sample="Orig. Study", lower = log10(1/(20*30)), upper = log10(1/(2*30)),
                                     Parameter="par_N", variable = "Median", Median = log10(1/150)))
## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% normal.samples.published.data  &
                                                                        neutral.parameters$Parameter=="par_lambda_ss",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% normal.samples.published.data &
                                                                        neutral.parameters$Parameter=="par_lambda_ss",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=unique(c("Orig. Study", "Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six"))))

ggplot(to.plot, 
       aes(x=Sample, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = 365*10^min, ymax = 365*10^max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)   


## Mutation rate
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_mu" ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% normal.samples.published.data ,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% normal.samples.published.data  &
                                                                        neutral.parameters$Parameter=="par_mu",]$lower, p = 0.025), 2),
                                max = rep(quantile(neutral.parameters[neutral.parameters$Sample %in% normal.samples.published.data  &
                                                                        neutral.parameters$Parameter=="par_mu",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=c("Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six")))

ggplot(to.plot, 
       aes(x=Sample, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  

dev.off()


####################################################################################################################################################
## Fig. 3h, predicted vs measured age for samples id2259 (Fabre et al.) and KX004 (Mitchell et al.)

fits <- read.csv(paste0(analysis.directory, "./Model_fits/Published_data/Fabre_et_al/id2259/Model_fit.csv"))

age <- sample.info.published.data[sample.info.published.data$SAMPLE=="id2259",]$Age*365

# the model converts the prior distribution for t_s and s into absolute values to match clones within the limits of min.prior.size and N. We here convert the parameter estimates accordingly.
fits$par_t_s_absolute <- apply(fits, 1, function(x){
  min.t.s <- 0
  max.t.s <- age - log(0.001*10^as.numeric(x["par_N"]))/10^as.numeric(x["par_lambda_ss"])
  t.s <- min.t.s + as.numeric(x["par_t_s"])*(max.t.s - min.t.s)
})

fits$par_s_absolute <- apply(fits, 1, function(x){
  ts <- as.numeric(x["par_t_s_absolute"])
  
  min.s <- (10^as.numeric(x["par_lambda_ss"])*(age-as.numeric(ts)) - log(10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*
                                                                                                        (age-as.numeric(ts)))
  if(min.s < 0){
    min.s <- 0
  }
  max.s <- (10^as.numeric(x["par_lambda_ss"])*(age-as.numeric(ts)) - log(0.001*10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*(age-as.numeric(ts)))
  
  s <- min.s + as.numeric(x["par_s"])*(max.s-min.s)
  
})

fits$growth_per_year <- exp(10^fits$par_lambda_ss*(1-fits$par_s_absolute)*365)-1

fits$age_of_clone <- (age - fits$par_t_s_absolute)/365


## compare our estimate of the clone age to the estimates in Fabre et al. and Mitchell et al. (manually extracted from their trees):

to.plot <- data.frame(ID = c(rep("id2259", 2), rep("KX004", 2)),
                      Age_of_clone_median = c(selected.parameters[selected.parameters$Parameter=="age_of_clone" &
                                                                    selected.parameters$Sample=="id2259",]$Median, 27.8,
                                              selected.parameters[selected.parameters$Parameter=="age_of_clone" &
                                                                    selected.parameters$Sample=="KX004",]$Median, 46),
                      Age_of_clone_min = c(selected.parameters[selected.parameters$Parameter=="age_of_clone" &
                                                                 selected.parameters$Sample=="id2259",]$lower, NA, NA,
                                           selected.parameters[selected.parameters$Parameter=="age_of_clone" &
                                                                 selected.parameters$Sample=="KX004",]$lower, NA),
                      Age_of_clone_max = c(selected.parameters[selected.parameters$Parameter=="age_of_clone" &
                                                                 selected.parameters$Sample=="id2259",]$upper, NA, NA,
                                           selected.parameters[selected.parameters$Parameter=="age_of_clone" &
                                                                 selected.parameters$Sample=="KX004",]$upper, NA),
                      Method = c("FLORENCE", "Phylodynamics", "FLORENCE", "Phylodynamics"))

to.plot$Driver <- sample.info.published.data[sapply(to.plot$ID, function(x){which(sample.info.published.data$SAMPLE==x)}),]$CHIP.mutation

pdf("./Figures/Figure_3h.pdf", width=4, height = 4)

## 2D-scatter plot:

to.plot.2 <- reshape(to.plot, v.names = c("Age_of_clone_min",
                                          "Age_of_clone_median", "Age_of_clone_max"), idvar = "ID", timevar = "Method", direction = "wide") 

ggplot(to.plot.2[to.plot.2$ID %in% selected.samples.published.data,], aes(y=Age_of_clone_median.FLORENCE, ymin=Age_of_clone_min.FLORENCE, ymax=Age_of_clone_max.FLORENCE,
                      x = Age_of_clone_median.Phylodynamics, col=Driver)) + geom_point() + geom_pointrange() + expand_limits(x = 0, y = 0) +
  geom_abline(slope=1, intercept = 0, linetype=2) + ggtitle("Selection with known driver") +
  scale_color_manual(values=CHIP.color) + scale_x_continuous(breaks = seq(0, 50, 25), 
                                                             labels=c("0", "25", "50"), limits=c(0, 50)) +
  scale_y_continuous(breaks= c(0, 25, 50), labels=c("0", "25", "50")) + theme(aspect.ratio = 1)

dev.off()


############################################################################################################################################
## Fig. 3i, j: Plot parameters across cohort: selected samples -- driver and evidence for selection

load("./RData/Cohort_parameters_published_data.RData")

pdf("Figures/Figure_3_i_j.pdf", width=3.5, height=2)

## Stem cell number
to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_N",], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selected.samples.published.data,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
to.plot$Driver <- sample.info.published.data[sapply(to.plot$Sample, function(x){which(sample.info.published.data$SAMPLE==x)}),]$CHIP.mutation

## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.samples.published.data &
                                                                         selected.parameters$Parameter=="par_N",]$lower, p = 0.025), 2),
                                max = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.samples.published.data & 
                                                                         selected.parameters$Parameter=="par_N",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=c("Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six")))

ggplot(to.plot, 
       aes(x=Sample, y=Median, ymin=lower, ymax=upper, col=Driver)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Stem cell number (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  



## lambda_ss
to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_lambda_ss",], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selected.samples.published.data,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
to.plot$Driver <- sample.info.published.data[sapply(to.plot$Sample, function(x){which(sample.info.published.data$SAMPLE==x)}),]$CHIP.mutation
## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.samples.published.data & 
                                                                         selected.parameters$Parameter=="par_lambda_ss",]$lower, p = 0.025), 2),
                                max = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.samples.published.data & 
                                                                         selected.parameters$Parameter=="par_lambda_ss",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=c("Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six")))

ggplot(to.plot, 
       aes(x=Sample, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col=Driver)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = 365*10^min, ymax = 365*10^max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)   


## Mutation rate
to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_mu",], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selected.samples.published.data,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
to.plot$Driver <- sample.info.published.data[sapply(to.plot$Sample, function(x){which(sample.info.published.data$SAMPLE==x)}),]$CHIP.mutation
## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.samples.published.data & 
                                                                selected.parameters$Parameter=="par_mu",]$lower, p = 0.025), 2),
                                max = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.samples.published.data & 
                                                                selected.parameters$Parameter=="par_mu",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=c("Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six")))

ggplot(to.plot, 
       aes(x=Sample, y=Median, ymin=lower, ymax=upper, col=Driver)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  


## Age of the selected clone

to.plot <- selected.parameters[selected.parameters$Parameter=="age_of_clone" & selected.parameters$Sample %in% selected.samples.published.data,]
to.plot$CHIP.mutation <- sapply(to.plot$Sample, function(x){
  sample.info.published.data[sample.info.published.data$ID==x,]$CHIP.mutation})

ggplot(to.plot, 
       aes(x=CHIP.mutation, y=Median, ymin=lower, ymax=upper, col=CHIP.mutation)) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of clone") +
  expand_limits(x = 0, y = 0) + theme(aspect.ratio = 1)


## Growth rate of the selected clone

to.plot <- selected.parameters[selected.parameters$Parameter=="growth_per_year" & selected.parameters$Sample %in% selected.samples.published.data,]
to.plot$CHIP.mutation <- sapply(to.plot$Sample, function(x){
  sample.info.published.data[sample.info.published.data$ID==x,]$CHIP.mutation})

ggplot(to.plot, 
       aes(x=CHIP.mutation, y=Median*100, ymin=lower*100, ymax=upper*100, col=CHIP.mutation)) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="% Growth per year") +
  expand_limits(x = 0, y = 0) + theme(aspect.ratio = 1)


dev.off()


############################################################################################################################################
## Fig. 3n,o Plot parameters across cohort: selected sample with unknown driver


load( "RData/Cohort_parameters_published_data.RData")

pdf("Figure/Figure_3n_o.pdf", width=3.5, height=2)

## Stem cell number
to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_N"  ,], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selected.no.driver.published.data,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.no.driver.published.data & 
                                                                         selected.parameters$Parameter=="par_N",]$lower, p = 0.025), 2),
                                max = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.no.driver.published.data & 
                                                                        selected.parameters$Parameter=="par_N",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=c("Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six")))

ggplot(to.plot, 
       aes(x=Sample, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Stem cell number (log10)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  



## lambda_ss
to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_lambda_ss",], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selected.no.driver.published.data,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.no.driver.published.data & 
                                                                         selected.parameters$Parameter=="par_lambda_ss",]$lower, p = 0.025), 2),
                                max = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.no.driver.published.data & 
                                                                        selected.parameters$Parameter=="par_lambda_ss",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=c("Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six")))

ggplot(to.plot, 
       aes(x=Sample, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = 365*10^min, ymax = 365*10^max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Division rate (1/y)") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)   


## Mutation rate
to.plot <- melt(selected.parameters[selected.parameters$Parameter=="par_mu",], 
                id.vars = c("Sample", "lower", "upper", "Parameter"), value.name = "Median")
to.plot <- to.plot[to.plot$Sample %in% selected.no.driver.published.data,]
to.plot$Sample <- replace(to.plot$Sample, to.plot$Sample=="Lee_Six_Mutect_Strelka", "Lee_Six")
## take medians of lower and upper quantiles as lower and upper bounds across the cohort
to.plot.quantiles <- data.frame(min = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.no.driver.published.data & 
                                                                         selected.parameters$Parameter=="par_mu",]$lower, p = 0.025), 2),
                                max = rep(quantile(selected.parameters[selected.parameters$Sample %in% selected.no.driver.published.data & 
                                                                        selected.parameters$Parameter=="par_mu",]$upper, p = 0.975), 2),
                                x=c(0, nrow(to.plot)))

to.plot$Sample <- factor(to.plot$Sample, levels=c("Lee_Six", setdiff(unique(to.plot$Sample),"Lee_Six")))

ggplot(to.plot, 
       aes(x=Sample, y=Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  geom_ribbon(data=to.plot.quantiles, aes(ymin = min, ymax = max, x=x), inherit.aes = F, fill="grey", alpha = 0.7) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  theme(aspect.ratio = 1)+ expand_limits(x = 0, y = 0)  



## Age of the selected clone

to.plot <- selected.parameters[selected.parameters$Parameter=="age_of_clone" & selected.parameters$Sample %in% selected.no.driver.published.data,]
to.plot$CHIP.mutation <- sapply(to.plot$Sample, function(x){
  sample.info.published.data[sample.info.published.data$ID==x,]$CHIP.mutation})
to.plot$CHIP.mutation <- replace(to.plot$CHIP.mutation, to.plot$CHIP.mutation=="healthy donor", "unknown driver")

ggplot(to.plot, 
       aes(x=Sample, y=Median, ymin=lower, ymax=upper, col=CHIP.mutation)) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of clone") +
  expand_limits(x = 0, y = 0) + theme(aspect.ratio = 1)


## Growth rate of the selected clone

to.plot <- selected.parameters[selected.parameters$Parameter=="growth_per_year" & selected.parameters$Sample %in% selected.no.driver.published.data,]
to.plot$CHIP.mutation <- sapply(to.plot$Sample, function(x){
  sample.info.published.data[sample.info.published.data$ID==x,]$CHIP.mutation})
to.plot$CHIP.mutation <- replace(to.plot$CHIP.mutation, to.plot$CHIP.mutation=="healthy donor", "unknown driver")

ggplot(to.plot, 
       aes(x=Sample, y=Median*100, ymin=lower*100, ymax=upper*100, col=CHIP.mutation)) +
  geom_pointrange() + scale_color_manual(values=CHIP.color) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="% Growth per year") +
  expand_limits(x = 0, y = 0) + theme(aspect.ratio = 1)

dev.off()

