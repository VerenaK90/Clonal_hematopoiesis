###### This scripts plots the output from the population genetics model for Brain WGS data from Bae et al. It produces
###### - For each patient the model fit vs data and the posterior probabilities of the parameters
###### - A summary of the 80%HDI estimates across the population
############################################################################################################################################
###### libraries and functions
library(cdata)
library(scales)
library(openxlsx)
library(SCIFER)

source("./Settings.R")

load("./RData/Bae_et_al/SNVs_brain.RData")
snvs.brain <- snvs

analyzed.samples <- read.xlsx("MetaData/Supplementary Tables.xlsx", startRow = 7, sheet = 8)
analyzed.samples$Age <- as.numeric(analyzed.samples$Age)
analyzed.samples$Average.coverage <- as.numeric(analyzed.samples$Average.coverage)
analyzed.samples$Sample.ID <- replace(analyzed.samples$Sample.ID, analyzed.samples$Sample.ID == "M3663M ", "M3663M")

############################################################################################################################################
####### Set global parameters

sex.colors <- c(Male = "black", Female = "grey")
phenotype.colors <- pal_jco(palette = c("default"), alpha = 1)(4)
region.colors <- brewer_pal(palette="Set1")(3)
selection.colors <- c(Selection = "firebrick", Neutral = "grey")

############################################################################################################################################
## Summarize the parameters across the analyzed cohort

parameters <- data.frame()
neutral.parameters <- data.frame()
selected.parameters <- data.frame()
model.support.selection <- matrix(NA, nrow=1, ncol=nrow(analyzed.samples),
                                dimnames = list("", analyzed.samples$Sample.ID))

use.sensitivity = F

plotlist.model.vs.data <- list()

for(sample.nr in 1:nrow(analyzed.samples)){
  
  print(sample.nr)
  
  donor <- analyzed.samples[sample.nr,]$Case.ID
  age <- analyzed.samples[sample.nr,]$Age*365
  depth <- analyzed.samples[sample.nr,]$Average.coverage
  sample <- analyzed.samples[sample.nr,]$Sample.ID
  
  # set model parameters according to sequencing depth
  if(depth < 150){
    min.vaf <- 0.05
    min.clone.size = 0.05
    min.prior.size=0.01
    
  }else{
    
    min.vaf <- 0.02
    min.clone.size = 0.01
    min.prior.size=0.001
  }
  

    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    ###### load observed data
    
    directory <- paste0(analysis.directory, "/Model_fits/Published_data/Bae_et_al/", patient.id)
  
    fits <- read.csv(paste0(directory, "/Model_fit.csv"))
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    #### Print parameter estimates
    
    # transform the parameters into absolute values and compute size and age of the selected clone.
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
    
    fits$growth_per_year <- exp(10^fits$par_lambda_ss*(1-fits$par_s_absolute)*365)-1
    
    fits$mutations_per_year <- fits$par_mu*10^fits$par_lambda_ss*365
    
    # N x tau per year
    
    fits$N_tau <- 10^fits$par_N / (10^fits$par_lambda_ss*365)
    
    snvs <- list(snvs.brain[[sample]])
    if(sample %in% names(snvs.brain)){
      snvs <- list(snvs.brain[[sample]])
    }else{
      snvs <- list(snvs.brain[[donor]])
    }
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
    #### Print parameter estimates
    
    pdf(paste0(directory, "/Parameter_estimates.pdf"), width=5, height=5)
    
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
      theme(aspect.ratio = 1)
    
    p[[2]] <- ggplot(to.plot[to.plot$variable=="par_mu",], aes(x=value,  y=0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_continuous(limits=c(0, 10), name="Number of SSNVs per division") + 
      theme(aspect.ratio = 1)
    
    p[[3]] <- ggplot(to.plot[to.plot$variable=="par_lambda_ss",], aes(x=value, y = 0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_continuous(limits=c(-3, -1), name="Division rate (1/d)") + 
      theme(aspect.ratio = 1)
    
    p[[4]] <- ggplot(to.plot[to.plot$variable=="size_of_clone",], aes(x=value, y = 0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_continuous(limits=c(0, 1), name="Size of selected clone") + 
      theme(aspect.ratio = 1)
    
    p[[5]] <- ggplot(to.plot[to.plot$variable=="age_of_clone",], aes(x=value,  y=0)) + 
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_continuous(limits=c(0, age/365), name="Age of selected clone (years)") + 
      theme(aspect.ratio = 1)
    
    p[[6]] <- ggplot(to.plot[to.plot$variable=="growth_per_year",], aes(x=value*100,  y = 0)) +
      stat_density_ridges(quantile_lines = TRUE, quantile_fun = hdi_custWidth, quantile = 0.8, vline_linetype = 2, alpha = 0.7, fill="grey") +
      scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")+
      scale_x_log10( name="Growth of selected clone (% per year)") + 
      theme(aspect.ratio = 1)
    
    
    print(ggarrange(plotlist=p, nrow=2, ncol=3, common.legend = T))
    
    hdinterval <- as.data.frame(t(hdi(fits, credMass=0.8)))
    hdinterval$Parameter <- rownames(hdinterval) 
    hdinterval$Median <- apply(fits, 2, function(x){median(as.numeric(x))})
    hdinterval$Sample <- sample

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
      scale_x_log10(name="mu*N (1/division)", limits=c(0.1*10^2.5, 10*10^8)) + scale_y_continuous(name="delta", limits = c(0, 0.75))+
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
    
    
    ## for selection associated values, store only the parameters associated with a clone >= 0.1
    
    fits.selection <- fits[fits$size_of_clone >= 2*0.05, ]
    fits.neutral <- fits[fits$size_of_clone < 2*0.05, ]
    
    if(nrow(fits.selection)>0){
      size.of.selected.clone <- median(fits$size_of_clone[fits$size_of_clone>=0.05])
      model.support.selection[1,sample] <- nrow(fits.selection)/10
    }else{
      model.support.selection[1,sample] <- 0
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
      hdinterval.selection <- as.data.frame(t(hdi(fits.selection, credMass=0.8)))
      hdinterval.selection <- as.data.frame(t(hdi(fits.selection, credMass=0.8)))
      hdinterval.selection$Parameter <- rownames(hdinterval.selection) 
      hdinterval.selection$Median <- apply(fits.selection, 2, function(x){median(as.numeric(x))})
      hdinterval.selection$Sample <- sample

      p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free") + ggtitle("All parameters")
      
      print(p)
      
      p <- ggplot(data=hdinterval.selection, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free") + ggtitle("Selection")
      
      print(p)
      
      parameters <- rbind(parameters, hdinterval)
      selected.parameters <- rbind(selected.parameters, hdinterval.selection)
      
      if(nrow(fits.neutral)>0){
        hdi.neutral <- as.data.frame(t(hdi(fits.neutral, credMass = 0.8)))
        hdi.neutral$Parameter <- rownames(hdi.neutral) 
        hdi.neutral$Median <- apply(fits.neutral, 2, function(x){median(as.numeric(x))})
        hdi.neutral$Sample <- sample

        p <- ggplot(data=hdi.neutral, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
          facet_wrap(~Parameter, scales="free") + ggtitle("Neutral")
        
        print(p)
        
        neutral.parameters <- rbind(neutral.parameters, hdi.neutral[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset", "N_tau",
                                                                      "mutations_per_year"),])
        
      }
      dev.off()
    }
    
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    ###### Plot fits for neutral and selected case
    
    use.sensitivity <- F
    
    source(paste0(custom.script.directory, "/Parameter_estimation/Bayesian_fit.R"))
    
    if( !file.exists(paste0(directory, "/Sim_trajectories.RData"))  ){
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
    
    max.y <-  max(to.plot$max.model, to.plot$mean + to.plot$sd)
    
    p <- ggplot(data=to.plot, aes(x=VAF, y=mean.data, ymin=mean.data-sd.data, ymax=mean.data+sd.data)) +
      geom_ribbon(data=to.plot, aes(x=VAF, y=mean.data, ymin=min.model,ymax=max.model), alpha=1, fill="purple") + 
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + scale_y_continuous(name="Cumulative # of mutations") + 
      scale_x_continuous(limits=c(0, 0.6)) 
    
    print(p)
    
    xlimits <- c(0,1/min.clone.size)
    
    p <- ggplot(data=to.plot, aes(x=1/VAF, y=mean.data, ymin=mean.data-sd.data, ymax=mean.data+sd.data)) +
      geom_ribbon(data=to.plot, aes(x=1/VAF, y=mean.data, ymin=min.model,ymax=max.model), alpha=1, fill="purple") + 
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + 
      scale_y_continuous(name="Cumulative # of mutations") + theme(aspect.ratio = 1) +
      scale_x_continuous(limits=xlimits) + coord_cartesian(ylim=c(0, max.y))
    
    
    ## if the selection model fits the data, illustrate the position of the selected clone
    if(nrow(fits.selection) > 150){
      p <- p + geom_ribbon(data= data.frame(x=2/(unlist(hdinterval.selection[hdinterval.selection$Parameter=="size_of_clone",c("lower", "upper")])),
                                            ymin=c(0,0), ymax=c(max.y, max.y)),
                           aes(x = x, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = F)
    }
    

    
    if(any(!is.na(snvs[[1]]$Cancer.Driver.Genes) & snvs[[1]]$`Consequence.(VEP)` %in%
           c("missense_variant", "splice_donor_variant", "stop_gained", "TF_binding_site_variant",
             "nonsense_variant", "stop_lost"))){
      
      driver.information <- snvs[[1]][!is.na(snvs[[1]]$Cancer.Driver.Genes) & snvs[[1]]$`Consequence.(VEP)` %in%
                                   c("missense_variant", "splice_donor_variant", "stop_gained", "TF_binding_site_variant",
                                     "nonsense_variant", "stop_lost"),]
      
      driver.information$upper <- driver.information$VAF + 1.96*sqrt(driver.information$VAF*(1-driver.information$VAF)/driver.information$Depth)
      driver.information$lower <- driver.information$VAF - 1.96*sqrt(driver.information$VAF*(1-driver.information$VAF)/driver.information$Depth)
      driver.information$lower[driver.information$lower < min.clone.size] <- min.clone.size
      
      driver.information$y <- sapply(driver.information$VAF, function(x){
        to.plot[min(which(to.plot$VAF >= x)),]$max.model
      })
      
      p <- p + geom_point(data = driver.information[driver.information$upper>=1/max(xlimits),], aes(x = 1/as.numeric(VAF), y = y), col="firebrick", inherit.aes = F) +
        geom_errorbarh(data =  driver.information[driver.information$upper>=1/max(xlimits),], aes(xmin = 1/as.numeric(upper), y = y,
                                                                                                  xmax = 1/as.numeric(lower)), height = max.y/50, col="firebrick", inherit.aes = F) +
        geom_text(data = driver.information[driver.information$upper>=1/max(xlimits),], aes( x=1/as.numeric(VAF), y=y*1.05), label=driver.information[driver.information$upper>=1/max(xlimits),'Gene.Symbol'],  color="firebrick",
                  size=6 ,  fontface="italic", inherit.aes = F)
    }
    if(sample == "LIBD82_Hippocampus"){
      p <- p + geom_point(data = data.frame(VAF = 0.15/2), aes(x=1/VAF, y = max.y/2), inherit.aes = F, col = "firebrick") +
        geom_text(data = data.frame(VAF = 0.15/2),  label = "Trisomy 7/ Monosomy 10", col = "firebrick")
    }
    
    print(p)
    plotlist.model.vs.data[[sample]] <- p
    
    dev.off()
}

parameters$ID <- gsub("-.*", "", parameters$Sample)
neutral.parameters$ID <- gsub("-.*", "", neutral.parameters$Sample)
selected.parameters$ID <- gsub("-.*", "", selected.parameters$Sample)

############################################################################################################################################
## Classify samples

## normal samples (no evidence for selection)
normal.samples.bae.et.al <- colnames(model.support.selection)[model.support.selection[1,]<15] ## require < 15 for clear-cut normals
## selected samples 
selected.samples.bae.et.al <-colnames(model.support.selection)[model.support.selection[1,]>=15] ## require >= 15 for selection

############################################################################################################################################
## Figur 8a-c: plot model fits for 2 positive controls and 2 cases with unknown drivers

pdf(paste0(analysis.directory, "/Figures/Figure_8a-c.pdf"), width=6, height=5)

ggarrange(plotlist=plotlist.model.vs.data[c("LIBD82_Hippocampus",
                                            "NC7-CX-OLI", "NC7-STR-INT",
                                            "LIBD87_DLPFC", "TS1-STR-Bulk")],
  nrow=3, ncol=2)

dev.off()

############################################################################################################################################
## Fig. 8d, plot model support for selection per sample and compare the number of selected cases across phenotypes, sex and region

to.plot <- data.frame(Sample = colnames(model.support.selection),
                      ID = sapply(colnames(model.support.selection), function(x){
                        analyzed.samples[analyzed.samples$Sample.ID==x,]$Case.ID
                      }),
                      P_selection = model.support.selection[1,])
to.plot <- to.plot[!is.na(to.plot$P_selection),]
to.plot$P_neutral <- 100 - to.plot$P_selection
to.plot <- melt(to.plot, value.name = "Posterior probability")

to.plot$Clone_size <- apply(to.plot, 1, function(x){
  res <- selected.parameters[selected.parameters$Sample==x[1] & selected.parameters$Parameter=="size_of_clone" 
                             & x[3]=="P_selection" ,]$Median
  if(length(res)==0){
    return(0)}else{
      return(res)
    }
})
to.plot$variable <- factor(to.plot$variable, levels=c("P_neutral", "P_selection"))
to.plot$Clone_size[to.plot$Clone_size==0] <- NA
# sort Sample by evidence
to.plot$Sample <- factor(to.plot$Sample, levels = unique(to.plot$Sample[order(to.plot$`Posterior probability`[to.plot$variable=="P_selection"])]))


pdf(paste0(analysis.directory, "/Figures/Figure_8d.pdf"), width=10, height=6)

# plot the model support per sample 
p <- ggplot(to.plot[!is.na(to.plot$`Posterior probability`),], aes(x=Sample, y=`Posterior probability`, fill=Clone_size)) + 
  geom_col(width=0.5, col=NA) + scale_x_discrete(name="") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  theme( strip.background = element_blank(), axis.text.x = element_text(angle=90) ) + 
  geom_hline(yintercept = 15, linetype=2) 

print(p)

dev.off()
############################################################################################################################################
## Extended Data Fig. 10b, compare the number of selected cases across phenotypes, sex and region


pdf(paste0(analysis.directory, "/Figures/Figure_S10_b.pdf"), width=4, height=3)

# plot number of individuals with evidence for clonal selection stratified by phenotype  

to.plot <- data.frame(Phenotype = rep(unique(analyzed.samples$Phenotype), 2),
                       Modeltype = rep(c("Neutral", "Selection"), each = length(unique(analyzed.samples$Phenotype))),
                       Number = c(sapply(unique(analyzed.samples$Phenotype), function(x){
                         tmp <- analyzed.samples[analyzed.samples$Phenotype==x,]
                         sum(sapply(unique(tmp$Sample.ID), function(y){
                           !any(tmp[tmp$Sample.ID==y,]$Sample %in% selected.samples.bae.et.al)
                         }))
                       }), sapply(unique(analyzed.samples$Phenotype), function(x){
                         tmp <- analyzed.samples[analyzed.samples$Phenotype==x,]
                         sum(sapply(unique(tmp$Sample.ID), function(y){
                           any(tmp[tmp$Sample.ID==y,]$Sample %in% selected.samples.bae.et.al)
                         }))                       })))
to.plot$Phenotype <- factor(to.plot$Phenotype, levels = c("Normal", "ASD", "SCZ", "TS"))

ggplot(to.plot, aes(x = Phenotype, y = Number, fill = Modeltype)) + geom_col(position = "dodge", col = "black") +
  scale_fill_manual(values = c(Neutral = "white", Selection = "black")) +
  scale_y_continuous('Number of individuals') + ggtitle('Evidence for selection in at least 1 sample')


## stratify by sex 

to.plot <- data.frame(Sex = rep(unique(analyzed.samples$Sex), 2),
                       Modeltype = rep(c("Neutral", "Selection"), each = length(unique(analyzed.samples$Sex))),
                       Number = c(sapply(unique(analyzed.samples$Sex), function(x){
                         tmp <- analyzed.samples[analyzed.samples$Sex==x,]
                         sum(sapply(unique(tmp$Sample.ID), function(y){
                           !any(tmp[tmp$Sample.ID==y,]$Sample %in% selected.samples.bae.et.al)
                         }))
                       }), sapply(unique(analyzed.samples$Sex), function(x){
                         tmp <- analyzed.samples[analyzed.samples$Sex==x,]
                         sum(sapply(unique(tmp$Sample.ID), function(y){
                           any(tmp[tmp$Sample.ID==y,]$Sample %in% selected.samples.bae.et.al)
                         }))                       })))

ggplot(to.plot, aes(x = Sex, y = Number, fill = Modeltype)) + geom_col(position = "dodge", col = "black") +
  scale_fill_manual(values = c(Neutral = "white", Selection = "black")) +
  scale_y_continuous('Number of individuals') + ggtitle('Evidence for selection in at least 1 sample')


# stratify by brain region, count samples rather than indviduals
to.plot <- data.frame(Region = rep(unique(analyzed.samples$Region), 2),
                       Modeltype = rep(c("Neutral", "Selection"), each = length(unique(analyzed.samples$Region))),
                       Number = c(sapply(unique(analyzed.samples$Region), function(x){
                         sum(analyzed.samples$Sample.ID[!analyzed.samples$Region==x] %in% selected.samples.bae.et.al)
                       }), sapply(unique(analyzed.samples$Region), function(x){
                         sum(analyzed.samples$Sample.ID[analyzed.samples$Region==x] %in% selected.samples.bae.et.al)
                       })))

ggplot(to.plot, aes(x = Region, y = Number, fill = Modeltype)) + geom_col(position = "dodge", col = "black") +
  scale_fill_manual(values = c(Neutral = "white", Selection = "black")) +
  scale_y_continuous('Number of samples')

dev.off()


############################################################################################################################################
## Figure 8e, compare incidence of selection across age groups

pdf(paste0(analysis.directory, "/Figures/Figure_8_e.pdf"), width=4, height=3)

# plot incidence of selected cases against age; bin by age group (10 years) 
to.plot <- data.frame(Age = seq(0, 100, 10))
to.plot$Incidence <- sapply(to.plot$Age, function(x){
  tmp <- analyzed.samples[round(analyzed.samples$Age/5)*5 == x,]
  if(nrow(tmp)==0){return(NA)}
  sum(sapply(unique(tmp$Case.ID), function(y){
    any(tmp[tmp$Case.ID==y,]$Sample %in% selected.samples.bae.et.al)
  }))/length(unique(tmp$Case.ID))
})

ggplot(to.plot, aes(x = Age, y = 100*Incidence)) + geom_point() + 
  scale_x_continuous(name = "Age (years)") + scale_y_continuous(name = "Incidence (%)") + 
  geom_smooth(col = "black") + coord_cartesian(ylim = c(0, 100))


# plot selected cases against age

to.plot <- data.frame(Selection = c(),
                       Age = c())

for(i in unique(analyzed.samples$Case.ID)){
  tmp <- analyzed.samples[analyzed.samples$Case.ID==i,]
  if(any(tmp$Sample.ID %in% selected.samples.bae.et.al)){
    to.plot <- rbind(to.plot, data.frame(Selection = T,
                                           Age = unique(analyzed.samples[analyzed.samples$Case.ID==i,]$Age)))
  }else{
    to.plot <- rbind(to.plot, data.frame(Selection = F,
                                           Age = unique(analyzed.samples[analyzed.samples$Case.ID==i,]$Age)))
  }
}

ggplot(to.plot, aes(x = Selection, y = Age)) + geom_boxplot() + geom_beeswarm() + 
  scale_y_continuous('Age (years)') + stat_compare_means()

dev.off()

############################################################################################################################################
## Figure 8f/g, Extended Data Fig. 10c: plot physiological parameters across cohort: number of stem cells, division rate and mutation rate

pdf(paste0(analysis.directory, "/Figures/Figure_8fg_S10c.pdf"), width=5, height=6)

## I. Stem cell number
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_N" & 
                                     neutral.parameters$Sample %in% normal.samples.bae.et.al,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "ID"), value.name = "Median")
to.plot$Selection <- "Neutral"
to.plot.selection <- melt(selected.parameters[selected.parameters$Parameter=="par_N" & 
                                                selected.parameters$Sample %in% selected.samples.bae.et.al,], 
                          id.vars = c("Sample", "lower", "upper", "Parameter", "ID"), value.name = "Median")
to.plot.selection$Selection <- "Selection"

to.plot <- rbind(to.plot, to.plot.selection)
to.plot$Region <- sapply(to.plot$Sample, function(x){
  analyzed.samples[analyzed.samples$Sample.ID==x,]$Region
})
to.plot$Age <- sapply(to.plot$Sample, function(x){
  analyzed.samples[analyzed.samples$Sample.ID==x,]$Age
})

## compare between regions

to.plot. <- melt(to.plot, id.vars = c("Sample", "Region"), measure.vars = c("lower", "upper"))


p1 <- ggplot(to.plot, aes(x = Region, y = Median, ymin = lower, ymax = upper, group = Region)) + 
  geom_boxplot(outliers = T,  fill = NA, linewidth = 1) + scale_y_continuous("Number of stem cells (log10)", limits = c(0, 8)) +
  stat_compare_means(size = 2) +  theme(aspect.ratio = 1) 


## plot across age

p2 <- ggplot(to.plot, aes(x = Age, y = Median, ymin = lower, ymax = upper)) + geom_pointrange(col = "darkgrey", alpha = 0.6) +
  scale_x_continuous(name = "Age (years)") + 
  scale_y_continuous(limits = c(0, 8), name = "Stem cell number (log10)") +
  geom_smooth()



## II. Division rate 
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_lambda_ss"& 
                                     neutral.parameters$Sample %in% normal.samples.bae.et.al,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "ID"), value.name = "Median")
to.plot$Selection <- "Neutral"
to.plot.selection <- melt(selected.parameters[selected.parameters$Parameter=="par_lambda_ss" & 
                                                selected.parameters$Sample %in% selected.samples.bae.et.al,], 
                          id.vars = c("Sample", "lower", "upper", "Parameter", "ID"), value.name = "Median")
to.plot.selection$Selection <- "Selection"
to.plot <- rbind(to.plot, to.plot.selection)
to.plot$Region <- sapply(to.plot$Sample, function(x){
  analyzed.samples[analyzed.samples$Sample.ID==x,]$Region
})

to.plot$Age <- sapply(to.plot$Sample, function(x){
  analyzed.samples[analyzed.samples$Sample.ID==x,]$Age
})


## compare between regions

to.plot. <- melt(to.plot, id.vars = c("Sample", "Region"), measure.vars = c("lower", "Median", "upper"))


p3 <- ggplot(to.plot, aes(x = Region, y = 10^Median*365, ymin = 10^lower*365, ymax = 10^upper*365, group = Region)) + 
  geom_boxplot(outliers = T, fill = NA, linewidth = 1) + scale_y_log10("Stem cell divisions per year") +
  stat_compare_means(size = 2) +  theme(aspect.ratio = 1) 

## plot across age

p4 <- ggplot(to.plot, aes(x = Age, y = 10^Median*365, ymin = 10^lower*365, ymax = 10^upper*365)) + geom_pointrange(col = "darkgrey", alpha = 0.6) +
  scale_x_continuous(name = "Age (years)") + 
  scale_y_log10(name = "Division rate (1/y)") +
  geom_smooth()


## III. Mutation rate
to.plot <- melt(neutral.parameters[neutral.parameters$Parameter=="par_mu" &
                                     neutral.parameters$Sample %in% normal.samples.bae.et.al,], 
                id.vars = c("Sample", "lower", "upper", "Parameter", "ID"), value.name = "Median")
to.plot$Selection <- "Neutral"
to.plot.selection <- melt(selected.parameters[selected.parameters$Parameter=="par_mu" & 
                                                selected.parameters$Sample %in% selected.samples.bae.et.al,], 
                          id.vars = c("Sample", "lower", "upper", "Parameter", "ID"), value.name = "Median")
to.plot.selection$Selection <- "Selection"
to.plot <- rbind(to.plot, to.plot.selection)
to.plot$Region <- sapply(to.plot$Sample, function(x){
  analyzed.samples[analyzed.samples$Sample.ID==x,]$Region
})

to.plot$Age <- sapply(to.plot$Sample, function(x){
  analyzed.samples[analyzed.samples$Sample.ID==x,]$Age
})


## compare between regions

to.plot. <- melt(to.plot, id.vars = c("Sample", "Region"), measure.vars = c("lower", "upper"))

p5 <- ggplot(to.plot, aes(x = Region, y = Median, ymin = lower, ymax = upper, group = Region)) +  
  geom_boxplot(outliers = T,  fill = NA, linewidth = 1) + scale_y_continuous("Number of SSNVs per division") +
  stat_compare_means(size = 2) +  theme(aspect.ratio = 1) 

## plot across age

p6 <- ggplot(to.plot, aes(x = Age, y = Median, ymin = lower, ymax = upper)) + geom_pointrange(col = "darkgrey", alpha = 0.6) +
  scale_x_continuous(name = "Age (years)") + 
  scale_y_continuous( name = "Number of SSNVs per division") +
  geom_smooth()

ggpubr::ggarrange(p1, p2, p3, p4, p5, p6 nrow=2)

dev.off()

############################################################################################################################################
## Figure 8h, Extended Data Fig. 10d, plot selection parameters across cohort

pdf(paste0(analysis.directory, "/Figures/Figure_S10d.pdf"), width=5, height=4)

## Age of the selected clone

to.plot <- selected.parameters[selected.parameters$Parameter=="age_of_clone" & selected.parameters$Sample %in% selected.samples.bae.et.al,]
to.plot$Age <- sapply(to.plot$Sample, function(x){
  analyzed.samples[analyzed.samples$Sample.ID==x,]$Age
})

## plot across age

p1 <- ggplot(to.plot, aes(x = Age, y = Median, ymin = lower, ymax = upper)) + geom_pointrange(col="darkgrey", alpha = 0.6) +
  scale_x_continuous(name = "Age (years)") + 
  scale_y_continuous( name = "Age of clone") +
  geom_smooth()


## compute age at driver acquisition
to.plot$Age_at_acquisition <- to.plot$Age - to.plot$Median
to.plot$Age_at_acquisition_lower <- to.plot$Age - to.plot$upper
to.plot$Age_at_acquisition_upper <- to.plot$Age - to.plot$lower

p2 <- ggplot(to.plot, aes(x = Age_at_acquisition, xmin = Age_at_acquisition_lower, xmax = Age_at_acquisition_upper)) +
  stat_ecdf() + stat_ecdf(aes(x=Age_at_acquisition_lower), linetype = 2) +
  stat_ecdf(aes(x=Age_at_acquisition_upper), linetype = 2) +
  scale_x_continuous(name = "Age at acquisition") + scale_y_continuous(name = "Fraction of samples")

## Growth rate of the selected clone

to.plot <- selected.parameters[selected.parameters$Parameter=="growth_per_year" & selected.parameters$Sample %in% selected.samples.bae.et.al,]
to.plot$Age <- sapply(to.plot$Sample, function(x){
  analyzed.samples[analyzed.samples$Sample.ID==x,]$Age
})
## plot across age

p3 <- ggplot(to.plot, aes(x = Age, y = 100*Median, ymin = 100*lower, ymax = 100*upper)) + geom_pointrange(col = "darkgrey", alpha = 0.6) +
  scale_x_continuous(name = "Age (years)") + 
  scale_y_log10( name = "% Growth per year") +
  geom_smooth()

ggpubr::ggarrange(p1, p2, p3, nrow=1, common.legend = T)


dev.off()



