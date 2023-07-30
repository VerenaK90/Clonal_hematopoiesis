source("./Settings.R")
######################### ######################### ######################### ######################### ######################### 
library(cdata)

# define where output is to be stored
study.directory <- paste0("Model_fits/Simulated_data/")
sample.info <- read.xlsx("MetaData/Sample_information_simulated_data.xlsx", sheet = 1)
rownames(sample.info) <- sample.info$SampleID

######################### ######################### ######################### ######################### ######################### 
## Analyze model fits

# collect highest density intervals of the parameters:
# .. overall ..
parameters <- data.frame()
# .. neutral fist only ..
neutral.parameters <- data.frame()
# .. selection fits only
selected.parameters <- data.frame()

# collect model support for selection (%)
model.support.selection <- matrix(NA, nrow=3, ncol=nrow(sample.info), dimnames = list(c("30", "90", "270"),
                                                                                      sample.info$SampleID))
## collect statistics of the fits - how many fits suggest selection, how many suggest neutral, what's the size of the true clone and what's the size of the inferred clone?
inference.stats <- as.data.frame(matrix(0, ncol=4, nrow = nrow(sample.info), 
                            dimnames = list(sample.info$SampleID, c("N_fits_sel", "N_fits_neutral", "Size_true_clone", "Inferred_size"))))
inference.stats$Size_true_clone <- sample.info$CH_VAF
## separate statistics for each simulated sequencing depth
inference.stats <- list(`30`=inference.stats, `90`=inference.stats, `270`=inference.stats)

use.sensitivity <- F
plotlist.model.vs.data <- list()

load(paste0(rdata.directory, "Simulated_data/SNVs_90x.RData"))
snvs.90 <- snvs
load(paste0(rdata.directory, "Simulated_data/SNVs_30x.RData"))
snvs.30 <- snvs
load(paste0(rdata.directory, "Simulated_data/SNVs_270x.RData"))
snvs.270 <- snvs

seq.type <- "bulk"

for(patient.id in sample.info$SampleID){
  print(patient.id)
  age <- sample.info[patient.id, ]$Age*365
  
  depths <- c(30, 90, 270)
  
  for(depth in depths){
    
    ## the parameters are specified in agreement with the fits (greater sensitivity with higher seq. depth)
    if(depth == 90){
      snvs <- list(snvs.90[[patient.id]])
      min.vaf <- 0.05 # minimal VAF analyzed
      min.clone.size = 0.05 # minimal size of clone evaluated with selection model
      min.prior.size = 0.01 # minimal clone size prior
    }else if(depth ==30){
      snvs <- list(snvs.30[[patient.id]])
      min.vaf <- 0.05
      min.clone.size = 0.05
      min.prior.size = 0.01
    }else if(depth==270){ # for 270x coverage, the lower bounds can be decreased
      snvs <- list(snvs.270[[patient.id]])
      min.vaf <- 0.01
      min.clone.size = 0.01
      min.prior.size = 0.001
    }
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    ###### load observed data
    
    directory <- paste0(study.directory, patient.id, "/", depth, "/")
    if(!file.exists(paste0(directory, "/Model_fit.csv"))){
      warning("No file for patient", patient.id)
      next
    }
    fits <- read.csv(paste0(directory, "/Model_fit.csv"))

    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    ###### Plot fits
    
    source("Parameter_estimation/Bayesian_fit.R")
    
    if( !file.exists(paste0(directory, "/Sim_trajectories.RData") )){
      sim <- matrix(0, nrow=100, ncol=length(mySumStatData$mutation.count[[1]]))
      
      # evaluate 100 instances of the model
      for(j in 1:100){
        print(j)
        
        parms <- list(mu=fits$par_mu[j], N=fits$par_N[j], delta_exp = fits$par_delta_exp[j], lambda_ss=fits$par_lambda_ss[j],
                      offset=fits$par_offset[j], t_s=fits$par_t_s[j], s=fits$par_s[j])
        
        model <- myModel(parms)
        sim[j,] <- model$modelResult[[1]]
        
      }
      
      # determine 95% prediction intervals
      max.pred <- apply(sim, 2, quantile, p=0.975)
      min.pred <- apply(sim, 2, quantile, p=0.025)
      
      
      data.vs.prediction <- data.frame(VAF=rep(vafs.of.interest), mean=mySumStatData$mutation.count[[1]],
                                       sd = mySumStatData$sampled.sd[[1]],
                                       min.model = min.pred, max.model=max.pred, 
                                       Age=mySumStatData$age/365)
      
      save(sim, data.vs.prediction, file=paste0(directory, "/Sim_trajectories.RData"))
      
    }else{
      load(paste0(directory, "/Sim_trajectories.RData"))
    }
    
    # plot model fits
    to.plot <- data.vs.prediction
    
    grDevices::pdf(paste0(directory, "/Model_fit.pdf"), width=3, height=2.5)
    
    max.y <- max(to.plot$max.model)
    
    p <- ggplot(data=to.plot, aes(x=VAF, y=mean, ymin=mean-sd, ymax=mean+sd)) +
      geom_ribbon(data=to.plot, aes(x=VAF, y=mean, ymin=min.model,ymax=max.model), alpha=1, fill=model.colors["selection"]) + ggtitle("Bone marrow")+
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + scale_y_continuous(name="Cumulative # of mutations") + 
      scale_x_continuous(limits=c(0, 0.6)) + geom_line(aes(x=VAF, y=best.fit.model.BM), inherit.aes = F, col="black")
    
    if(sample.info[patient.id,]$CH_VAF>0){
      p <- p + geom_vline(xintercept = sample.info[patient.id,]$CH_VAF)
    }
    
    print(p)
    
    if(depth==270){
      xlimits=c(0,100)
      xbreaks=c(100, 50, 25)
      xlabels=c("0.01", "0.02", "0.04")
    }else{
      xlimits <- c(0,20)
      xbreaks=c(20, 10, 5)
      xlabels=c("0.05", "0.1", "0.2")
    }
    
    # plot against 1/VAF
    p <- ggplot(data=to.plot, aes(x=1/VAF, y=mean, ymin=mean-sd, ymax=mean+sd)) +
      geom_ribbon(data=to.plot, aes(x=1/VAF, y=mean, ymin=min.model,ymax=max.model), alpha=1, fill=model.colors["selection"]) + ggtitle("Bone marrow")+
      geom_pointrange(lwd=0.25, shape=1, fatten=1) + 
      scale_y_continuous(name="Cumulative # of mutations") + 
      theme(aspect.ratio = 1) +
      scale_x_continuous(limits=xlimits, breaks=xbreaks, labels=xlabels) + coord_cartesian(ylim=c(0, max.y))
    
    if(sample.info[patient.id,]$CH_VAF>0){
      p <- p + geom_vline(xintercept = 1/sample.info[patient.id,]$CH_VAF)
    }
    
    print(p)
    plotlist.model.vs.data[[patient.id]] <- p
    
    dev.off()
    
    
    ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
    #### Print parameter estimates
    
    pdf(paste0(directory,  "/Parameter_estimates.pdf"), width=5, height=5)
    
    # the model converts the prior distribution for t_s and s into absolute values to match clones within the limits of min.prior.size and N. We here convert the parameter estimates accordingly.
    fits$par_t_s_absolute <- apply(fits, 1, function(x){
      min.t.s <- 0
      max.t.s <- age - log(max(1, min.prior.size*10^as.numeric(x["par_N"])))/10^as.numeric(x["par_lambda_ss"])
      t.s <- min.t.s +as.numeric(x["par_t_s"])*(max.t.s - min.t.s)
    })
    
    fits$par_s_absolute <- apply(fits, 1, function(x){
      ts <- as.numeric(x["par_t_s_absolute"])
      
      min.s <- (10^as.numeric(x["par_lambda_ss"])*(age-ts) - log(10^as.numeric(x["par_N"])))/(10^as.numeric(x["par_lambda_ss"])*
                                                                                                (age-ts))
      if(min.s < 0){
        min.s <- 0
      }
      max.s <- (10^as.numeric(x["par_lambda_ss"])*(age-ts) - log(max(1, min.prior.size*10^as.numeric(x["par_N"]))))/(10^as.numeric(x["par_lambda_ss"])*(age-ts))
      
      s <- min.s + as.numeric(x["par_s"])*(max.s-min.s)
      
    })
    
    # Compute the clone size according to exponential growth
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
    to.plot$par_lambda_ss <- 10^to.plot$par_lambda_ss*365
    
    to.plot <- melt(to.plot)
    
    p <- ggplot(to.plot, aes(x=value, y= 0, fill = stat(quantile))) + 
      geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
      #geom_density(fill="grey") + 
      facet_wrap("variable", scales="free") +
      scale_fill_manual(values = c("transparent", "grey", "transparent"), guide = "none")+
    scale_y_continuous(labels=NULL, breaks=NULL, name="Posterior probability")
    
    print(p)
    
    p <- ggplot(fits, aes(x=par_N, y=par_lambda_ss)) + geom_point() + geom_density2d() + 
      scale_x_continuous(name="Number of stem cells (log10)") + scale_y_continuous(name="Division rate (log10)")
    
    print(p)
    
    # compute highest density intervals for the parameters, 80% credibility mass
    hdinterval <- as.data.frame(t(hdi(fits, credMass=0.8)))
    hdinterval$Parameter <- rownames(hdinterval) 
    hdinterval$Median <- apply(fits, 2, function(x){median(as.numeric(x))})
    hdinterval$Sample <- patient.id
    hdinterval$depth <- depth
    
    ### 2D-correlations
    
    ## parameters to plot
    parameter.names <- c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", 
                         "par_t_s_absolute", "par_s_absolute", "size_of_clone", 
                         "growth_per_year")
    
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
    
    ggarrange(plotlist=p, nrow=10, ncol=10, align="hv")
    
    ## for selection associated values, store only the parameters associated with a clone >= the minimal clone size considered by the model
    
    fits.selection <- fits[fits$size_of_clone >= min.clone.size, ]
    fits.neutral <- fits[fits$size_of_clone < min.clone.size, ]
    
    # store the % of fits that support selection
    model.support.selection[as.character(depth),patient.id] <- nrow(fits.selection)/10
    
    if(any(fits$size_of_clone >= min.clone.size)){
      # compute the estimated clone size and store it
      size.of.selected.clone <- median(fits$size_of_clone[fits$size_of_clone>=min.clone.size])
      # store the number of fits supporting the selected clone
      inference.stats[[as.character(depth)]][patient.id, "N_fits_sel"] <-  sum(fits$size_of_clone>= min.clone.size)
      
      ## round the clone size to 5% if it's at least 5% and in 1% steps below
      if(min.clone.size >= 0.05){
        rounded.clone.size <- round(size.of.selected.clone/0.05)*5
      }else{
        rounded.clone.size <- ifelse(size.of.selected.clone >= 0.05,  round(size.of.selected.clone/0.05)*5,
                                     round(size.of.selected.clone/0.01))
      }
 
      inference.stats[[as.character(depth)]][patient.id, "Inferred_size"] <-   rounded.clone.size
      
    }
    # store the number of fits supporting neutral evolution
    inference.stats[[as.character(depth)]][patient.id, "N_fits_neutral"] <-  1000 - sum(fits$size_of_clone>= min.clone.size)
    
    
    if(nrow(fits.selection)<=1){
      
      parameters <- rbind(parameters, hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset",
                                                   "mutations_per_year"),])
      
      neutral.parameters <- rbind(neutral.parameters,hdinterval[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset",
                                                                  "mutations_per_year"),])
      
      p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
        facet_wrap(~Parameter, scales="free") 
      
      print(p)
      
      
      dev.off()
      next
    }
    
    # compute highest densitiy intervals for parameters associated with selection only
    hdinterval.selection <- as.data.frame(t(hdi(fits.selection, credMass=0.8)))
    hdinterval.selection$Parameter <- rownames(hdinterval.selection) 
    hdinterval.selection$Median <- apply(fits.selection, 2, function(x){median(as.numeric(x))})
    hdinterval.selection$Sample <- patient.id
    hdinterval.selection$depth <- depth
    
    
    p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
      facet_wrap(~Parameter, scales="free") 
    
    print(p)
    
    
    parameters <- rbind(parameters, hdinterval)
    
    selected.parameters <- rbind(selected.parameters, hdinterval.selection)
    
    # compute highest density intervals for parameters associated with neutral evolution only
    if(nrow(fits.neutral)>=0){
      hdi.neutral <- as.data.frame(t(hdi(fits.neutral, credMass = 0.8)))
      hdi.neutral$Parameter <- rownames(hdi.neutral) 
      hdi.neutral$Median <- apply(fits.neutral, 2, function(x){median(as.numeric(x))})
      hdi.neutral$Sample <- patient.id
      hdi.neutral$depth <- depth
      
      neutral.parameters <- rbind(neutral.parameters, hdi.neutral[c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "par_offset",
                                                                    "mutations_per_year"),])
      
    }
    
    dev.off()
    
  }
  
}

# store the inferred parameters and model statistics
save(parameters, neutral.parameters, selected.parameters, plotlist.model.vs.data, model.support.selection, inference.stats,
     file="RData/Simulated_data/Cohort_parameters.RData")


######################### ######################### ######################### ######################### ######################### 
### Figure 1g: Posterior probability for the selected clone

to.plot <- melt(t(model.support.selection), value.name = "P_selection")
colnames(to.plot)[c(1,2)] <- c("Patient", "Depth")
to.plot$P_neutral <- 100 - to.plot$P_selection
to.plot <- melt(to.plot, value.name = "Posterior probability", measure.vars=c("P_selection", "P_neutral"))
to.plot$Patient <- as.character(to.plot$Patient)
to.plot <- to.plot[!to.plot$Patient %in% paste0("STN", 11:20),]
to.plot$Clone_size <- apply(to.plot, 1, function(x){
  res <- selected.parameters[selected.parameters$Sample==x[1] & selected.parameters$Parameter=="size_of_clone" &
                                as.character(x[3])=="P_selection" & selected.parameters$depth == as.numeric(x["Depth"]),]$Median
  if(length(res)==0){
    return(0)}else{
      return(res)
    }
})

to.plot$variable <- factor(to.plot$variable, levels=c("P_neutral", "P_selection"))
to.plot$Clone_size[to.plot$Clone_size==0] <- NA

## plot posterior probabilities supporting the neutral and the selection model for a subset of examples:

pdf(paste0(analysis.directory, "/Figures/Figure_1_g.pdf"), width=6, height=6)

to.plot.selection <- to.plot[to.plot$Depth==90 & !is.na(to.plot$`Posterior probability`) &
                               to.plot$Patient %in% c("STN8", "STS25", "STS35", "STS46", "STS53", "STS84"),]

ggplot(to.plot.selection, 
       aes(x=Patient, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("90x")+ theme( strip.background = element_blank() )

dev.off()

######################### ######################### ######################### ######################### ######################### 
### Figure 1h/S1d-f: ROC curve for varying thresholds

# 90x seq depth
to.plot <- inference.stats[["90"]]
colnames(to.plot) <- c("Selection", "Neutral", "Size_true_clone", "Inferred_size")
to.plot <- to.plot[to.plot$Selection + to.plot$Neutral>0,]
to.plot <- to.plot[!rownames(to.plot) %in% paste0("STN", 11:20),]
## to build the ROC curve, use different cutoffs for the evidence, ranging from 5% to 100% in 5% steps:

roc <- data.frame()

for(cutoff in seq(1, 100, 1)){
  for(size in setdiff(unique(to.plot$Size_true_clone), 0)){
    
    ## I. True positive rate: true positives/(true positives + false negatives)
    ## all cases with this clone size
    all.cases.w.clone.size <- sum(to.plot$`Size_true_clone`==size)
    if(all.cases.w.clone.size==0){next}
    true.cases.w.clone.size <- sum(to.plot$`Size_true_clone`==size & to.plot$Inferred_size<=1.5*size & 
                                     to.plot$Inferred_size >= 0.5*size & to.plot$Selection/10 >= cutoff)
    
    tp <- true.cases.w.clone.size/all.cases.w.clone.size
    
    ## II. False positive rate: false positives/(false positives + true negatives)
    false.cases.w.clone.size <- sum(to.plot$`Size_true_clone`==0 & to.plot$Inferred_size==size & to.plot$Selection/10 >= cutoff )
    true.negatives.w.clone.size <- sum((to.plot$`Size_true_clone`==0 & to.plot$Inferred_size!=size  )|
                                         (to.plot$`Size_true_clone`==0 & to.plot$Inferred_size==size & to.plot$Selection/10 < cutoff))
    
    fp <- false.cases.w.clone.size/(false.cases.w.clone.size + true.negatives.w.clone.size)
    
    roc <- rbind(roc, data.frame(TP=tp, FP=fp, Cutoff = cutoff, size = size ))
    
  }
}

roc <- roc[order(roc$TP),]
roc <- roc[order(roc$FP),]

# plot TP against FP
p1 <- ggplot(roc, aes(x=FP, y=TP, col=size, group=size)) +  geom_line() + 
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(5, 10, 25, 50), trans="log10", limits=c(2.5, 100)) +
  geom_abline(slope = 1, intercept = 0, linetype=2) + scale_x_continuous(limits=c(0,1)) + ggtitle("90x") +
  geom_point(data=roc[roc$Cutoff==15,], aes(x=FP, y=TP), col="firebrick")

# plot TP against the cutoff
p2 <- ggplot(roc[roc$size>=5,], aes(x=Cutoff, y=TP, col=size, group=size)) + geom_line() + 
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(5, 10, 25, 50), trans="log10", limits=c(2.5, 100)) +  ggtitle("90x")+
  geom_point(data=roc[roc$Cutoff==15 & roc$size>=5,], aes(x=Cutoff, y=TP), col="firebrick")

# plot FP against the cutoff
p3 <- ggplot(roc[roc$size>=5,], aes(x=Cutoff, y=FP, col=size, group=size)) +  geom_line() + 
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(5, 10, 25, 50), trans="log10", limits=c(2.5, 100)) +  ggtitle("90x")+
  geom_point(data=roc[roc$Cutoff==15 & roc$size>=5,], aes(x=Cutoff, y=FP), col="firebrick")

# plot the difference between TP and FP
p4 <- ggplot(roc[roc$size>=5,], aes(x=Cutoff, y=TP - FP, col=size, group=size)) + geom_line() + geom_line(stat = "summary", linetype=2, aes(x=Cutoff, y=TP - FP), inherit.aes = F) + #geom_smooth(aes(x=Cutoff, y=TP - FP), inherit.aes = F) +
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(5, 10, 25, 50), trans="log10", limits=c(2.5, 100)) +  ggtitle("90x")+
  geom_vline(xintercept = 15, col="firebrick")
  
# plot the average difference between TP and FP - choose the operating point at the maximum
p5 <- ggplot(roc[roc$size>=5,], aes(x=Cutoff, y=TP - FP, group=size)) + geom_line(stat = "summary", linetype=2, aes(x=Cutoff, y=TP - FP), inherit.aes = F) + #geom_smooth(aes(x=Cutoff, y=TP - FP), inherit.aes = F) +
  ggtitle("90x") + coord_cartesian(ylim=c(0,1))
  geom_vline(xintercept = 15, col="firebrick")

pdf(paste0(analysis.directory, "/Figures/Figure_1_h.pdf"), width=4, height = 3)

print(p5)

dev.off()

pdf(paste0(analysis.directory, "/Figures/Figure_1_i.pdf"), width=4, height = 3)

print(p1)

dev.off()


pdf(paste0(analysis.directory, "/Figures/Figure_S1_d.pdf"), width=4, height = 3)

print(p2)

dev.off()

pdf(paste0(analysis.directory, "/Figures/Figure_S1_e.pdf"), width=4, height = 3)

print(p3)

dev.off()

pdf(paste0(analysis.directory, "/Figures/Figure_S1_f.pdf"), width=4, height = 3)

print(p4)

dev.off()

######################### ######################### ######################### ######################### ######################### 
### Figure 1j, compute AUC; add point 1/1 to every combination

auc <- data.frame()

for(size in unique(roc.90$size)){
  tmp <- unique(roc.90[roc.90$size==size,c("TP", "FP"),drop=F])
  tmp <- unique(rbind(tmp, 1, 1))
  
  auc <- rbind(auc,
               data.frame(AUC=sum(rowMeans(cbind(tmp$TP[-1], tmp$TP[-length(tmp$TP)])) * (tmp$FP[-1] - tmp$FP[-length(tmp$FP)])),
                          size = size, depth=90))
}

for(size in unique(roc.30$size)){
  tmp <- unique(roc.30[roc.30$size==size,c("TP", "FP"),drop=F])
  tmp <- unique(rbind(tmp, 1, 1))
  
  auc <- rbind(auc,
               data.frame(AUC=sum(rowMeans(cbind(tmp$TP[-1], tmp$TP[-length(tmp$TP)])) * (tmp$FP[-1] - tmp$FP[-length(tmp$FP)])),
                          size = size, depth=30))
}

for(size in unique(roc.270$size)){
  tmp <- unique(roc.270[roc.270$size==size,c("TP", "FP"),drop=F])
  tmp <- unique(rbind(tmp, 1, 1))
  
  auc <- rbind(auc,
               data.frame(AUC=sum(rowMeans(cbind(tmp$TP[-1], tmp$TP[-length(tmp$TP)])) * (tmp$FP[-1] - tmp$FP[-length(tmp$FP)])),
                          size = size, depth=270))
}

pdf(paste0(analysis.directory, "/Figures/Figure_1_j.pdf"), width=4, height = 3)

ggplot(auc, aes(x=size, y=AUC, group=depth, linetype=as.character(depth))) + geom_line() + scale_x_log10(limits=c(1, 100)) +
  scale_y_continuous(limits = c(0.5,1)) + geom_vline(xintercept = 10, linetype=2)

dev.off()

