###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### libraries and functions
library(cdata)
library(scales)
library(openxlsx)
library(ggridges)
library(SCIFER)

source("./Settings.R")

## mutation data
patient.id <- 'KX004'
patient.nr <- which(sample.info.published.data$ID == patient.id)
load(paste0("./RData/", sample.info.published.data[patient.nr,]$SNVs_file_name))
############################################################################################################################################

age <- sample.info.published.data[patient.nr, ]$Age*365
ncells <- sample.info.published.data[patient.nr,]$ncells

depth=10000

min.vaf <- 0.02
min.clone.size = 0.01
min.prior.size=0.001

use.sensitivity = F

seq.type="sc"

snvs <- list(data.frame(VAF=snvs[[patient.id]], Depth=100, varCounts=snvs[[patient.id]]*100))

mother.daughter <- matrix(c(1, 2, 1,3), byrow=T, ncol=2)

directory <- paste0(analysis.directory, "/Model_fits_2_clones/Published_data/", sample.info.published.data[patient.nr,]$Path)


##########################################################################################################
## read in fits and analyze the parameters

fits <- read.csv(paste0(directory, "/Model_fit.csv"))

# compute selective advantage from size of the selected clone
fits$s1 <- 1 - log(10^fits$par_size1*10^fits$par_N)/(10^fits$par_lambda_ss*(age - fits$par_ts1))
fits$s2 <- 1 - log(10^fits$par_size2*10^fits$par_N)/(10^fits$par_lambda_ss*(age - fits$par_ts2))


fits$mutations_per_year <- fits$par_mu*10^fits$par_lambda_ss*365

true_clone_sizes <- apply(fits, 1, function(x){
  res <- SCIFER:::.compute_actual_size (t.s = c(0, as.numeric(x["par_ts1"]), as.numeric(x["par_ts2"])), mother.daughter = mother.daughter, 
                                        N = 10^as.numeric(x["par_N"]), lambda.ss = 10^as.numeric(x["par_lambda_ss"]), 
                                        delta.ss = 10^as.numeric(x["par_delta_ss"]),
                                        lambda.exp = 1, delta.exp = as.numeric(x["par_delta_exp"]), 
                                        size = c(10^as.numeric(x["par_size1"]), 10^as.numeric(x["par_size2"])), t.end = age)
  res[-1]/10^as.numeric(x["par_N"])
})

fits$true_size_1 <- true_clone_sizes[1,]
fits$true_size_2 <- true_clone_sizes[2,]

fits$age_of_clone_1 <- (age - fits$par_ts1)/365
fits$age_of_clone_2 <- (age - fits$par_ts2)/365

## sort clones by size
order.fits <- t(apply(fits, 1, function(x){
  order(x[c("true_size_1", "true_size_2")], decreasing = T)
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
fits[,c("true_size_1", "true_size_2")] <- t(apply(cbind(fits[,c("true_size_1", "true_size_2")], order.fits), 1, function(x){
  x[c(1,2)][x[c(3,4)]]
}))
fits[,c("par_ts1", "par_ts2")] <- t(apply(cbind(fits[,c("par_ts1", "par_ts2")], order.fits), 1, function(x){
  x[c(1,2)][x[c(3,4)]]
}))

# to compute the average growth per year, go from 1 cell to the final clone size
fits$growth_per_year1 <- apply(fits, 1, function(x){
  # effective selective advantage
  seff = 1 - log(as.numeric(x["true_size_1"])*10^as.numeric(x["par_N"]))/(10^as.numeric(x["par_lambda_ss"])*(age - as.numeric(x["par_ts1"])))
  # effective growth rate per year
  exp(10^as.numeric(x["par_lambda_ss"])*(1 - seff)*365) - 1
})
fits$growth_per_year2 <- apply(fits, 1, function(x){
  # effective selective advantage
  seff = 1 - log(as.numeric(x["true_size_2"])*10^as.numeric(x["par_N"]))/(10^as.numeric(x["par_lambda_ss"])*(age - as.numeric(x["par_ts2"])))
  # effective growth rate per year
  exp(10^as.numeric(x["par_lambda_ss"])*(1 - seff)*365) - 1
})

## compute Nxtau

fits$N_tau <- 10^fits$par_N / (365 * 10^fits$par_lambda_ss)

modelResult <- list()


to.plot <- fits[,c("par_N", "par_delta_exp", "par_lambda_ss", "par_mu", "N_tau", "par_offset", "par_ts1",
                   "par_ts2", "s1", "s2", "true_size_1", "true_size_2", "growth_per_year1", "growth_per_year2",
                   "mutations_per_year")]

# highest density interval, all parameters
hdinterval <- as.data.frame(t(hdi(fits, credMass=0.8)))
hdinterval$Parameter <- rownames(hdinterval)
hdinterval$Median <- apply(fits, 2, function(x){median(as.numeric(x))})


### selection parameters
fits.selection.1 <- fits[fits$true_size_1 > 2*min.vaf,]
fits.selection.2 <- fits[fits$true_size_2 > 2*min.vaf,]

hdinterval.selection.1 <- as.data.frame(t(hdi(fits.selection.1, credMass=0.8)))
hdinterval.selection.1$Parameter <- rownames(hdinterval.selection.1)
hdinterval.selection.1$Median <- apply(fits.selection.2, 2, function(x){median(as.numeric(x))})

hdinterval.selection.2 <- as.data.frame(t(hdi(fits.selection.2, credMass=0.8)))
hdinterval.selection.2$Parameter <- rownames(hdinterval.selection.2)
hdinterval.selection.2$Median <- apply(fits.selection.2, 2, function(x){median(as.numeric(x))})

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

p9 <- ggplot(fits, aes(x=10^par_N, y=10^par_size1)) +
  geom_density_2d_filled(col=NA, contour_var = "ndensity",  aes( fill = ..level..)) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
  scale_x_log10(name="N", limits=c(10^2.5, 10^8)) + scale_y_continuous(name="size1", limits = c(0.01, 1))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")


print(ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow=3, ncol=3))



parms.of.interest <- c("par_N", "par_lambda_ss", "par_delta_exp", "par_mu",
                       "par_size1", "par_size2", "par_ts1", "par_ts2")

parameters.to.plot <- data.frame(parms_i = rep(parms.of.interest, each = length(parms.of.interest)),
                                 parms_j = rep(parms.of.interest, length(parms.of.interest)),
                                 i = rep(1:length(parms.of.interest), each = length(parms.of.interest)),
                                 j = rep(1:length(parms.of.interest), length(parms.of.interest)))

plot.design <- data.frame(Parameter = parms.of.interest,
                          Min = c(2.5, -3, 0, 0.1, -3, -3, 0, 0),
                          Max = c(8, -1, 0.75, 10, 0, 0, age, age),
                          Label = c("N (log10)", "lambda_ss (log10)", "delta/lambda",
                                    "mu (1/div)", "size of clone 1", "size of clone2",
                                    "age of clone1 (years)", "age of clone2 (years)"))


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

p <- ggplot(data=hdinterval, aes(x=Parameter, y = Median, ymin=lower, ymax=upper)) + geom_pointrange() +
  facet_wrap(~Parameter, scales="free")

print(p)


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Extended Data Fig. 3e: plot model fit

source(paste0(custom.script.directory, "/Parameter_estimation/Bayesian_fit_multiclone.R"))

mother.daughter <- matrix(c(1,2,
                            1,3), byrow=T, ncol = 2)

if( !file.exists(paste0(directory, "/Sim_trajectories.RData"))){
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
  
  save(sim, data.vs.prediction, file=paste0(directory, "/Sim_trajectories.RData"))
  
}else{
  load(paste0(directory, "/Sim_trajectories.RData"))
}

to.plot <- data.vs.prediction


pdf(paste0(analysis.directory, "/Figures/Ext_Data_Fig_3j_fit.pdf"), width=3.5, height=3.5)

max.y <- max(to.plot$max.model)

xlimits <- c(0,1/min.vaf)

p <- ggplot(data=to.plot, aes(x=1/VAF, y=mean.data, ymin=mean.data-sd.data, ymax=mean.data+sd.data)) +
  geom_ribbon(data=to.plot, aes(x=1/VAF, y=mean.data, ymin=min.model,ymax=max.model), alpha=1,
              fill="violet") +
  geom_pointrange(lwd=0.25, shape=1, fatten=1) +
  scale_y_continuous(name="Cumulative # of mutations") +
  theme(aspect.ratio = 1) +
  scale_x_continuous( breaks = c(5, 10, 20, 50, 100), labels = c("0.2", "0.1", "0.05", "0.02", "0.01"), name = "Variant allele frequency") +
  coord_cartesian(ylim=c(0, max.y*1.05), xlim=xlimits) 


p <- p + geom_ribbon(data = data.frame(x = unlist(hdinterval.selection.1[hdinterval.selection.1$Parameter=="true_size_1",c("lower", "upper")]),
                                       ymin = c(0,0),
                                       ymax = c(max.y, max.y)), aes(x=2/x, ymin=ymin, ymax = ymax), inherit.aes = F, fill="darkblue", alpha = 0.5)


p <- p + geom_ribbon(data = data.frame(x = unlist(hdinterval.selection.2[hdinterval.selection.2$Parameter=="true_size_2",c("lower", "upper")]),
                                       ymin = c(0,0),
                                       ymax = c(max.y, max.y)), aes(x=2/x, ymin=ymin, ymax = ymax), inherit.aes = F, fill="darkgreen", alpha = 0.5)

print(p)

dev.off()

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Extended Data Fig. 3j: posterior probability for 2nd clone

# Model evidence for 1st and 2nd clone
pdf(paste0(analysis.directory, "/Figures/Ext_Data_Fig_3j_posterior.pdf"), width=3.5, height=3.5)

to.plot <- data.frame(Clone = rep(c(1,2), each = 2),
                      Probability = c(nrow(fits.selection.1)/10, 100 - nrow(fits.selection.1)/10,
                                      nrow(fits.selection.2)/10, 100 - nrow(fits.selection.2)/10),
                      Mode = rep(c("selected", "neutral"), 2),
                      Size = c(hdinterval.selection.1["true_size_1",]$Median, NA,
                               hdinterval.selection.2["true_size_2",]$Median, NA))

ggplot(to.plot, aes(x = Clone, y = Probability, fill = Size)) + geom_col() + 
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  geom_hline(yintercept = 15, linetype = 2)

dev.off()

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Extended Data Fig. 3k, compare parameter estimates with single-clone model

to.plot <- data.frame(Parameter = c("par_N", "par_lambda_ss", "par_mu", "age_of_clone", "growth_of_clone"),
                      Model = "2-clone",
                      Median = hdinterval.selection.2[c("par_N", "par_lambda_ss", "par_mu",
                                                        "age_of_clone_1", "growth_per_year1"),]$Median,
                      Min = hdinterval.selection.2[c("par_N", "par_lambda_ss", "par_mu",
                                                     "age_of_clone_1", "growth_per_year1"),]$lower,
                      Max = hdinterval.selection.2[c("par_N", "par_lambda_ss", "par_mu",
                                                     "age_of_clone_1", "growth_per_year1"),]$upper)

load(paste0(rdata.directory, "Cohort_parameters_published_data.RData"))

to.plot <- rbind(to.plot,
                 data.frame(Parameter = c("par_N", "par_lambda_ss", "par_mu", "age_of_clone", "growth_of_clone"),
                            Model = "1-clone",
                            Median = sapply(c("par_N", "par_lambda_ss", "par_mu", "age_of_clone", "growth_per_year"),function(x){
                              selected.parameters[selected.parameters$Parameter==x &
                                                    selected.parameters$Sample==patient.id &
                                                    selected.parameters$Resolution==0.01,]$Median
                            }),
                            Min = sapply(c("par_N", "par_lambda_ss", "par_mu", "age_of_clone", "growth_per_year"),function(x){
                              selected.parameters[selected.parameters$Parameter==x &
                                                    selected.parameters$Sample==patient.id&
                                                    selected.parameters$Resolution==0.01,]$lower
                            }),
                            Max = sapply(c("par_N", "par_lambda_ss", "par_mu", "age_of_clone", "growth_per_year"),function(x){
                              selected.parameters[selected.parameters$Parameter==x &
                                                    selected.parameters$Sample==patient.id&
                                                    selected.parameters$Resolution==0.01,]$upper
                            })))

pdf(paste0(analysis.directory, "/Figures/Ext_Data_Fig_3k.pdf"), width=3.5, height=2.5)

ggplot(to.plot[to.plot$Parameter=="par_N",], aes(x = Model, y = Median, ymin = Min, ymax = Max)) +
  geom_pointrange() + scale_y_continuous(name = "N (log10)", limits = c(0, 8))+ theme(aspect.ratio = 1.5)

ggplot(to.plot[to.plot$Parameter=="par_lambda_ss",], aes(x = Model, y = 365*10^Median, ymin = 365*10^Min, ymax = 365*10^Max)) +
  geom_pointrange() + scale_y_log10(name = "Division rate per year", limits = c(0.3, 30))+ theme(aspect.ratio = 1.5)

ggplot(to.plot[to.plot$Parameter=="par_mu",], aes(x = Model, y = Median, ymin = Min, ymax = Max)) +
  geom_pointrange() + scale_y_continuous(name = "SSNVs/division")+ theme(aspect.ratio = 1.5)

ggplot(to.plot[to.plot$Parameter=="age_of_clone",], aes(x = Model, y = Median, ymin = Min, ymax = Max)) +
  geom_pointrange() + scale_y_continuous(name = "Age of leading clone (years)", limits = c(0, 100))+ theme(aspect.ratio = 1.5)

ggplot(to.plot[to.plot$Parameter=="growth_of_clone",], aes(x = Model, y = 100*Median, ymin = 100*Min, ymax = 100*Max)) +
  geom_pointrange() + scale_y_continuous(name = "clonal growth per year (%, leading clone)") + theme(aspect.ratio = 1.5) + expand_limits(y=0)

ggplot(hdinterval.selection.2["age_of_clone_2",], aes(x = "2-clone", y = Median, ymin = lower, ymax = upper)) +
  geom_pointrange() + scale_y_continuous(name = "Age of second clone (years)", limits = c(0, 100)) + theme(aspect.ratio = 2)

ggplot(hdinterval.selection.2["growth_per_year2",], aes(x = "2-clone", y = 100*Median, ymin = 100*lower, ymax = 100*upper)) +
  geom_pointrange() + scale_y_log10(name = "clonal growth per year (%, leading clone)")+ theme(aspect.ratio = 2) +
  expand_limits(y = 0)

dev.off()
