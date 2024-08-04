source("./Settings.R")

####################################################################################################################################################
## In the paper, we show parameter estimates obtained with both the one-clone model and the 2-clone model.
## We first source the file "Assess_fits_heme_WGS_data.R" to get the parameters obtained with the 1-clone model.
## In a second step, we source the file "Assess_fits_heme_WGS_data_2_clone_model.R" to load the parameters obtained with the 2-clone model.

## Collect all parameters from the one-clone and two-clone model. For each case
## distinguish whether it was neutrally evolving, had one selected clone or 2 selected clones and store the parameters corresponding to the respective model.

## Parameters estimated with the one-clone model
source(paste0(custom.script.directory, "/Analysis_and_figures/Assess_fits_heme_WGS_data.R"))

## Posterior probability for selection, 1-clone model
one.clone.model.support.selection <- model.support.selection

## Cases classified as selected 
one.clone.model.parameters.selected <- selected.parameters.1cm[(selected.parameters.1cm$Paper_ID %in% selected.samples.90 &
                                                                  selected.parameters.1cm$Sort == "CD34") |
                                                             (selected.parameters.1cm$Paper_ID %in% selected.samples.270 &
                                                                selected.parameters.1cm$Sort == "CD34_deep"),]
one.clone.model.parameters.selected$Subclones <- "1 clone"
one.clone.model.parameters.selected$Mode <- ""


## Cases classified as neutrally evolving
one.clone.model.parameters.neutral <- neutral.parameters.1cm[(neutral.parameters.1cm$Paper_ID %in% neutral.samples.90 &
                                                                neutral.parameters.1cm$Sort == "CD34") |
                                                           (neutral.parameters.1cm$Paper_ID %in% neutral.samples.270 &
                                                              neutral.parameters.1cm$Sort == "CD34_deep"),]
one.clone.model.parameters.neutral$Subclones <- "Neutral"
one.clone.model.parameters.neutral$Mode <- ""

## Parameters estimated with the 2-clone model 
source(paste0(custom.script.directory, "/Analysis_and_figures/Assess_fits_heme_WGS_data_2_clone_model.R"))

## Cases where a 2nd clone was identified
two.clone.model.parameters <- selected.parameters.2cm.2clones[paste(selected.parameters.2cm.2clones$Paper_ID, selected.parameters.2cm.2clones$Mode, sep="_") %in%
                                                                colnames(model.support.selection.2_clones[,model.support.selection.2_clones["Clone 2",] > 15]),]
two.clone.model.parameters$Subclones <- "2 clones"
## posterior probability of 2nd clone
two.clone.model.parameters$Posterior.2 <- model.support.selection.2_clones["Clone 2",paste(two.clone.model.parameters$Paper_ID, two.clone.model.parameters$Mode, sep="_")]


## Now, collect all parameters
## - for those cases with 2 clones: take the topology with higher posterior probability; if both topologies fit equally well, take the average
## - for those cases with only 1 clone: take the parameters estimated with the one-clone model
## - for those cases which evolve neutrally: take the parameters estimated with the neutral model

all.cd34.parameters <- data.frame()

for(i in unique(two.clone.model.parameters$Paper_ID)){
  for(j in unique(two.clone.model.parameters$Parameter)){
    tmp <- two.clone.model.parameters[two.clone.model.parameters$Paper_ID == i &
                                        two.clone.model.parameters$Parameter==j,]
    if(nrow(tmp)==0){next}
    if(nrow(tmp)==2){
      if(tmp$Posterior.2[1]== tmp$Posterior.2[2]){
        all.cd34.parameters <- rbind(all.cd34.parameters,
                                data.frame(lower = mean(tmp$lower),
                                           upper = mean(tmp$upper),
                                           Parameter = j,
                                           Median = mean(tmp$Median),
                                           Paper_ID = i,
                                           Subclones = tmp$Subclones[1]))
      }else{
        all.cd34.parameters <- rbind(all.cd34.parameters,
                                data.frame(lower = tmp$lower[which.max(tmp$Posterior.2)],
                                           upper = tmp$upper[which.max(tmp$Posterior.2)],
                                           Parameter = j,
                                           Median = tmp$Median[which.max(tmp$Posterior.2)],
                                           Paper_ID = i,
                                           Subclones = tmp$Subclones[which.max(tmp$Posterior.2)]))
      }
    }else{
      all.cd34.parameters <- rbind(all.cd34.parameters,
                              data.frame(lower = mean(tmp$lower),
                                         upper = mean(tmp$upper),
                                         Parameter = j,
                                         Median = mean(tmp$Median),
                                         Paper_ID = i,
                                         Subclones = tmp$Subclones[1]))
    }
    
  }
}

all.cd34.parameters$Depth <- 270


all.cd34.parameters <- rbind(all.cd34.parameters, one.clone.model.parameters.neutral[,colnames(all.cd34.parameters)],
                        one.clone.model.parameters.selected[!(one.clone.model.parameters.selected$Paper_ID %in% all.cd34.parameters$Paper_ID &
                                                                one.clone.model.parameters.selected$Sort=="CD34_deep"),colnames(all.cd34.parameters)]) # take 2-clone model parameters wherever applicable

all.cd34.parameters$Depth <- replace(all.cd34.parameters$Depth, all.cd34.parameters$Depth==300, 270) # combine cases with 270x/300x WGS
all.cd34.parameters$Depth <- replace(all.cd34.parameters$Depth, all.cd34.parameters$Depth==120, 90) # combine cases with 90x/120x WGS

rm(two.clone.model.parameters)
rm(one.clone.model.parameters.neutral)
rm(one.clone.model.parameters.selected)

####################################################################################################################################################
## Figures 4a/b, 5d, Extended Data Fig. 7a-c, 8a/b: plot the model fits stratified by type

# no known CH driver, 90x:
pdf(paste0(analysis.directory, "/Figures/Ext_Data_Fig_8a.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[c("1-N CD34+", "2-U CD34+", "3-N CD34+", 
                                            "4-N CD34+", "8-UU CD34+", "11-UU CD34+",
                                            "13-U CD34+", "14-U CD34+", "15-N CD34+", "16-UU CD34+")],
          nrow=3, ncol=3)

dev.off()

# no known CH driver, 270x:
pdf(paste0(analysis.directory, "/Figures/Figure_5d_Ext_Data_Fig_8b.pdf"), width=8, height=8)

# cases with 1 selected clone
ggarrange(plotlist=plotlist.model.vs.data[c("1-N CD34+", "3-N CD34+", "4-N CD34+", 
                                            "13-U CD34+", "14-U CD34+", "15 N_CD34+")],
          nrow=3, ncol=3)

# cases with 2 selected clones
ggarrange(plotlist=plotlist.model.vs.data.2cm[c("8-UU_branched", "11-UU_branched", "16-UU_branched")],
          nrow=3, ncol=3)

dev.off()


# known CH driver, 90x:
pdf(paste0(analysis.directory, "/Figures/Ext_Data_Fig_7a.pdf"), width=8, height=8)

ggarrange(plotlist=plotlist.model.vs.data[c("5-DU CD34+", "6-AU CD34+", "7-T CD34+", 
                                            "8-UU CD34+", "9-DU CD34+", "10-D CD34+", 
                                            "12-AT CD34+", "17-T CD34+", "18-DU CD34+", 
                                            "19-D CD34+", "20-UT CD34+", "21-DU CD34+", "22-TU CD34+")],
  nrow=3, ncol=3)

dev.off()

# known CH driver, 270x
pdf(paste0(analysis.directory, "/Figures/Figure_4ab_Ext_Data_Fig_7b.pdf"), width=8, height=8)

# cases with 1 selected clone
ggarrange(plotlist=plotlist.model.vs.data[c("7-T CD34+_deep", "10-D CD34+_deep", "17-T CD34+_deep")],
  nrow=3, ncol=3)

# cases with 2 selected clones
ggarrange(plotlist=plotlist.model.vs.data.2cm[c("5-DU_branched", "6-AU_branched", "9-DU_linear", "12-AT_linear",
                                            "18-DU_branched", "20-UT_branched", "22-TU_linear")],
          nrow=3, ncol=3)

dev.off()

## for sample 21-DU, plot a zoom in:

pdf(paste0(analysis.directory, "/Figures/Ext_Data_Fig_7c.pdf"), width=3.5, height=3.5)

plotlist.model.vs.data$`21-DU CD34+` +
  scale_x_continuous( breaks = c(1/0.5, 1/0.4, 1/0.3, 1/0.2), labels = c("0.5", "0.4", "0.3", "0.2"), name="Variant allele frequency") +
  scale_y_continuous( breaks=seq(0,15), labels=c("0", "", "", "", "", "5", "", "", "", "", "10", "", "", "", "", "15"), name = "Number of SSNVs") +
  coord_cartesian(ylim=c(0, 15), xlim=c(1/0.5, 1/0.2))

dev.off()

## estimate time point of clone emergence: ~5-6 SSNVs

n.div <- 6/selected.parameters.1cm[selected.parameters.1cm$Paper_ID=="21-DU" & selected.parameters.1cm$Sort=="CD34" &
                                     selected.parameters.1cm$Parameter=="par_mu" ,c("lower", "Median", "upper")]

############################################################################################################################################
## Figs. 5a: plot the SNV data for the cases without known CH driver

to.plot <- data.frame(VAF= c(), MutationCount=c(), Sample=c(), Age=c())
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
                                       Depth = 90))
}

## 270x
vafs.of.interest <- seq(0.02, 1, 0.01)

for(i in sample.info[sample.info$CH.driver.found=="no" &
                     sample.info$`Coverage.WGS.CD34+.2`>0,]$Paper_ID){
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

pdf(paste0(analysis.directory, "/Figures/Figure_5a.pdf"), width=6, height=5)

# stratify by age
to.plot$BinAge <- floor(to.plot$Age/10)*10
ggplot(to.plot, aes(x=1/VAF, y=MutationCount, col=Age, group=paste(Paper_ID, Depth), 
                    linetype = factor(Depth))) + geom_point() + geom_line() +
  scale_color_gradientn(colors=hcl.colors(n=7, palette="Zissou 1"), limits=c(25, 80)) +
  scale_x_continuous(breaks = c(1/0.2, 1/0.1, 1/0.05, 1/0.02), labels = c("0.2", "0.1", "0.05", "0.02"), name="Variant allele frequency") +
  theme(aspect.ratio = 0.02/0.05)+ ggtitle("CD34") + facet_wrap(~BinAge, scales = "free_y", ncol = 2)


dev.off()


############################################################################################################################################
## Figs. 4c/d/e, 5b/c: plot the posterior probability for the leading and second selected clone, stratified by coverage

####  I. 1st clone 
# Get the posterior probability of selection from the 1-clone model
to.plot.1st.clone <- melt(t(one.clone.model.support.selection), 
                value.name = "P_selection")
colnames(to.plot.1st.clone)[c(1,2)] <- c("Paper_ID", "Sort")
to.plot.1st.clone$P_neutral <- 100 - to.plot.1st.clone$P_selection
to.plot.1st.clone <- melt(to.plot.1st.clone, value.name = "Posterior probability", variable.name = "Type")
to.plot.1st.clone$Sort <- as.character(to.plot.1st.clone$Sort)
to.plot.1st.clone$Depth <- ifelse(to.plot.1st.clone$Sort == "CD34_deep", 270, 90)
to.plot.1st.clone$Paper_ID <- as.character(to.plot.1st.clone$Paper_ID)
to.plot.1st.clone$Clone_size <- apply(to.plot.1st.clone, 1, function(x){

  res <- selected.parameters.1cm[selected.parameters.1cm$Paper_ID==x["Paper_ID"] & selected.parameters.1cm$Parameter=="size_of_clone" &
                                    x["Type"]=="P_selection" &
                                   selected.parameters.1cm$Sort == x["Sort"],]$Median

  if(length(res)==0){
    return(0)}else{
      return(res)
    }
})

to.plot.1st.clone$Sort <- factor(to.plot.1st.clone$Sort, levels = c("CD34", "CD34_deep", "MNC-T", "MNC", "PB"))
to.plot.1st.clone$Type <- factor(to.plot.1st.clone$Type, levels=c("P_neutral", "P_selection"))
to.plot.1st.clone$Clone_size[to.plot.1st.clone$Clone_size==0] <- NA
to.plot.1st.clone$Known_CHIP_driver <- sample.info[as.character(to.plot.1st.clone$Paper_ID),]$CH.driver.found

to.plot.1st.clone <- to.plot.1st.clone[order(to.plot.1st.clone$Sort, decreasing = T),]
to.plot.1st.clone <- to.plot.1st.clone[!is.na(to.plot.1st.clone$`Posterior probability`),]


#### II., 2nd clone
to.plot.2nd.clone <- melt(t(model.support.selection.2_clones), value.name = "P_selection")
colnames(to.plot.2nd.clone)[c(1,2)] <- c("Fit_ID", "Clone")
to.plot.2nd.clone$Paper_ID <- gsub("_.*", "", to.plot.2nd.clone$Fit_ID)
## take maximum support for second clone per ID
to.plot.2nd.clone <- to.plot.2nd.clone[to.plot.2nd.clone$Clone == "Clone 2",]
to.plot.2nd.clone <- to.plot.2nd.clone[order(to.plot.2nd.clone$P_selection, decreasing = T),]
to.plot.2nd.clone <- to.plot.2nd.clone[order(to.plot.2nd.clone$Paper_ID),]
to.plot.2nd.clone <- to.plot.2nd.clone[!duplicated(to.plot.2nd.clone$Paper_ID),]
to.plot.2nd.clone$P_neutral <- 100 - to.plot.2nd.clone$P_selection
to.plot.2nd.clone <- melt(to.plot.2nd.clone, value.name = "Posterior probability", id.vars = c("Fit_ID", "Clone", "Paper_ID"),
                          measure.vars = c("P_selection", "P_neutral"), variable.name = "Type")
to.plot.2nd.clone$Clone_size <- apply(to.plot.2nd.clone, 1, function(x){
  mode <- gsub(".*_", "", x[1])
  res <- selected.parameters.2cm.2clones[selected.parameters.2cm.2clones$Paper_ID == x["Paper_ID"] &
                                       selected.parameters.2cm.2clones$Mode == mode &
                                       selected.parameters.2cm.2clones$Parameter=="clone_size_2" &
                                       x["Type"] == "P_selection",]$Median
  if(length(res)==0){
    return(0)
  }else{
    return(res)
  }
})
to.plot.2nd.clone$Type <- factor(to.plot.2nd.clone$Type, levels=c("P_neutral", "P_selection"))
to.plot.2nd.clone$Clone_size[to.plot.2nd.clone$Clone_size==0] <- NA
to.plot.2nd.clone$Known_CHIP_driver <- sample.info[as.character(to.plot.2nd.clone$Paper_ID),]$CH.driver.found
to.plot.2nd.clone <- to.plot.2nd.clone[!grepl("N", as.character(to.plot.2nd.clone$Paper_ID)),]
to.plot.2nd.clone$Paper_ID <- factor(to.plot.2nd.clone$Paper_ID, 
                                     levels = sort(sample.info$Paper_ID[!grepl("N", sample.info$Paper_ID)]))

# Figure 4c/d/e: all cases with known CH drivers
pdf(paste0(analysis.directory, "/Figures/Figures_4cde.pdf"), width=6, height=3)

to.plot. <- to.plot.1st.clone[to.plot.1st.clone$Known_CHIP_driver == "yes" &
                                grepl("CD34", to.plot.1st.clone$Sort),]
to.plot.$Paper_ID <- factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(as.numeric(gsub("-.*", "", to.plot.$Paper_ID)))]))

# 90x, leading clone
ggplot(to.plot.[to.plot.$Depth==90,],
       aes(x=Paper_ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, CH driver, 90x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2) +
  scale_y_continuous("Posterior probability leading clone (%)")

# 270x, leading clone
ggplot(to.plot.[to.plot.$Depth==270,],
       aes(x=Paper_ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, CH driver, 270x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2) +
  scale_y_continuous("Posterior probability leading clone (%)") + scale_x_discrete(drop = F)

## 270x, 2nd clone

to.plot. <- to.plot.2nd.clone[to.plot.2nd.clone$Known_CHIP_driver == "yes",]
to.plot.$Paper_ID <- factor(to.plot.$Paper_ID, 
                         levels = unique(sample.info$Paper_ID[sample.info$CH.driver.found == "yes"][
                           order(as.numeric(gsub("-.*", "", sample.info$Paper_ID[sample.info$CH.driver.found == "yes"])))]))

ggplot(to.plot.,
       aes(x=Paper_ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") + 
  scale_x_discrete(drop = F) +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  scale_y_continuous("Posterior probability 2nd clone (%)") +
  ggtitle("CD34, CH driver, 270x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

dev.off()

# Figure 5b/c: all cases without known CH drivers
pdf(paste0(analysis.directory, "/Figures/Figures_5bc.pdf"), width=6, height=3)

to.plot. <- to.plot.1st.clone[to.plot.1st.clone$Known_CHIP_driver == "no" &
                                grepl("CD34", to.plot.1st.clone$Sort),]
to.plot.$Paper_ID <- factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(as.numeric(gsub("-.*", "", to.plot.$Paper_ID)))]))

# 90x, leading clone
ggplot(to.plot.[to.plot.$Depth==90,],
       aes(x=Paper_ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, no CH driver, 90x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2) +
  scale_y_continuous("Posterior probability (%)")

# 270x, leading clone
ggplot(to.plot.[to.plot.$Depth==270,],
       aes(x=Paper_ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  ggtitle("CD34, no CH driver, 270x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2) +
  scale_y_continuous("Posterior probability (%)") + scale_x_discrete(drop = F)

## 270x, 2nd clone

to.plot. <- to.plot.2nd.clone[to.plot.2nd.clone$Known_CHIP_driver == "no",]
to.plot.$Paper_ID <- factor(to.plot.$Paper_ID, 
                         levels = unique(sample.info$Paper_ID[sample.info$CH.driver.found == "no" & !grepl("N", sample.info$Paper_ID)][
                           order(as.numeric(gsub("-.*", "", sample.info$Paper_ID[sample.info$CH.driver.found =="no" & !grepl("N", sample.info$Paper_ID)])))]))

ggplot(to.plot.,
       aes(x=Paper_ID, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") + 
  scale_x_discrete(drop = F) +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  scale_y_continuous("Posterior probability 2nd clone (%)") +
  ggtitle("CD34, no CH driver, 270x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)

dev.off()



# Supplementary Figure 3: all cases profiled for multiple sorts


pdf(paste0(analysis.directory, "/Figures/Suppl_Fig_3.pdf"), width=6, height=6)

samples.w.all.sorts <- sample.info[sample.info$Coverage.WGS.PB.granulocytes!="Not tested",]$Paper_ID
ggplot(to.plot.1st.clone[to.plot.1st.clone$Paper_ID %in% samples.w.all.sorts,],
       aes(x=Sort, y=`Posterior probability`, fill=Clone_size)) + geom_col(width=0.5, col="black") + 
  scale_x_discrete(drop = F) + facet_wrap(~Paper_ID) +
  scale_fill_gradientn(colors=hcl.colors(n = 7, palette="Zissou 1"), breaks=c(0.05, 0.10, 0.25, 0.50), trans="log10", limits=c(0.025, 1)) +
  scale_y_continuous("Posterior probability 2nd clone (%)") +
  ggtitle("CD34, no CH driver, 270x")+ theme( strip.background = element_blank() ) + geom_hline(yintercept = 15, linetype=2)


dev.off()

############################################################################################################################################
## Fig. 4g/h, Fig. 5e-h Plot the physiological parameters (N, lambda, mu) estimated with 90x/270x from CD34 cells side-by-side

to.plot <- all.cd34.parameters[all.cd34.parameters$Parameter %in% c("N_tau", "par_N", "par_mu", "par_lambda_ss"),]
to.plot$Known_CHIP.mutation <- sample.info[as.character(to.plot$Paper_ID),]$CH.driver.found
to.plot$Leading.clone <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit
to.plot$Age <- floor(sample.info[to.plot$Paper_ID,]$Age)
to.plot <- to.plot[order(as.numeric(gsub("-.*", "", to.plot$Paper_ID))),]
to.plot <- to.plot[order(to.plot$Age),]
to.plot$Paper_ID <- factor(to.plot$Paper_ID,
                        levels=sort(unique(sample.info$Paper_ID[order(sample.info$Age)])))



pdf(paste0(analysis.directory, "/Figures/Figures_4gh_5e-h.pdf"), width=2.5, height=2.5)

## N x tau
to.plot. <- to.plot[to.plot$Parameter == "N_tau",]
to.plot.$Linetype <- factor(to.plot.$Depth, levels = c("90", "270"))

# neutrally evolving samples
to.plot.. <- to.plot.[grepl("N", to.plot.$Paper_ID),]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=log10(Median), ymin=log10(lower), ymax=log10(upper), 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="N x Tau (log10)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

# samples with unknown driver
to.plot.. <- to.plot.[grepl("U", to.plot.$Paper_ID) &
                        to.plot.$Known_CHIP.mutation == "no",]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=log10(Median), ymin=log10(lower), ymax=log10(upper), 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="N x Tau (log10)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

# samples with known CH drivers
to.plot.. <- to.plot.[to.plot.$Known_CHIP.mutation == "yes",]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=log10(Median), ymin=log10(lower), ymax=log10(upper), 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="N x Tau (log10)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)


## N
to.plot. <- to.plot[to.plot$Parameter == "par_N",]
to.plot.$Linetype <- factor(to.plot.$Depth, levels = c("90", "270"))

# neutrally evolving samples
to.plot.. <- to.plot.[grepl("N", to.plot.$Paper_ID),]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=Median, ymin=lower, ymax=upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of stem cells (log10)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

# samples with unknown driver
to.plot.. <- to.plot.[grepl("U", to.plot.$Paper_ID) &
                        to.plot.$Known_CHIP.mutation == "no",]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=Median, ymin=lower, ymax=upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of stem cells (log10)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

# samples with known CH drivers
to.plot.. <- to.plot.[to.plot.$Known_CHIP.mutation == "yes",]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=Median, ymin=lower, ymax=upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of stem cells (log10)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)


## lambda

to.plot. <- to.plot[to.plot$Parameter == "par_lambda_ss",]
to.plot.$Linetype <- factor(to.plot.$Depth, levels = c("90", "270"))

# neutrally evolving samples
to.plot.. <- to.plot.[grepl("N", to.plot.$Paper_ID),]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Division rate (1/y)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

# samples with unknown driver
to.plot.. <- to.plot.[grepl("U", to.plot.$Paper_ID) &
                        to.plot.$Known_CHIP.mutation == "no",]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Division rate (1/y)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

# samples with known CH drivers
to.plot.. <- to.plot.[to.plot.$Known_CHIP.mutation == "yes",]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Division rate (1/y)") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

## mu

to.plot. <- to.plot[to.plot$Parameter == "par_mu",]
to.plot.$Linetype <- factor(to.plot.$Depth, levels = c("90", "270"))

# neutrally evolving samples
to.plot.. <- to.plot.[grepl("N", to.plot.$Paper_ID),]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=Median, ymin=lower, ymax=upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

# samples with unknown driver
to.plot.. <- to.plot.[grepl("U", to.plot.$Paper_ID) &
                        to.plot.$Known_CHIP.mutation == "no",]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=Median, ymin=lower, ymax=upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

# samples with known CH drivers
to.plot.. <- to.plot.[to.plot.$Known_CHIP.mutation == "yes",]
to.plot..$x <- as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2* 
  as.numeric(factor(to.plot..$Depth, levels = c(90, 270))) -1 

ggplot(to.plot..,
       aes(x=x, y=Median, ymin=lower, ymax=upper, 
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  expand_limits(x = 0, y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot..$Paper_ID, levels = unique(to.plot..$Paper_ID[order(to.plot..$Age)])))+0.2 - 1),
                     labels = (unique(to.plot..$Paper_ID)), name = "") + theme(aspect.ratio = 1)

dev.off()

############################################################################################################################################
## Fig. 6a-c, Extended Data Fig. 9a-d Plot the selection parameters (age and growth) estimated with 90x/270x from CD34 cells side-by-side

to.plot <- all.cd34.parameters[all.cd34.parameters$Parameter %in% c("age_of_clone", "age_of_clone_1", "age_of_clone_2",
                                                          "growth_per_year", "growth_per_year1", "growth_per_year2") &
                                 !grepl("N", all.cd34.parameters$Paper_ID),]
to.plot$Known_CHIP.mutation <- sample.info[as.character(to.plot$Paper_ID),]$CH.driver.found
to.plot$Leading.clone <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit
to.plot$Second.clone <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit.2
to.plot$Age <- floor(sample.info[to.plot$Paper_ID,]$Age)
to.plot <- to.plot[order(as.numeric(gsub("-.*", "", to.plot$Paper_ID))),]
to.plot <- to.plot[order(to.plot$Age),]
to.plot$Paper_ID <- factor(to.plot$Paper_ID,
                           levels=sort(unique(sample.info$Paper_ID[order(sample.info$Age)])))


pdf(paste0(analysis.directory, "/Figures/Figures_6abc_ExtDataFig_9.pdf"), width=3.5, height=2.5)

## age of 1st selected clone
to.plot. <- to.plot[to.plot$Parameter %in% c("age_of_clone", "age_of_clone_1"),]
to.plot.$Linetype <- factor(to.plot.$Depth, levels = c("90", "270"))
to.plot.$x <- as.numeric(factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(to.plot.$Age)])))+0.2* 
  as.numeric(factor(to.plot.$Depth, levels = c(90, 270))) -1 

ggplot(to.plot.,
       aes(x=x, y=Median, ymin=lower, ymax=upper, col = Leading.clone,
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of leading clone (years)") +
  expand_limits(y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(to.plot.$Age)])))+0.2 - 1),
                     labels = (unique(to.plot.$Paper_ID)), name = "") 

# summarize 270x where applicable in boxplots per driver
to.plot. <- to.plot.[order(to.plot.$Depth, decreasing = T),]
to.plot. <- to.plot.[-which(duplicated(to.plot.$Paper_ID)),]

ggplot(to.plot.,
       aes(x=Leading.clone, y=Median, ymin=lower, ymax=upper, col = Leading.clone)) + geom_boxplot(outliers = F) +
  geom_pointrange(fatten = 2, position = position_jitter()) + scale_color_manual(values=CHIP.color) +
  scale_y_continuous(name ="Age of leading clone (years)") +
  expand_limits(y = 0) + theme(legend.position = "bottom") 


## age of 2nd selected clone
to.plot. <- to.plot[to.plot$Parameter %in% c("age_of_clone_2"),]
to.plot.$x <- as.numeric(factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(to.plot.$Age)])))+0.2* 
  as.numeric(factor(to.plot.$Depth, levels = c(90, 270))) -1 

ggplot(to.plot.,
       aes(x=x, y=Median, ymin=lower, ymax=upper, col = Second.clone)) +
  geom_pointrange(fatten = 2) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Age of second clone (years)") +
  expand_limits(y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(to.plot.$Age)])))+0.2 - 1),
                     labels = (unique(to.plot.$Paper_ID)), name = "")

# summarize  in boxplots per driver

ggplot(to.plot.,
       aes(x=Second.clone, y=Median, ymin=lower, ymax=upper, col = Second.clone)) + geom_boxplot(outliers = F) +
  geom_pointrange(fatten = 2, position = position_jitter()) + scale_color_manual(values=CHIP.color) +
  scale_y_continuous(name ="Age of second clone (years)") +
  expand_limits(y = 0) + theme(legend.position = "bottom") 


## growth of 1st selected clone
to.plot. <- to.plot[to.plot$Parameter %in% c("growth_per_year", "growth_per_year1"),]
to.plot.$Linetype <- factor(to.plot.$Depth, levels = c("90", "270"))
to.plot.$x <- as.numeric(factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(to.plot.$Age)])))+0.2* 
  as.numeric(factor(to.plot.$Depth, levels = c(90, 270))) -1 

ggplot(to.plot.,
       aes(x=x, y=100*Median, ymin=100*lower, ymax=100*upper, col = Leading.clone,
           linetype=Linetype, shape = Linetype)) +
  geom_pointrange(fatten = 2) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Annual growth leading clone (%)") +
  expand_limits(y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(to.plot.$Age)])))+0.2 - 1),
                     labels = (unique(to.plot.$Paper_ID)), name = "")


## growth of 2nd selected clone
to.plot. <- to.plot[to.plot$Parameter %in% c("growth_per_year2"),]
to.plot.$x <- as.numeric(factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(to.plot.$Age)])))+0.2* 
  as.numeric(factor(to.plot.$Depth, levels = c(90, 270))) -1 

ggplot(to.plot.,
       aes(x=x, y=100*Median, ymin=100*lower, ymax=100*upper, col = Second.clone)) +
  geom_pointrange(fatten = 2) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Annual growth second clone (%)") +
  expand_limits( y = 0) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = unique(as.numeric(factor(to.plot.$Paper_ID, levels = unique(to.plot.$Paper_ID[order(to.plot.$Age)])))+0.2 - 1),
                     labels = (unique(to.plot.$Paper_ID)), name = "")



## growth of both selected clone
to.plot. <- to.plot[to.plot$Parameter %in% c("growth_per_year", "growth_per_year1", "growth_per_year2"),]
to.plot.$Linetype <- factor(to.plot.$Depth, levels = c("90", "270"))

# summarize 270x where applicable in boxplots per driver
to.plot. <- to.plot.[order(to.plot.$Depth, decreasing = T),]
to.plot.$Parameter[to.plot.$Parameter=="growth_per_year"] <- "growth_per_year1"
to.plot. <- to.plot.[-which(duplicated(paste(to.plot.$Paper_ID, to.plot.$Parameter))),]
to.plot.$Driver <- to.plot.$Leading.clone
to.plot.$Driver[to.plot.$Parameter == "growth_per_year2"] <- to.plot.$Second.clone[to.plot.$Parameter=="growth_per_year2"]

ggplot(to.plot.,
       aes(x=Driver, y=100*Median, ymin=100*lower, ymax=100*upper, col = Driver)) + geom_boxplot(outliers = F) +
  geom_pointrange(fatten = 2, position = position_jitter()) + scale_color_manual(values=CHIP.color) +
  scale_y_log10(name ="Annual growth per clone (%)") +
  expand_limits(y = 0) + theme(legend.position = "bottom") 


dev.off()

############################################################################################################################################
## Figure 5i-k: Compare physiological parameters (N, lambda, mu) between neutral evolution and clonal selection

to.plot <- all.cd34.parameters[all.cd34.parameters$Parameter %in% c("par_N", "par_mu", "par_lambda_ss"),]
# use 270x data wherever possible
to.plot <- to.plot[order(to.plot$Depth, decreasing = T),]
to.plot <- to.plot[-which(duplicated(paste(to.plot$Paper_ID, to.plot$Parameter))),]

pdf(paste0(analysis.directory, "/Figures/Figures_5i-k.pdf"), width=3.5, height=2.5)

to.plot$Selection <- ifelse(grepl("N", to.plot$Paper_ID),
                            "Neutral evolution", "Clonal selection")
to.plot$Selection <- factor(to.plot$Selection, levels = c("Neutral evolution",
                                                          "Clonal selection"))
to.plot$Age <- sample.info[to.plot$Paper_ID,]$Age
# make ages unique
to.plot$Age[duplicated(to.plot$Age)] <- to.plot$Age[duplicated(to.plot$Age)] +
  rnorm(1, mean = 0, sd = 0.5)

## N
to.plot. <- to.plot[to.plot$Parameter == "par_N",]

ggplot(to.plot.,
       aes(x=Selection, y=Median, ymin=lower, ymax=upper, col = Selection)) +
  geom_boxplot(outliers = F) +
  geom_pointrange(fatten = 2, position = position_jitter()) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of stem cells (log10)") +
  expand_limits(y = 0) + theme(legend.position = "bottom") +
  scale_color_manual(values = c("Neutral evolution" = "grey", "Clonal selection" = "red"))


## lambda
to.plot. <- to.plot[to.plot$Parameter == "par_lambda_ss",]

ggplot(to.plot.,
       aes(x=Selection, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col = Selection)) +
  geom_boxplot(outliers = F) +
  geom_pointrange(fatten = 2, position = position_jitter()) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Divisions per year") +
  expand_limits(y = 0) + theme(legend.position = "bottom") +
  scale_color_manual(values = c("Neutral evolution" = "grey", "Clonal selection" = "red"))

# against age
ggplot(to.plot.,
       aes(x=Age, y=365*10^Median, ymin=365*10^lower, ymax=365*10^upper, col = Selection)) +
  geom_pointrange(fatten = 2,) + 
  theme(axis.text.x=element_text(angle=90)) + scale_y_log10(name ="Divisions per year") +
  expand_limits(y = 0) + theme(legend.position = "bottom") + 
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("Neutral evolution" = "grey", "Clonal selection" = "red"))


## mu
to.plot. <- to.plot[to.plot$Parameter == "par_mu",]

ggplot(to.plot.,
       aes(x=Selection, y=Median, ymin=lower, ymax=upper, col = Selection)) +
  geom_boxplot(outliers = F) +
  geom_pointrange(fatten = 2, position = position_jitter()) + scale_color_manual(values=CHIP.color) +
  theme(axis.text.x=element_text(angle=90)) + scale_y_continuous(name ="Number of SSNVs per division") +
  expand_limits(y = 0) + theme(legend.position = "bottom") +
  scale_color_manual(values = c("Neutral evolution" = "grey", "Clonal selection" = "red"))

dev.off()


############################################################################################################################################
## Fig. 4f: Compare the estimated to the measured clone size for the CH cases with known drivers

## compare with driver vaf

to.plot <- all.cd34.parameters[all.cd34.parameters$Parameter %in% c("size_of_clone", "clone_size_1", "clone_size_2"),]
# use 270x data wherever possible
to.plot$Parameter[to.plot$Parameter=="size_of_clone"] <- "clone_size_1"
to.plot <- to.plot[order(to.plot$Depth, decreasing = T),]
to.plot <- to.plot[-which(duplicated(paste(to.plot$Paper_ID, to.plot$Parameter))),]
to.plot$CHIP.mutation <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit
to.plot$CHIP.mutation <- factor(to.plot$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))
to.plot$CHIP.mutation.2 <- sample.info[to.plot$Paper_ID,]$CHIP.mutation.associated.with.fit.2
to.plot$CHIP.mutation.2 <- factor(to.plot$CHIP.mutation.2, levels=c("ASXL1", "DNMT3A", "TET2", "unknown driver", "no selected clone"))


to.plot <- to.plot[(to.plot$CHIP.mutation!="no selected clone" & !is.na(to.plot$CHIP.mutation) &
                      to.plot$CHIP.mutation != "unknown driver") |
                     (to.plot$CHIP.mutation.2 != "unknown driver" & !is.na(to.plot$CHIP.mutation.2)),]

to.plot$Driver.mean <- apply(to.plot, 1, function(x){
  
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  driver.information <- driver.information[which.max(driver.information$VAF),]
  vaf=as.numeric(driver.information["VAF"])
  if(is.na(vaf)){
    vaf <- 0
  }
  return(vaf)
})
to.plot$Driver.min <- apply(to.plot, 1, function(x){
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  
  driver.information <- driver.information[which.max(driver.information$VAF),]
  vaf=as.numeric(driver.information["VAF"])
  DP=as.numeric(driver.information["Depth"])
  print(vaf)
  if(length(vaf)==0 || is.na(vaf)){
    vaf <- 0
    DP <- 0
  }
  lower=vaf-1.96*sqrt(vaf*(1-vaf)/DP) ## Wald approximation for 95% binomial CI
  return(lower)
})
to.plot$Driver.max <- apply(to.plot, 1, function(x){
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  
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

# for case 17-T there's a second TET2 mutation in the same clone; also show their VAF

tmp <- to.plot[to.plot$Paper_ID=="17-T",]
tmp$Driver.mean <- apply(tmp, 1, function(x){
  
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  driver.information <- driver.information[order(driver.information$VAF, decreasing = T)[2],]
  vaf=as.numeric(driver.information["VAF"])
  if(is.na(vaf)){
    vaf <- 0
  }
  return(vaf)
})
tmp$Driver.min <- apply(tmp, 1, function(x){
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  
  driver.information <- driver.information[order(driver.information$VAF, decreasing = T)[2],]
  vaf=as.numeric(driver.information["VAF"])
  DP=as.numeric(driver.information["Depth"])
  print(vaf)
  if(length(vaf)==0 || is.na(vaf)){
    vaf <- 0
    DP <- 0
  }
  lower=vaf-1.96*sqrt(vaf*(1-vaf)/DP) ## Wald approximation for 95% binomial CI
  return(lower)
})
tmp$Driver.max <- apply(tmp, 1, function(x){
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  
  driver.information <- driver.information[order(driver.information$VAF, decreasing = T)[2],]
  vaf=as.numeric(driver.information["VAF"])
  DP=as.numeric(driver.information["Depth"])
  
  if(length(vaf)==0 || is.na(vaf)){
    vaf <- 0
    DP <- 0
  }
  lower=vaf+1.96*sqrt(vaf*(1-vaf)/DP) ## Wald approximation for 95% binomial CI
  return(lower)
})

to.plot <- rbind(to.plot, tmp)
rm(tmp)


# 2nd driver

to.plot$Driver.mean.2 <- apply(to.plot, 1, function(x){
  
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit.2"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  driver.information <- driver.information[which.max(driver.information$VAF),]
  vaf=as.numeric(driver.information["VAF"])
  if(is.na(vaf)){
    vaf <- 0
  }
  return(vaf)
})
to.plot$Driver.min.2 <- apply(to.plot, 1, function(x){
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit.2"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  
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
to.plot$Driver.max.2 <- apply(to.plot, 1, function(x){
  column <- colnames(putative.drivers)[gsub("_.*", "", colnames(putative.drivers))==x["Paper_ID"]]
  column <- column[grepl("CD34", column)]
  if(any(grepl("deep", column))){
    column <- column[grepl("deep", column)]
  }
  
  driver.information <- putative.drivers[,c("CHROM", "POS", "FORMAT", "REF", "ALT", "GENE", "AAchange",
                                            column)]
  driver.information <- driver.information[!is.na(driver.information[,column]),]
  driver.information$VAF <- Extract.info.from.vcf(driver.information, info="VAF", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information$Depth <- Extract.info.from.vcf(driver.information, info="depth", mutationcaller = "mpileup", sample.col.mpileup = column)
  driver.information <- driver.information[driver.information$VAF>0 & driver.information$GENE==as.character(sample.info[as.character(x["Paper_ID"]),"CHIP.mutation.associated.with.fit.2"]),,drop=F]
  driver.information <- driver.information[driver.information$VAF < 0.9,] # exclude germline SNPs
  
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
to.plot$Driver.mean.2[to.plot$Driver.mean.2==0] <- NA


pdf(paste0(analysis.directory, "/Figures/Figures_4f.pdf"), width=4, height=4)

ggplot(to.plot[to.plot$Driver.mean!=0 & to.plot$Parameter %in% c("size_of_clone", "clone_size_1"),], aes(x=Median/2, xmin = lower/2, xmax = upper/2, color=CHIP.mutation,
                                                                                                        y=Driver.mean, ymin = Driver.min, ymax=Driver.max)) + 
  geom_point(data = to.plot[to.plot$Driver.mean!=0 & to.plot$Parameter %in% c("size_of_clone", "clone_size_1"),]) +
  geom_errorbar() + geom_errorbarh() + scale_color_manual(values=CHIP.color) +
  geom_abline(slope = 1, intercept = 0, linetype=2)  + scale_x_continuous("VAF population genetics model") +
  scale_y_continuous("VAF WGS") + theme(aspect.ratio = 1) + 
  geom_point(data = to.plot[to.plot$Driver.mean.2!=0 & to.plot$Parameter %in% "clone_size_2",],
             aes(x=Median/2,  color=CHIP.mutation.2,
                 y=Driver.mean.2), inherit.aes = F) + 
  geom_errorbar(data = to.plot[to.plot$Driver.mean.2!=0 & to.plot$Parameter %in% "clone_size_2",],
                aes(x=Median/2, xmin = lower/2, xmax = upper/2, color=CHIP.mutation.2,
                    y=Driver.mean.2, ymin = Driver.min.2, ymax=Driver.max.2), width = 0, inherit.aes = F) +
  geom_errorbarh(data = to.plot[to.plot$Driver.mean.2!=0 & to.plot$Parameter %in% "true_size_2",],
                 aes(xmin = lower/2, xmax = upper/2, color=CHIP.mutation.2,
                     y=Driver.mean.2), inherit.aes = F, height = 0)

dev.off()

############################################################################################################################################
## Fig. 6d: Plot incidence of selected clone emergence across age

to.plot.driver <- all.cd34.parameters[all.cd34.parameters$Parameter %in% c("par_ts1", "par_ts2", "par_t_s_absolute"),]
# use 270x data wherever possible
to.plot.driver$Parameter[to.plot.driver$Parameter=="par_t_s_absolute"] <- "par_ts1"
to.plot.driver <- to.plot.driver[order(to.plot.driver$Depth, decreasing = T),]
to.plot.driver <- to.plot.driver[-which(duplicated(paste(to.plot.driver$Paper_ID, to.plot.driver$Parameter))),]


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
#to.plot.driver$CHIP.mutation <- sample.info[to.plot.driver$Sample,]$CHIP.mutation.associated.with.fit
to.plot.driver$CHIP.mutation <- factor(to.plot.driver$CHIP.mutation, levels=c("ASXL1", "DNMT3A", "TET2", "unknown CHIP mutation", "healthy donor"))
to.plot.driver$Paper_ID <- factor(to.plot.driver$Paper_ID,
                               levels=sort(unique(sample.info$Paper_ID)))

## cohort summary
to.plot. <- data.frame(Age = seq(0, ceiling(max(to.plot.driver$upper)/365/25)*25),
                       Incidence = sapply(seq(0,ceiling(max(to.plot.driver$upper)/365/25)*25), function(a){sum(to.plot.driver$Median/365<=a)})/nrow(to.plot.driver),
                       Lower = sapply(seq(0, ceiling(max(to.plot.driver$upper)/365/25)*25), function(a){sum(to.plot.driver$upper/365<=a)})/nrow(to.plot.driver),
                       Upper = sapply(seq(0,ceiling(max(to.plot.driver$upper)/365/25)*25), function(a){sum(to.plot.driver$lower/365<=a)})/nrow(to.plot.driver))


pdf(paste0(analysis.directory, "/Figures/Figures_6d.pdf"), width=3.5, height=2.5)

ggplot(to.plot., aes(x=Age, y=Incidence, ymin=Lower, ymax=Upper)) + geom_ribbon(fill="grey", alpha=0.5) + geom_line() +
  geom_point(data=to.plot.driver, aes(x=Median/365, y=y, col=CHIP.mutation), inherit.aes = F) +
  geom_errorbarh(data=to.plot.driver, aes(xmin=lower/365, xmax=upper/365, col=CHIP.mutation, y=y), inherit.aes = F, height=0) +
  scale_color_manual(values=CHIP.color) + scale_x_continuous(breaks = seq(0, 75, 25)) + theme(aspect.ratio = 1)

dev.off()


############################################################################################################################################
## Extended Data Fig. 7g: Compare the estimated parameters between one-clone and two-clone model

cases.w.2.clones <- colnames(model.support.selection.2_clones)[model.support.selection.2_clones["Clone 2",]>= 15]
cases.w.2.clones <- gsub("_.*", "", cases.w.2.clones)
to.plot.1cm <- selected.parameters.1cm[selected.parameters.1cm$Paper_ID %in% cases.w.2.clones &
                                         selected.parameters.1cm$Sort=="CD34_deep",]
to.plot.2cm <- selected.parameters.2cm.2clones[selected.parameters.2cm.2clones$Paper_ID %in% cases.w.2.clones &
                                                 model.support.selection.2_clones["Clone 2", paste(selected.parameters.2cm.2clones$Paper_ID,
                                                                                                   selected.parameters.2cm.2clones$Mode, sep="_")]>=15,]
## posterior probability of 2nd clone
to.plot.2cm$Posterior.2 <- model.support.selection.2_clones["Clone 2",paste(to.plot.2cm$Paper_ID, to.plot.2cm$Mode, sep="_")]


# if both topologies fit, take parameters associated with higher posterior; else take average

to.plot.2cm. <- data.frame()

for(i in unique(to.plot.2cm$Paper_ID)){
  for(j in unique(to.plot.2cm$Parameter)){
    tmp <- to.plot.2cm[to.plot.2cm$Paper_ID == i &
                         to.plot.2cm$Parameter==j,]
    if(nrow(tmp)==0){next}
    if(nrow(tmp)==2){
      if(tmp$Posterior.2[1]== tmp$Posterior.2[2]){
        to.plot.2cm. <- rbind(to.plot.2cm.,
                                     data.frame(lower = mean(tmp$lower),
                                                upper = mean(tmp$upper),
                                                Parameter = j,
                                                Median = mean(tmp$Median),
                                                Paper_ID = i))
      }else{
        to.plot.2cm. <- rbind(to.plot.2cm.,
                                     data.frame(lower = tmp$lower[which.max(tmp$Posterior.2)],
                                                upper = tmp$upper[which.max(tmp$Posterior.2)],
                                                Parameter = j,
                                                Median = tmp$Median[which.max(tmp$Posterior.2)],
                                                Paper_ID = i))
      }
    }else{
      to.plot.2cm. <- rbind(to.plot.2cm.,
                                   data.frame(lower = mean(tmp$lower),
                                              upper = mean(tmp$upper),
                                              Parameter = j,
                                              Median = mean(tmp$Median),
                                              Paper_ID = i))
    }
    
  }
}

to.plot.2cm <- to.plot.2cm.
rm(to.plot.2cm.)

to.plot <- data.frame(Paper_ID = to.plot.1cm$Paper_ID,
                      Parameter = to.plot.1cm$Parameter,
                      One_clone.median = to.plot.1cm$Median,
                      One_clone.min = to.plot.1cm$lower,
                      One_clone.max = to.plot.1cm$upper)

to.plot$Parameter <- replace(to.plot$Parameter, to.plot$Parameter=="growth_per_year", "growth_per_year1")
to.plot$Parameter <- replace(to.plot$Parameter, to.plot$Parameter=="age_of_clone", "age_of_clone_1")
to.plot <- to.plot[to.plot$Parameter %in% c("par_N", "par_lambda_ss",
                                            "par_mu", "age_of_clone_1", "growth_per_year1"),]

to.plot$Two_clone.median <- apply(to.plot, 1, function(x){
  to.plot.2cm[to.plot.2cm$Paper_ID==x["Paper_ID"] &
                to.plot.2cm$Parameter==x["Parameter"],"Median"]
})

to.plot$Two_clone.min <- apply(to.plot, 1, function(x){
  to.plot.2cm[to.plot.2cm$Paper_ID==x["Paper_ID"] &
                to.plot.2cm$Parameter==x["Parameter"],"lower"]
})

to.plot$Two_clone.max <- apply(to.plot, 1, function(x){
  to.plot.2cm[to.plot.2cm$Paper_ID==x["Paper_ID"] &
                to.plot.2cm$Parameter==x["Parameter"],"upper"]
})


pdf(paste0(analysis.directory, "/Figures/Ext_data_Fig_7g.pdf"), width=3.5, height=3.5)

ggplot(to.plot[to.plot$Parameter=="par_N",], 
       aes(x = One_clone.median, xmin = One_clone.min, xmax = One_clone.max,
           y = Two_clone.median, ymin = Two_clone.min, ymax = Two_clone.max)) +
  geom_point() + geom_errorbar(width = 0) + geom_errorbarh(height = 0) +
  scale_x_continuous(name = "One-clone model") +
  scale_y_continuous(name = "Two-clone model") + ggtitle("Number of stem cells (log10)") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + expand_limits(x=0, y=0)

ggplot(to.plot[to.plot$Parameter=="par_lambda_ss",], 
       aes(x = 365*10^One_clone.median, xmin = 365*10^One_clone.min, xmax = 365*10^One_clone.max,
           y = 365*10^Two_clone.median, ymin = 365*10^Two_clone.min, ymax = 365*10^Two_clone.max)) +
  geom_point() + geom_errorbar(width = 0) + geom_errorbarh(height = 0) +
  scale_x_log10(name = "One-clone model") +
  scale_y_log10(name = "Two-clone model") + ggtitle("Division rate (1/y)")+
  geom_abline(slope = 1, intercept = 0, linetype = 2) + expand_limits(x=0, y=0)

ggplot(to.plot[to.plot$Parameter=="par_mu",], 
       aes(x = One_clone.median, xmin = One_clone.min, xmax = One_clone.max,
           y = Two_clone.median, ymin = Two_clone.min, ymax = Two_clone.max)) +
  geom_point() + geom_errorbar(width = 0) + geom_errorbarh(height = 0) +
  scale_x_continuous(name = "One-clone model") +
  scale_y_continuous(name = "Two-clone model") + ggtitle("SSNVs per division")+
  geom_abline(slope = 1, intercept = 0, linetype = 2) + expand_limits(x=0, y=0)

ggplot(to.plot[to.plot$Parameter=="age_of_clone_1",], 
       aes(x = One_clone.median, xmin = One_clone.min, xmax = One_clone.max,
           y = Two_clone.median, ymin = Two_clone.min, ymax = Two_clone.max)) +
  geom_point() + geom_errorbar(width = 0) + geom_errorbarh(height = 0) +
  scale_x_continuous(name = "One-clone model") +
  scale_y_continuous(name = "Two-clone model") + ggtitle("Age of leading clone (years)")+
  geom_abline(slope = 1, intercept = 0, linetype = 2) + expand_limits(x=0, y=0)

ggplot(to.plot[to.plot$Parameter=="growth_per_year1",], 
       aes(x = 100*One_clone.median, xmin = 100*One_clone.min, xmax = 100*One_clone.max,
           y = 100*Two_clone.median, ymin = 100*Two_clone.min, ymax = 100*Two_clone.max)) +
  geom_point() + geom_errorbar(width = 0) + geom_errorbarh(height = 0) +
  scale_x_log10(name = "One-clone model") +
  scale_y_log10(name = "Two-clone model") + ggtitle("Annual growth leading clone (%)")+
  geom_abline(slope = 1, intercept = 0, linetype = 2) + expand_limits(x=0, y=0)


dev.off()

####################################################################################################################################################
## In Extended Data Fig. 7f we compare parameter estimates obtained with the one-clone model assuming size compensation for the selected clone and the one-clone model without size-compensation
## To load the parameters obtained without size compensation, we source the file "Assess_fits_heme_WGS_data_no_size_compensation.R"

source(paste0(custom.script.directory, "/Analysis_and_figures/Assess_fits_heme_WGS_data_no_size_compensation.R"))

## parameters obtained without size compensation
tmp.1 <- selected.parameters.1cm.nsc
tmp.1$Size_compensation <- "no"

## parameters obtained with size compensation
tmp.2 <- selected.parameters.1cm
tmp.2$Size_compensation <- "yes"

to.plot <- rbind(tmp.2[tmp.2$Paper_ID %in% tmp.1$Paper_ID & tmp.2$Sort == "CD34" & tmp.2$Depth %in% c(90, 120),],
                 tmp.1)

rm(tmp.1, tmp.2)


pdf(paste0(analysis.directory, "/Figures/Ext_Data_Fig_7f.pdf"), width=3.5, height=2.5)


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




