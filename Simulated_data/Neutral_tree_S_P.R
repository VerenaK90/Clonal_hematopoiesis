###########################################################################################################################################
## Simulate a neutrally evolving tree with a progenitor compartment

library(SCIFER)
library(doParallel)
library(foreach)
library(parallel)

### simulate trees with a stem cell count of 1000 that divide once  a year and acquire 1 mutation per division.
### simulate progenitor cells that divide 5 times a year and differentiate 4.6 times a year. Hence there are 2.5x as many progenitors as stem cells
### the simulated trees evolve neutrally, hence, set the driver mutation rate to zero.

N = 1000
NP = 2500
mut.rate = 1
time.max = 75
parms.exp=c(lambda.s=1, delta.s=0, lambda.p=1, delta.p=0, alpha.s=1)
parms.steady=c(lambda.s=1, delta.s=0, alpha.s=1, lambda.p=4.6, delta.p=5.0)
driver.mode="fixed_time"
t.driver = 100
mut.rate.D = 0
mutation.mode="constant" # model neutral mutation acquisition at constant rate
s.shape = 1.5
s.rate = 70

## parallelize the runs
n.cores <- 7
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)


trees <- foreach(sim.nr = 1:10) %dopar% {
  set.seed(sim.nr*354)
  gillespie.sim.s.p(parms.exp, parms.steady, time.max=75, time.samples=seq(0, 75, 25), N=N, NP=NP,
                    mut.rate = mut.rate, mutation.mode="Binomial", driver.mode = driver.mode, t.driver = t.driver, 
                    mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, 
                    tau = N/10, report.at.f = seq(0.05, 1, 0.05))
}
save(trees, file=paste0("RData/Simulated_data/Neutral_tree_1000_HSCs_2500_Ps.RData"))


###########################################################################################################################################
################ Plot the simulated trees, the corresponding VAF histograms and compute the distance between the cumulative VAF distribution computed from all cells to stem cells

tree.plots <- list()
vaf.plots <- list()
distance.to.hsc <- data.frame()

## Iterate through the simulations
for(i in 1:length(trees)){

  ## For each instance we reported 4 time points
  
  ## simulate VAF for the full trees..
  vafs.all <- lapply(trees[[i]], function(x){
    vafs <- get_vaf_from_tree(x)
    vafs
  })
  ## ... the stem cell trees only ...
  hsc.trees <- lapply(trees[[i]], function(x){
    .subset.cell.type(x, cell.type = 1)
  })
  vafs.hsc <- lapply(hsc.trees, function(x){
    vafs <- get_vaf_from_tree(x)
    vafs
  })
  ## ... and the progenitor trees only!
  p.trees <- lapply(trees[[i]], function(x){
    .subset.cell.type(x, cell.type = 2)
  })
  vafs.p <- lapply(p.trees, function(x){
    vafs <- get_vaf_from_tree(x)
    vafs
  })
  
  save(vafs.all, vafs.hsc, vafs.p, file=paste0("./RData/Simulated_data/True_VAFs_neutral_1000_HSCs_2500_Ps_simnr_", i, ".RData")) ## pure VAFs, no sequencing simulated on top
  
  for(j in 1:length(trees[[i]])){
    print(j)
    ## full tree
    tree. <- trees[[i]][[j]]
    ## annotate tips
    meta.data <- data.frame(node = as.numeric(tree.$tip.label), Cell.type = as.character(tree.$tip.class))
    tree.plots[[length(tree.plots)+1]] <- ggtree(tree.) + ggtitle(paste("Sim. nr. ", i, "; age = ", names(trees[[i]])[j], "; all cells")) 
    tree.plots[[length(tree.plots)]] <- tree.plots[[length(tree.plots)]] %<+% meta.data + 
      geom_tippoint(aes(col=Cell.type), size = 0.5) + scale_color_manual(values = c("1" = "purple", "2" = "orange"))
    ## sample tree down to 10%
    subsampled.tree <- .simulate.sampling(trees[[i]][[j]], sample.size = 350)
    meta.data <- data.frame(node = as.numeric(subsampled.tree$tip.label), Cell.type = as.character(subsampled.tree$tip.class))
    tree.plots[[length(tree.plots)+1]] <- ggtree(subsampled.tree) + ggtitle(paste("Sim. nr. ", i, "; age = ", names(trees[[i]])[j], "; 10% sampled cells")) 
    tree.plots[[length(tree.plots)]] <- tree.plots[[length(tree.plots)]] %<+% meta.data + 
      geom_tippoint(aes(col=Cell.type), size = 0.5) + scale_color_manual(values = c("1" = "purple", "2" = "orange"))
    
    ## stem cells only
    hsc.tree <- hsc.trees[[j]]
    ## progenitors only
    p.tree <- p.trees[[j]]
    tree.plots[[length(tree.plots)+1]] <- ggtree(hsc.tree) + ggtitle(paste("Sim. nr. ", i, "; age = ", names(trees[[i]])[j], "; HSCs"))
    tree.plots[[length(tree.plots)+1]] <- ggtree(p.tree) + ggtitle(paste("Sim. nr. ", i, "; age = ", names(trees[[i]])[j], "; Progenitors"))
    
    ## plot the VAF distribution for the whole tree, stem cells only, progenitors only
    to.plot <- rbind(data.frame(VAF=seq(0.01, 1, 0.01),
                                M = sapply(seq(0.01, 1, 0.01), function(x){
                                  sum(vafs.all[[j]] >= x)}),
                                Cells = "HSC + P"
    ),
    data.frame(VAF=seq(0.01, 1, 0.01),
               M = sapply(seq(0.01, 1, 0.01), function(x){
                 sum(vafs.hsc[[j]] >= x)}),
               Cells = "HSC"
    ),
    data.frame(VAF=seq(0.01, 1, 0.01),
               M = sapply(seq(0.01, 1, 0.01), function(x){
                 sum(vafs.p[[j]] >= x)}),
               Cells = "P"
    ))
    
    vaf.plots[[length(vaf.plots) + 1]] <- ggplot(to.plot, aes(x=1/VAF, y=M, col=Cells)) + geom_point(size = 0.5)+
      ggtitle(paste("Sim. nr. ", i, "; age = ", names(trees[[i]])[j]))
    vaf.plots[[length(vaf.plots) + 1]] <- ggplot(to.plot[to.plot$VAF>=0.05,], aes(x=1/VAF, y=M, col=Cells)) + geom_point()+
      ggtitle(paste("Sim. nr. ", i, "; age = ", names(trees[[i]])[j]))
    
    ## compute the distance between the VAF distribution measured by the full sample and by sorting for Ps to the VAF distribution in HSCs
    ## take the absolute distance at each evaluated VAF and normalize by the number of variants in HSCs with this VAF.
    ## To avoid infinite values due to zero counts, add a pseudocount of 1
    to.plot$M <- to.plot$M + 1
    distance.to.hsc <- rbind(distance.to.hsc,
                             data.frame(Sample = "HSC + P",
                                        Distance = sum(abs(to.plot$M[to.plot$Cells=="HSC + P" & to.plot$VAF >= 0.05] -
                                                             to.plot$M[to.plot$Cells=="HSC" & to.plot$VAF >= 0.05])/
                                                         to.plot$M[to.plot$Cells=="HSC" & to.plot$VAF >= 0.05], na.rm = T)/length(unique(to.plot$VAF[to.plot$VAF>=0.05])),
                                        Resolution = "5%",
                                        Age = names(trees[[i]])[j]),
                             data.frame(Sample = "HSC + P",
                                        Distance = sum(abs(to.plot$M[to.plot$Cells=="HSC + P"] -
                                                             to.plot$M[to.plot$Cells=="HSC"])/
                                                         to.plot$M[to.plot$Cells=="HSC"], na.rm = T)/length(unique(to.plot$VAF)),
                                        Resolution = "1%",
                                        Age = names(trees[[i]])[j]),
                             data.frame(Sample = "P",
                                        Distance = sum(abs(to.plot$M[to.plot$Cells=="P" & to.plot$VAF >= 0.05] -
                                                             to.plot$M[to.plot$Cells=="HSC" & to.plot$VAF >= 0.05])/
                                                         to.plot$M[to.plot$Cells=="HSC" & to.plot$VAF >= 0.05], na.rm = T)/length(unique(to.plot$VAF[to.plot$VAF>=0.05])),
                                        Resolution = "5%",
                                        Age = names(trees[[i]])[j]),
                             data.frame(Sample = "P",
                                        Distance = sum(abs(to.plot$M[to.plot$Cells=="P"]-
                                                             to.plot$M[to.plot$Cells=="HSC"])/
                                                         to.plot$M[to.plot$Cells=="HSC"], na.rm = T)/length(unique(to.plot$VAF)),
                                        Resolution = "1%",
                                        Age = names(trees[[i]])[j]))
  }
}

pdf("./Neutral_S_P_trees.pdf", width = 5, height = 4)
print(tree.plots)
dev.off()

pdf("./Neutral_S_P_VAF_distr.pdf", width = 6, height = 4)
ggarrange(plotlist=vaf.plots, nrow = 2, ncol=4, common.legend = T)
dev.off()

pdf("./Neutral_difference_P_S.pdf", width = 6, height = 4)

ggplot(distance.to.hsc[distance.to.hsc$Resolution=="5%",], aes(x=Age, y=Distance, col=Sample)) + geom_point() +
  stat_summary(aes(y = Distance,group=Sample, col=Sample), fun.y=mean, geom="line") +
  scale_y_continuous(name="Mean relative distance") + ggtitle("Resolution = 5%")

ggplot(distance.to.hsc[distance.to.hsc$Resolution=="1%",], aes(x=Age, y=Distance, col=Sample)) + geom_point() +
  stat_summary(aes(y = Distance,group=Sample, col=Sample), fun.y=mean, geom="line") +
  scale_y_continuous(name="Mean relative distance") + ggtitle("Resolution = 1%")

dev.off()

