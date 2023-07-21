## use the data from Fabre et al. and generate pseudo-bulks

library(ggplot2)
library(ggpubr)
library(phangorn)
library(RRphylo)
library(cgwtools)
library(phytools)

source("Simulated_data/Tree_post_processing.R") ## source modalities to extract information from trees
folders <- list.files("./Fabre_et_al/", pattern = "id")

driver.colors <- c("SF3B1" = "orange", "CBL" = "purple", CTCF = "purple", TET2 = "firebrick", U2AF1 = "green",
                   SRSF2 = "firebrick", PPM1D = "firebrick", TP53 = "firebrick") 
snvs <- list()

for(i in folders){
  
  tree.file <- list.files(paste0("./Fabre_et_al/", i), pattern = "tree", full.names = T)
  mut.file <- list.files(paste0("./Fabre_et_al/", i), pattern = "details", full.names = T)
  
  load(tree.file)
  load(mut.file)
  
  muts <- get(lsdata(mut.file))
  ## subset the drivers reported in Fabre et al.:
  drivers <- muts[muts$coding_change_CHgeneInModel=="Coding change" & muts$Gene %in% c("SF3B1", "CBL", "TET2", "CTCF", "U2AF1", "SRSF2", "PPM1D", "TP53"),]
  
  ## compute branch length in actual # mutations using mutation file
  adjusted.length <- rep(NA, nrow(tree_SNV_c_ultra$edge))
  
  for(j in 1:nrow(tree_SNV_c_ultra$edge)){
    adjusted.length[j] <- sum(muts$node==tree_SNV_c_ultra$edge[j,2])
  }
  
  tree_SNV_c_ultra$edge.length <- adjusted.length
  tree_SNV_c_ultra$edge.color <- rep("black", length(tree_SNV_c_ultra$edge.length))
  # order drivers by hierarchy
  drivers <- drivers[c(which(drivers$node >= tree$Nnode)[order(drivers$node[drivers$node >= tree$Nnode])], 
                       which(drivers$node < tree$Nnode)[order(drivers$node[drivers$node < tree$Nnode])]), ]
  
  for(j in 1:nrow(drivers)){
    tree_SNV_c_ultra$edge.color[tree$edge[,2] == drivers[j,]$node] <- driver.colors[drivers[j,]$Gene]
    children <- getDescendants(tree_SNV_c_ultra, node = drivers[j,]$node) ## also color downstream of the node
    for(k in children){
      tree_SNV_c_ultra$edge.color[tree$edge[,2] == k] <- driver.colors[drivers[j,]$Gene]
    }
  }
  
  ## get vafs from tree
  vafs <- get_vaf_from_tree(tree_SNV_c_ultra)
  
  
  p <- list()
  
  p[[1]] <- ggplot(data=data.frame(VAF=vafs[vafs >= 0.01]), aes(x=VAF)) + geom_histogram(binwidth = 0.01) + 
    scale_x_continuous(limits = c(0,1), name="VAF") + scale_y_continuous(name="# SSNVs") 
 
  # cumulative distribution:
  to.plot <- data.frame(VAF=seq(0.01, 1, 0.01), M=sapply(seq(0.01, 1, 0.01), function(x){
    sum(vafs >=x)
  }))
  
  p[[2]] <- ggplot(data=to.plot, aes(x=1/VAF, y=M)) + geom_point() + geom_line() +
    scale_x_continuous(limits = c(1,1/min(to.plot$VAF)), name="VAF") + scale_y_continuous(name="# SSNVs")+
    ggtitle("True VAFs")

  if(!dir.exists(paste0("./Fabre_et_al/", i))){
    dir.create(paste0("./Fabre_et_al/", i))
  }
  
  
  pdf(paste0("./Fabre_et_al/", i, "/VAFs.pdf"), width=6, height = 6)
  
  print(ggarrange(plotlist=p, nrow=2, ncol=2))
  
  dev.off()
  
  snvs[[i]] <- vafs
  
  ## Read in the tree file for this tumor and plot

  tree <- tree_SNV_c_ultra
  
  pdf(paste0("./Fabre_et_al/", i, "/Trees.pdf"), width=6, height = 6)
  
  plot.phylo(tree, use.edge.length = TRUE, show.tip.label = F, direction = "downwards", edge.color = tree$edge.color)
  axisPhylo(side=4, backward = F)
 
  dev.off()
  
  ## length of the tree
  length.of.tree <- c(length.of.tree, max(node.depth.edgelength(tree)))
  
  ntips <- length(tree$tip.label)
 
}

save(snvs, file="./RData/Fabre_et_al/SNVs.RData")




