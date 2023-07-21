## use the data from Mitchell et al. and generate pseudo-bulks

library(ggplot2)
library(ggpubr)
library(phangorn)
library(RRphylo)

folders <- list.files("./Published_data/Mitchell_et_al/", pattern = "00")

folders <- setdiff(folders, "KX007") # no data available

snvs <- list()

for(i in folders){
  
  mut.file <- list.files(paste0("./Published_data/Mitchell_et_al/", i),
                         pattern = "annotated", full.names = T)
  
  load(mut.file)
  ## save mutations as vcf file for input for annovar to search for germline variants with gnomad
  vcf <- data.frame("CHROM" = sapply(rownames(filtered_muts$Genotype_shared_bin), function(x){
    as.numeric(strsplit(x, split="-")[[1]][1])
  }), POS = sapply(rownames(filtered_muts$Genotype_shared_bin), function(x){
    as.numeric(strsplit(x, split="-")[[1]][2])
  }), ID=".",
  REF = sapply(rownames(filtered_muts$Genotype_shared_bin), function(x){
    strsplit(x, split="-")[[1]][3]
  }), ALT = sapply(rownames(filtered_muts$Genotype_shared_bin), function(x){
    strsplit(x, split="-")[[1]][4]
  }),
  QUAL=10,
  FILTER="PASS",
  INFO="ID='.'")
  
  colnames(vcf)[1] <- "#CHROM"
  
  ## order by position
  vcf <- vcf[order(vcf$POS),]
  vcf <- vcf[order(vcf$`#CHROM`),]
  
  ## to match input format required by vcfanno
  write(x = "##fileformat=VCFv4.2", file=paste0("./Published_data/Mitchell_et_al/", i,
                               "/SNVs_", i, ".vcf"))
  write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">"', file=paste0("./Published_data/Mitchell_et_al/", i,
                                                  "/SNVs_", i, ".vcf"), append=T)
  write.table(vcf, file=paste0("./Published_data/Mitchell_et_al/", i,
                               "/SNVs_", i, ".vcf"), quote = F, row.names = F, sep="\t", append = T)

  vafs <- rowSums(filtered_muts$Genotype_shared_bin!=0)/ncol(filtered_muts$Genotype_shared_bin!=0)/2
  
  p <- list()

  p[[1]] <- ggplot(data=data.frame(VAF=vafs[vafs>=0.01]), aes(x=VAF)) + geom_histogram(binwidth = 0.01) + 
    scale_x_continuous(limits = c(0,1), name="VAF") + scale_y_continuous(name="# SSNVs") 
 
  to.plot <- data.frame(VAF=seq(0.01, 1, 0.01), M=sapply(seq(0.01, 1, 0.01), function(x){
    sum(vafs >=x)
  }))
  
  p[[2]] <- ggplot(data=to.plot, aes(x=1/VAF, y=M)) + geom_point() + geom_line() +
    scale_x_continuous(limits = c(1,1/min(to.plot$VAF)), name="VAF") + scale_y_continuous(name="# SSNVs")
  
 
  pdf(paste0("./Published_data/Mitchell_et_al/", i, "/VAFs.pdf"),
      width=6, height = 6)
  
  print(ggarrange(plotlist=p, nrow=2, ncol=2))
  
  dev.off()
  
  snvs[[i]] <- vafs
  
  sim.snvs[[i]] <- sim.vafs
  
  ## Read in the tree file for this tumor and plot
  
  tree.file <- list.files(paste0("./Published_data/Mitchell_et_al/", i),
                         pattern = "tree", full.names = T)
  
  tree <- read.tree(tree.file)
  
  pdf(paste0("./Published_data/Mitchell_et_al/", i, "/Trees.pdf"),
      width=6, height = 6)
  

  plot.phylo(tree, use.edge.length = TRUE, show.tip.label = F, direction = "downwards")
  axisPhylo(side=4, backward = F)
  
  dev.off()
  
}

save(snvs, sim.snvs, file="./RData/Mitchell_et_al/SNVs.RData")
