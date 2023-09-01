library(ggVennDiagram)
library(ggplot2); theme_set(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=8, color="black"),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.75),
                                  axis.ticks = element_line(color = "black", size=0.75),
                                  axis.text = element_text(size=8, color="black")))
############################################################################################################################################
## Read in data
lee.six.mutect.strelka <- read.delim("./Lee-Six_et_al/Mutect_Strelka/Mut_table_mutect_strelka.csv", sep=",")
identifiers.mutect.strelka <- paste(lee.six.mutect.strelka$chr, lee.six.mutect.strelka$pos,  sep=".")
lee.six.caveman <- read.delim("./Lee-Six_et_al/Caveman/Shearwater_calls_FDR0.95_all_muts.txt", sep="\t")
lee.six.caveman <- lee.six.caveman[,-ncol(lee.six.caveman)]
## subset on SNVs
lee.six.caveman <- lee.six.caveman[lee.six.caveman$ALT %in% c("A", "C", "G", "T") &
                                     lee.six.caveman$REF %in% c("A", "C", "G", "T"),]
identifiers.caveman <- paste(lee.six.caveman$X.CHROM, lee.six.caveman$POS, sep=".")
############################################################################################################################################
## Fig. S2b,c: summary statistics

# S2c: number of SNVs per cell; number of private SNVs per cell

to.plot <- data.frame(Counts = c(colSums(lee.six.mutect.strelka[,-c(1:2)]!=0, na.rm = T), 
                                 colSums(lee.six.mutect.strelka[rowSums(lee.six.mutect.strelka[,-c(1:2)]!=0)>1,-c(1:2)], na.rm = T),
                                 colSums(lee.six.mutect.strelka[rowSums(lee.six.mutect.strelka[,-c(1:2)]!=0)==1,-c(1:2)], na.rm = T),
                                 colSums(lee.six.caveman[, -c(1:4)], na.rm=T),
                                 colSums(lee.six.caveman[rowSums(lee.six.caveman[,-c(1:4)])>1,-c(1:4)], na.rm=T),
                                 colSums(lee.six.caveman[rowSums(lee.six.caveman[,-c(1:4)])==1,-c(1:4)], na.rm=T)),
                      Type=c(rep("All", ncol(lee.six.mutect.strelka)-2), 
                             rep("Shared", ncol(lee.six.mutect.strelka) - 2),
                             rep("Private", ncol(lee.six.mutect.strelka)-2),
                             rep("All", ncol(lee.six.caveman)-4),
                             rep("Shared", ncol(lee.six.caveman) - 4),
                             rep("Private", ncol(lee.six.caveman)-4)),
                      Calling = c(rep("Mutect/Strelka", ncol(lee.six.mutect.strelka)*3-6),
                                  rep("Caveman", ncol(lee.six.caveman)*3-12)))


pdf(paste0(analysis.directory, "Figures/Figure_S2c.pdf"), width=7, height=3)

ggplot(to.plot, aes(x=Calling, y=Counts)) + ggbeeswarm::geom_beeswarm() + facet_wrap(~Type)

dev.off()


## S2b: concordance between the 2 callers

pdf(paste0(analysis.directory, "Figures/Figure_S2b.pdf"), width=3.5, height=3)

ggVennDiagram(list(Mutect_Stelka=paste(lee.six.mutect.strelka$chr, lee.six.mutect.strelka$pos,  sep="."),
            Caveman=paste(lee.six.caveman$X.CHROM, lee.six.caveman$POS, sep=".")), label_alpha=0, edge_lty=0)

dev.off()


############################################################################################################################################
## Fig. 3b: generate a pseudo-bulk from the Mutect/Strelka calls 

lee.six.pb <- rowSums(lee.six.mutect.strelka[,-c(1:2)]!=0, na.rm=T)/(rowSums(!is.na(lee.six.mutect.strelka[,-c(1:2)]!=0)))/2
snvs <- list(Lee_Six = lee.six.pb)
save(snvs, file=paste0(data.directory, "/Lee-Six_et_al/Mutect_Strelka/SNVs.RData"))


pdf(paste0(analysis.directory, "Figures/Figure_3b.pdf"), width=4, height=3)

to.plot <- data.frame(VAF=lee.six.pb)

ggplot(to.plot[to.plot$VAF>=0.01,], aes(x=VAF)) + geom_histogram(binwidth = 0.01) + 
  scale_x_continuous(name="Variant allele frequency", limits = c(0, 1)) +
  scale_y_continuous(name="Number of SSNVs")

dev.off()



############################################################################################################################################
# Fig. S3d: Generate a pseudo-bulk from the original Caveman calls 

lee.six.pb <- rowSums(lee.six.caveman[,-c(1:4)]!=0, na.rm=T)/(rowSums(!is.na(lee.six.caveman[,-c(1:4)]!=0)))/2


pdf(paste0(analysis.directory, "Figures/Figure_S3d.pdf"), width=4, height=3)
to.plot <- data.frame(VAF=lee.six.pb)

ggplot(to.plot[to.plot$VAF>=0.01,], aes(x=VAF)) + geom_histogram(binwidth = 0.01) + 
  scale_x_continuous(name="Variant allele frequency", limits = c(0, 1)) +
  scale_y_continuous(name="Number of SSNVs") 


dev.off()

snvs <- list(Lee_Six = lee.six.pb)
save(snvs, file="./Lee-Six_et_al/Caveman/SNVs.RData")

