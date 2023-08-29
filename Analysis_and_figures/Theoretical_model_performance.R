###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Illustrate the theoretical model behaviour

source("./Settings.R")

####################################################################################################################
## Figure 1c: how does the VAF distribution shaped by neutral evolution change over time when analyzing different stem cell numbers and a constant division rate of 1/y?

## set lambda and mu constant
lambda = 1/365
mu <- 10

to.plot <- data.frame()

for(N in c(5*10^3, 5*10^4, 5*10^5)){
  for(age in seq(0, 100, 10)){ # simulate for varying agess

    to.plot <- rbind(to.plot,
                     data.frame(VAF=seq(0.05, 1, 0.01), Number_of_Variants = sapply(2*seq(0.05, 1, 0.01), function(x){
                       mutational.burden(mu=mu, N=N, lambda.exp = 1, delta.exp = 0.2, lambda.ss = lambda, t.end = age*365, b = x*N)
                     }), Age = age, Stem_cell_count=N))

  }
}

pdf(paste0(analysis.directory, "/Figures/Figure_1_c.pdf"), width=6.5, height=6.5, useDingbats = F)

ggplot(to.plot, aes(x=1/VAF, y=Number_of_Variants, col=Age, group=Age)) +
  geom_line() +
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1"), limits = c(0, 100)) +
  scale_y_continuous(  name="Cumulative # of SSNVs") +
  theme(aspect.ratio = 1, legend.position = "bottom") + facet_wrap(~Stem_cell_count, scale="free")

dev.off()


####################################################################################################################
## Fig. S1a: How does the model behave with varying N, while leaving N/lambda constant?

age <- 100 # evaluate at fixed age
mu <- 1
to.plot <- data.frame()

for(Ntau in c(10^3, 10^4, 10^5)){
  print(Ntau)
  for(N in c(1000, 10000, 100000, 10^6)){
    lambda <- N/Ntau
    #  mu <- mu_per_year/lambda
    to.plot <- rbind(to.plot,
                     data.frame(VAF=seq(0.05, 1, 0.01), Number_of_Variants = sapply(2*seq(0.05, 1, 0.01), function(x){
                       mutational.burden(mu=mu, N=N, lambda.exp = 1, delta.exp = 0.2, lambda.ss = lambda, t.end = age, b = x*N)
                     }),  Stem_cell_count=N, Division_rate = lambda, Mutation_rate=mu, Ntau=Ntau, Phases = "both"),
                     data.frame(VAF=seq(0.05, 1, 0.01), Number_of_Variants = sapply(2*seq(0.05, 1, 0.01), function(x){
                       mutational.burden(mu=mu, N=N, lambda.exp = 1, delta.exp = 0.2, lambda.ss = lambda, t.end = age, b = x*N, phase="early")
                     }),  Stem_cell_count=N, Division_rate = lambda, Mutation_rate=mu, Ntau=Ntau, Phases = "expansion"),
                     data.frame(VAF=seq(0.05, 1, 0.01), Number_of_Variants = sapply(2*seq(0.05, 1, 0.01), function(x){
                       mutational.burden(mu=mu, N=N, lambda.exp = 1, delta.exp = 0.2, lambda.ss = lambda, t.end = age, b = x*N, phase="homeostasis")
                     }),  Stem_cell_count=N, Division_rate = lambda, Mutation_rate=mu, Ntau=Ntau, Phases = "homeostasis"))

  }
}

pdf(paste0(analysis.directory, "/Figures/Figure_S1_1.pdf"), width=9, height=9, useDingbats = F)

ggplot(to.plot, aes(x=1/VAF, y=Number_of_Variants, col=log10(Stem_cell_count), group=Stem_cell_count)) +  geom_line() +
  scale_color_gradientn(colors=rev(hcl.colors(n = 7, palette="Zissou 1")), limits=c(3, 8))+
  scale_y_continuous(  name="Cumulative # of SSNVs") + facet_wrap(~Ntau*Phases, scales = "free_y") +
  #scale_color_continuous(trans="log10")+
  theme(aspect.ratio = 1, legend.position = "bottom")

# normalize to 1
to.plot.2 <- to.plot
for(phase in unique(to.plot.2$Phases)){
  for(ntau in unique(to.plot.2$Ntau)){
    for(N in unique(to.plot.2[to.plot.2$Ntau==ntau,]$Stem_cell_count)){
      to.plot.2[to.plot.2$Stem_cell_count==N & to.plot.2$Ntau==ntau & to.plot.2$Phases==phase,]$Number_of_Variants <- to.plot.2[to.plot.2$Stem_cell_count==N & to.plot.2$Ntau==ntau & to.plot.2$Phases==phase,]$Number_of_Variants/max(to.plot[to.plot$Stem_cell_count==N & to.plot$Ntau==ntau & to.plot$Phases=="both",]$Number_of_Variants)
    }
  }
}

ggplot(to.plot.2, aes(x=1/VAF, y=Number_of_Variants, col=log10(Stem_cell_count), group=Stem_cell_count)) + geom_line() +
  scale_color_gradientn(colors=rev(hcl.colors(n = 7, palette="Zissou 1")), limits=c(3, 8))+
  scale_y_continuous(  name="Cumulative # of SSNVs") + facet_wrap(~Ntau*Phases, scales = "free") +
  theme(aspect.ratio = 1, legend.position = "bottom")

dev.off()


##############################################################################################################################################
## Figure 1g: predicted selection over time - 25,000 stem cells that divide 10 times per year. A selected clone is introduced at 20 years and grows with a selective advantage of 0.98

vafs.of.interest <- seq(0.01, 1, 0.005)

## parameters
lambda.exp <- 1
delta.exp <- 0
lambda.ss <- 10/365

N <- 25000
mu <- 1

## selective advantage (as a reduction in cell death)
s <- 0.98
## age at which selective advantage is acquired
t.s <- 20*365

## simulate the burden with these parameters under selection
simulated.burden.selection <- c()
for(age in c(40,50, 52, 54, 56, 58, 60, 70)){
  simulated.burden.selection <- rbind(simulated.burden.selection,
                                      data.frame(VAF= vafs.of.interest,
                                                 SSNVs = sapply(2*vafs.of.interest, function(b){
                                                   mutational.burden.with.selection(mu, N, lambda.exp, delta.exp, lambda.ss, t.end=age*365, t.s = t.s, s=s, b= b*N, min.clone.size = 0.01)
                                                 }), Time=age, Clone.size =  min(N, exp(lambda.ss*(1-s)*(age*365-t.s)))))

}

pdf(paste0(analysis.directory, "/Figures/Figure_1_g.pdf"), width=3.5, height=3.5, useDingbats = F)

ggplot(simulated.burden.selection[simulated.burden.selection$VAF>=0.05,], aes(x=1/VAF, y=SSNVs, col=Time, group=Time)) +
  geom_line() +
  scale_color_gradientn(colors=rev(hcl.colors(n = 7, palette="RdYlBu")))+
  scale_y_continuous(name="Number of SSNVs")  + theme(aspect.ratio = 1)

dev.off()

##############################################################################################################################################
## Figure 1f: Predicted selection for varying t.s but at the same clone size

vafs.of.interest <- seq(0.01, 1, 0.005)

## parameters
lambda.exp <- 1
delta.exp <- 0
lambda.ss <- 10/365

N <- 25000
mu <- 1

## selective advantage (as a reduction in cell death)
s <- 0.98

## simulate clones of equal size but with start and end at different time points; in total 45 years expansion, amounting to a clone size of 32%
growth.time <- 45

simulated.burden.selection <- c()
for(age in c(50, 60, 70, 80, 90)){
  simulated.burden.selection <- rbind(simulated.burden.selection,
                                      data.frame(VAF= vafs.of.interest,
                                                 SSNVs = sapply(2*vafs.of.interest, function(b){
                                                   mutational.burden.with.selection(mu, N, lambda.exp, delta.exp, lambda.ss, t.end=age*365,
                                                                                    t.s = (age-growth.time)*365, s=s, b= b*N, min.clone.size = 0.01)
                                                 }), Time=(age-growth.time)*365, Clone.size =  min(N, exp(lambda.ss*(1-s)*(growth.time*365)))))

}


pdf(paste0(analysis.directory, "/Figures/Figure_1_f.pdf"), width=3.5, height=3.5, useDingbats = F)

ggplot(simulated.burden.selection[simulated.burden.selection$VAF>=0.05,], aes(x=1/VAF, y=SSNVs, col=Time/365, group=Time)) +
  geom_line() +
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1"))+
  scale_y_continuous(name="Number of SSNVs")  + theme(aspect.ratio = 1)

dev.off()


##############################################################################################################################################
## Sigure 1h: Predicted selection for varying s but at the same age of 70 years

vafs.of.interest <- seq(0.01, 1, 0.005)

## parameters
lambda.exp <- 1
delta.exp <- 0
lambda.ss <- 10/365

N <- 25000
mu <- 1

t.s <- 50
age <- 70

## simulate the burden with these parameters under selection
simulated.burden.selection <- c()
for(s in seq(0.95, 0.99, 0.002)){
  simulated.burden.selection <- rbind(simulated.burden.selection,
                                      data.frame(VAF= vafs.of.interest,
                                                 SSNVs = sapply(2*vafs.of.interest, function(b){
                                                   mutational.burden.with.selection(mu, N, lambda.exp, delta.exp, lambda.ss, t.end=age*365, t.s = t.s*365, s=s, b= b*N, min.clone.size = 0.01)
                                                 }), s=s, Clone.size =  min(N, exp(lambda.ss*(1-s)*(age*365-t.s*365)))))

}
simulated.burden.selection$growth_per_year <- (exp(lambda.ss*(1-simulated.burden.selection$s)*365)-1)*100


pdf(paste0(analysis.directory, "/Figures/Figure_1_h.pdf"), width=3.5, height=3.5, useDingbats = F)

ggplot(simulated.burden.selection[simulated.burden.selection$VAF>=0.05,], aes(x=1/VAF, y=SSNVs, col=growth_per_year, group=1-s)) +
  geom_line() + scale_x_continuous(breaks=c(5, 10, 20), labels = c("0.2", "0.1", "0.05"),
                                   name="Variant allele frequency distribution")+
  scale_color_gradientn(colors=rev(hcl.colors(n = 7, palette="Zissou1")))+
  scale_y_continuous(name="Number of SSNVs")  + theme(aspect.ratio = 1)

dev.off()


##############################################################################################################################################
## Figure 1i: Predicted selection for varying t.s and evaluating at the same age

vafs.of.interest <- seq(0.01, 1, 0.005)

## parameters
lambda.exp <- 1
delta.exp <- 0
lambda.ss <- 10/365

N <- 25000
mu <- 1

## selective advantage (as a reduction in cell death)
s <- 0.98

## simulate clones of equal size but with start and end at different time points; in total 45 years expansion, amounting to a clone size of 32%
t.s <- seq(20, 30)
age <- 70

simulated.burden.selection <- c()
for(t.s in t.s){
  simulated.burden.selection <- rbind(simulated.burden.selection,
                                      data.frame(VAF= vafs.of.interest,
                                                 SSNVs = sapply(2*vafs.of.interest, function(b){
                                                   mutational.burden.with.selection(mu, N, lambda.exp, delta.exp, lambda.ss, t.end=age*365,
                                                                                    t.s = t.s*365, s=s, b= b*N, min.clone.size = 0.01)
                                                 }), Time=t.s*365, Clone.size =  min(N, exp(lambda.ss*(1-s)*(age*365-t.s*365)))))

}


pdf(paste0(analysis.directory, "/Figures/Figure_1_i.pdf"), width=3.5, height=3.5, useDingbats = F)

ggplot(simulated.burden.selection[simulated.burden.selection$VAF>=0.05,], aes(x=1/VAF, y=SSNVs, col=Time/365, group=Time)) +
  geom_line() +
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1"))+
  scale_y_continuous(name="Number of SSNVs")  + theme(aspect.ratio = 1)

dev.off()


