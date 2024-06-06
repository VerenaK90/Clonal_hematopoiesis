###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Illustrate the theoretical model behaviour

source("./Settings.R")

####################################################################################################################
## Figure 1e-g / Extended Data Fig. 1a: how does the VAF distribution change over time in the absence of selection, for different stem cell numbers and a constant division rate of 1/y?

vafs.of.interest <- 10^seq(log10(0.01), 0, length.out=100) 

## set lambda and mu constant
lambda = 5/365
mu <- 10

to.plot <- data.frame()

for(N in c(5*10^3, 10^4, 2*10^4, 5*10^4, 10^5, 2*10^5, 5*10^5, 10^6)){
  for(age in seq(0, 100, 10)){ # simulate for varying agess

    to.plot <- rbind(to.plot,
                     data.frame(VAF=vafs.of.interest, Number_of_Variants = sapply(2*vafs.of.interest, function(x){
                       mutational.burden(mu=mu, N=N, lambda.exp = 1, delta.exp = 0.2, lambda.ss = lambda, t.end = age*365, b = x*N)
                     }), Age = age, Stem_cell_count=N))

  }
}

pdf(paste0(analysis.directory, "/Figures/Figure_1efg_Extended_Data_Figure_1a.pdf"), width=6.5, height=6.5, useDingbats = F)

ggplot(to.plot, aes(x=1/VAF, y=Number_of_Variants, col=Age, group=Age)) +
  geom_line() +
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1"), limits = c(0, 100)) +
  scale_y_continuous(  name="Cumulative # of SSNVs") +
  theme(aspect.ratio = 1, legend.position = "bottom") + facet_wrap(~Stem_cell_count, scale="free")

dev.off()


####################################################################################################################
## Supplementary fig. 1: How does the SFS of clones generated during homeostasis behave in the limit of long times?

# exact and approximate analytical solutions for the SFS of variants generated during homeostasis:
ana <- function(mu, lambda, t, f, N){
  mu/f*(lambda*t/(1+lambda*t))^(f*N)
}
ana.approx <- function(mu, lambda, t, f, N){
  mu/f*exp(-f*N/(lambda*t))
}

# exact and approximate solutions
to.plot <- data.frame()
# exact cumulative distribution
to.plot.cum <- data.frame()

# model parameters
N <- 10^3
mu <- 1
lambda <- 1

for(t in c(2, 10, 100, 1000)){
  to.plot <- rbind(to.plot, data.frame(time=t, S=ana(mu, lambda, t, vafs.of.interest, N),
                                       f=vafs.of.interest, Method="Exact"))
  to.plot <- rbind(to.plot, data.frame(time=t, S=ana.approx(mu, lambda, t, vafs.of.interest, N),
                                       f=vafs.of.interest, Method="Approx"))
  to.plot.cum <- rbind(to.plot.cum, data.frame(time=t, M=mutational.burden(mu, N, 1, 0, lambda, t, vafs.of.interest*N,
                                                                           phase="homeostasis"),
                                               VAF = vafs.of.interest))
}

# same time, same N/lambda, varying N
t <- 100
Nlambda <- 10^4
to.plot.3 <- data.frame()
for(N in c(10^3, 10^4, 10^5, 10^6)){
  lambda <- N/Nlambda
  to.plot.3 <- rbind(to.plot.3, data.frame(VAF=vafs.of.interest,
                                           M=mutational.burden(mu, N, 1, 0, lambda, t, vafs.of.interest*N),
                                           N=N))
}


pdf(paste0(analysis.directory, "/Figures/Supp_1.pdf"), width=5, height=3.5)

p1 <- ggplot(to.plot, aes(x=f, y=S, col=time, group=interaction(time,Method), linetype=Method)) + geom_line() + 
  scale_x_log10(name="f") + scale_y_log10(limits=c(0.01, 500), name="S(f)") + 
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1"), trans="log") +
  geom_line(data = data.frame(f=vafs.of.interest, S=mu/vafs.of.interest), aes(x=f, y=S), inherit.aes=F,
            col="black")

p2 <- ggplot(to.plot.cum, aes(x=VAF, y=M, col=time, group=time)) + geom_line() + 
  scale_x_log10(name="f") + scale_y_log10(limits=c(1, 10000), name="Cumulative number of variants") + 
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1"), trans="log") +
  geom_line(data = data.frame(f=vafs.of.interest, S=mu*N*log(1/vafs.of.interest)), aes(x=f, y=S), inherit.aes=F,
            col="black")

p3 <- ggplot(to.plot.3[to.plot.3$VAF==0.01,], aes(x=N, y=M, col=N, group=N)) + geom_point() + 
  scale_x_log10(name="N") + scale_y_log10( name="Cumulative number of variants")+ 
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1"), trans="log") 

print(p1)
print(p2)
print(p3)

dev.off()


##############################################################################################################################################
## Figure 1k: Predicted selection for varying t.s but at the same clone size


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


pdf(paste0(analysis.directory, "/Figures/Figure_1k.pdf"), width=3.5, height=3.5, useDingbats = F)

ggplot(simulated.burden.selection[simulated.burden.selection$VAF>=0.05,], aes(x=1/VAF, y=SSNVs, col=Time/365, group=Time)) +
  geom_line() +
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1"))+
  scale_y_continuous(name="Number of SSNVs")  + theme(aspect.ratio = 1)

dev.off()


##############################################################################################################################################
## Figure 1l: predicted selection over time - 25,000 stem cells that divide 10 times per year. A selected clone is introduced at 20 years and grows with a selective advantage of 0.98

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

pdf(paste0(analysis.directory, "/Figures/Figure_1l.pdf"), width=3.5, height=3.5, useDingbats = F)

ggplot(simulated.burden.selection[simulated.burden.selection$VAF>=0.05,], aes(x=1/VAF, y=SSNVs, col=Time, group=Time)) +
  geom_line() +
  scale_color_gradientn(colors=rev(hcl.colors(n = 7, palette="RdYlBu")))+
  scale_y_continuous(name="Number of SSNVs")  + theme(aspect.ratio = 1)

dev.off()

##############################################################################################################################################
## Extended Data Fig. 1b: predicted selection over time - 25,000 stem cells that divide 10 times per year. A selected clone is introduced at 20 years and grows with a selective advantage of 0.98 until it reaches a clone size of VAF 5%, 10%, ..., 30%

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

simulated.burden.selection <- c()
for(N in c(5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)){
  print(N)
  ages = t.s + log(c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)*2*N)/((1-s)*lambda.ss)
  for(age in ages){
    simulated.burden.selection <- rbind(simulated.burden.selection,
                                        data.frame(VAF= vafs.of.interest,
                                                   SSNVs = sapply(2*vafs.of.interest, function(b){
                                                     mutational.burden.with.selection(mu, N, lambda.exp, delta.exp, lambda.ss, t.end=age, t.s = t.s, s=s, b= b*N, min.clone.size = 0.01)
                                                   }), Time=age/365, N = N, Clone.size =  min(N, exp(lambda.ss*(1-s)*(age-t.s)))))
    
  }
  
}

pdf(paste0(analysis.directory, "/Figures/Extended_Data_Figure_1b.pdf"), width=6.5, height=6.5, useDingbats = F)

ggplot(simulated.burden.selection[simulated.burden.selection$VAF>=0.05,], aes(x=1/VAF, y=SSNVs, col=Time, group=Time)) +
  geom_line() + facet_wrap(~N, scales = "free") +
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1")) +
  scale_y_continuous(name="Number of SSNVs")  + theme(aspect.ratio = 1)

dev.off()

##############################################################################################################################################
## Extended Data Fig. 1c: predicted selection over time - 25,000 stem cells that divide 10 times per year. A selected clone is introduced at 20 years and grows with a varying selective advantage until it reaches 10% VAF

## parameters
lambda.exp <- 1
delta.exp <- 0
lambda.ss <- 10/365

mu <- 1

## age at which selective advantage is acquired
t.s <- 20*365

simulated.burden.selection <- c()
for(N in c(5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)){
  for(s in seq(0.05, 0.95, 0.1)){
    print(N)
    age = t.s + log(0.1*2*N)/((1-s)*lambda.ss)
    simulated.burden.selection <- rbind(simulated.burden.selection,
                                        data.frame(VAF= vafs.of.interest,
                                                   SSNVs = sapply(2*vafs.of.interest, function(b){
                                                     mutational.burden.with.selection(mu, N, lambda.exp, delta.exp, lambda.ss, t.end=age, t.s = t.s, s=s, b= b*N, min.clone.size = 0.01)
                                                   }), Time=age/365, N = N, s = s))
    
  }
}

pdf(paste0(analysis.directory, "/Figures/Extended_Data_Figure_1c.pdf"), width=6.5, height=6.5)

ggplot(simulated.burden.selection[simulated.burden.selection$VAF>=0.05,], aes(x=1/VAF, y=SSNVs, col=1-s, group=Time)) +
  geom_line() + facet_wrap(~N, scales = "free") +
  scale_color_gradientn(colors=hcl.colors(n = 7, palette="Zissou1")) +
  scale_y_continuous(name="Number of SSNVs")  + theme(aspect.ratio = 1)

dev.off()
