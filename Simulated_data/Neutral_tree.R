#################################################################################################################################################
############ Simulate neutral trees 

source("./Simulated_data/Simulate_trees.R")
source("./Simulated_data/Post_processing.R")

library(doParallel)
library(foreach)
library(parallel)


### simulate trees with a stem cell count of 25 000 that divide ten times a year and acquire 1 mutation per division.
### simulate stem cells only, thus set progenitor parameters to zero.
### The simulated trees evolve neutrally, hence, set the driver mutation rate to zero

N = 25000
NP = 0
mut.rate = 1
time.max = 75
parms.exp=c(lambda.s=1, delta.s=0, lambda.p=0, delta.p=0, alpha.s=0)
parms.steady=c(lambda.s=10, delta.s=10, alpha.s=0, lambda.p=0, delta.p=0)
driver.mode="fixed_time"
t.driver = 100
mut.rate.D = 0
mutation.mode="constant" # model neutral mutation acquisition at constant rate
time.samples = seq(0, 75, 25)

## parallelize the runs
n.cores <- 7
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

trees <- foreach(sim.nr = 1:10) %dopar% {
  set.seed(sim.nr*354)
  print(sim.nr)
  source("./Simulated_data/Simulate_trees.R")
  
  gillespie.sim.s.p(parms.exp, parms.steady, time.max=75, time.samples=time.samples, N=N, NP=NP,
                    mut.rate = mut.rate, mutation.mode="Binomial", driver.mode = driver.mode, t.driver = t.driver, 
                    mut.rate.D = mut.rate.D, tau = N/100)
}
save(trees, file="RData/Neutral_tree_25000_HSCs.RData")

#################################################################################################################################################
############ Compute the VAFs by downsampling the trees to 10,000 cells

for(sim.nr in 1:length(trees)){
  print(sim.nr)
  print("...")
  pdf(paste0("Neutral_25000_HSCs_Simulation_nr_", sim.nr, ".pdf"))
  par(mfrow=c(2,2))
  ## downsampling to 10000 cells
  sampled.trees <- lapply(trees[[sim.nr]], function(x){
    .simulate.sampling(x, min(10000, length(x$tip.class)))
  })
  
  ## compute VAF
  vafs <- lapply(sampled.trees, function(x){
    vafs <- get_vaf_from_tree(x)
    vafs
  })
  
  save(vafs, file=paste0("./True_VAFs_neutral_25000_N_simnr_", sim.nr, ".RData")) ## pure VAFs, no sequencing simulated on top
  
}
  
 