#################################################################################################################################################
############ Simulate selected trees 

library(FLORENCE)

library(doParallel)
library(foreach)
library(parallel)

### Simulate trees with a stem cell count of 25 000 with a driver introduced at 20 years and read out at 5, 10, ..., % clone size
###  stem cells divide ten times a year and acquire 1 mutation per division. We simulate stem cells only, thus set progenitor parameters to zero.

N = 25000
NP = 0
mut.rate = 1
time.max = 75
parms.exp=c(lambda.s=1, delta.s=0, lambda.p=0, delta.p=0, alpha.s=0)
parms.steady=c(lambda.s=10, delta.s=10, alpha.s=0, lambda.p=0, delta.p=0)
driver.mode="fixed_time"
t.driver = 20
mut.rate.D = 0
mutation.mode="constant"
s.shape = 1.4
s.rate = 70
time.samples = seq(0, 75, 25)
report.at.f <- c(0.01, 0.02, 0.04, 0.05, 0.075, seq(0.1, 1, 0.05))

## parallelize the runs
n.cores <- 7
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

trees <- foreach(sim.nr = 1:10) %dopar% {
  set.seed(sim.nr*354)
  source("../Functions/Simulate_trees.R")
  print(sim.nr)
  gillespie.sim.s.p(parms.exp, parms.steady, time.max=75, time.samples=time.samples, N=N, NP=NP,
                    mut.rate = mut.rate, mutation.mode="Binomial", driver.mode = driver.mode, t.driver = t.driver,
                    mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate,
                    tau = N/10, report.at.f = report.at.f)
}
save(trees, file=paste0("RData/Simulated_data/Selected_tree_25000_HSCs_t_s_20_sel_0.02.RData"))


#################################################################################################################################################
############ Compute the VAFs


for(sim.nr in 1:length(trees)){
  print("...")
  print(sim.nr)
  ## downsampling to 10000 cells
  sampled.trees <- lapply(trees[[sim.nr]], function(x){
    .simulate.sampling(x, min(10000, length(x$tip.class)))
  })
  
  ## simulate VAF
  vafs <- lapply(sampled.trees, function(x){
    vafs <- get_vaf_from_tree(x)
    vafs
  })
  
  save(vafs, file=paste0("./RData/Simulated_data/True_VAFs_selection_", N, "_N_simnr_", sim.nr, ".RData")) ## no sequencing simulated on top
}

