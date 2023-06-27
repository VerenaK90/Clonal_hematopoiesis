library("rstudioapi")   
setwd(dirname(getActiveDocumentContext()$path)) 
source("Simulate_trees.R")
source("Post_processing.R")

## set parameters for test run
N = 50
NP = 100
mut.rate = 1
time.max = 10
parms.exp=c(lambda.s=1, delta.s=0, lambda.p=1, delta.p=0, alpha.s=0.5)
parms.steady=c(lambda.s=1, delta.s=1, alpha.s=0, lambda.p=1.023, delta.p=1)

###########################################################################################################################################
#### test initialisation

set.seed(310322)

test <- .initialize.sim.s.p()

#### then let the stem cells divide 10 times

for(i in 1:10){
  test <- .division.s.p(cell.type = 1, tree = test, mutation.mode = "constant", mut.rate = mut.rate)
}
## after 10 divisions, the expected cell number is 11
length(test$tip.class)==11
## after 10 divisions, the expected mutation count is 2*20 + 1(from the original cell)
sum(test$edge.length)==(2*length(test$tip.class)-1)+1


#### now, let the stem cells die 5 times

for(i in 1:5){
  test <- .loss(cell.type = 1,tree = test, nr = 1)
}
## after 5 losses, the expected cell number is 6
length(test$tip.class)==6


###########################################################################################################################################
## this looks all okay. Let's now try to simulate an expanding population of up to 1,000 cells with selection

## set parameters for test run
N = 1000
NP = 0
mut.rate = 1
time.max = 10
parms.exp=c(lambda.s=1, delta.s=0, lambda.p=1, delta.p=0, alpha.s=0.5)
parms.steady=c(lambda.s=1, delta.s=1, alpha.s=0, lambda.p=1, delta.p=1.1)
mut.rate.D = 0.002
mutation.mode="constant"
s.shape = 1.5
s.rate = 35


test = gillespie.sim.s.p(parms.exp = parms.exp, parms.steady = parms.steady, time.max = time.max, time.samples = c(0), N = N,
                         NP = NP, mut.rate = mut.rate, mutation.mode = mutation.mode, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate)

plot(test[[2]], show.tip.label = F)

## increase tree to 10^5 cells
N = 10^3

test = gillespie.sim.s.p(parms.exp = parms.exp, parms.steady = parms.steady, time.max = time.max, time.samples = c(0), N = N,
                         NP = NP, mut.rate = mut.rate, mutation.mode = mutation.mode, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate)

## sample 100 cells from this tree
sample <- .simulate.sampling(test[[2]], 100)
plot(sample, show.tip.label = F)

###########################################################################################################################################
## Finally, test for asymmetric division

## set parameters for test run
N = 1000
NP = 0
mut.rate = 1
time.max = 10
parms.exp=c(lambda.s=1, delta.s=0, lambda.p=1, delta.p=0, alpha.s=0.5)
parms.steady=c(lambda.s=1, delta.s=1, lambda.a=4, alpha.s=0, lambda.p=1, delta.p=1.2)
mut.rate.D = 0.002
mutation.mode="constant"
s.shape = 1.5
s.rate = 35


test = gillespie.sim.s.p(parms.exp = parms.exp, parms.steady = parms.steady, time.max = time.max, time.samples = c(0), N = N,
                         NP = NP, mut.rate = mut.rate, mutation.mode = mutation.mode, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate)

plot(test[[2]], show.tip.label = F)

table(test[[2]]$tip.class)
## subset on the stem cells

stem.cells <- .subset.cell.type(test[[2]], 1)

plot(stem.cells, show.tip.label=F)
plot(.simulate.sampling(stem.cells, 100), show.tip.label=F)
## expected nr of mutations per cell: lambda*mutation.rate*time

muts.per.stem.cell <- get_mutations_per_tip(stem.cells)

3*1*10*1000
sum(stem.cells$edge.length)


progenitor.cells <- .subset.cell.type(test[[2]], 2)

plot(progenitor.cells, show.tip.label=F)

## expected nr of mutations per cell: lambda*mutation.rate*time

muts.per.progenitor.cell <- get_mutations_per_tip(progenitor.cells)
