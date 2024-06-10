###### load observed data

## specify the VAFs at which model and data are compared; min.vaf must be given in the Run_model.script or defaults to 0.05

if(!exists("min.vaf")){
	min.vaf <- 0.05
}
## should the sensitivity model be used?
if(!exists("use.sensitivity")){
	use.sensitivity <- T
}
## what lower limit for the clone size should be used to learn VAFs?
if(!exists("min.clone.size")){
	min.clone.size=0.05
}
## what lower limit for the prior clone size should be used?
if(!exists("min.prior.size")){
	min.prior.size=0.01
}
if(!exists("depth")){
  depth=90
}

if(!exists("seq.type")){
  seq.type="bulk" ## default the sequencing analysis to bulkWGS
}

if(!exists("ncells")){
  ncells=100 ## default the number of analyzed cells to 100 (only used in sc-seq analysis)
}


vafs.of.interest <- seq(min.vaf, 1, 0.01)


cumulated.vaf.count.mean <- list()
cumulated.vaf.count.sd <- list()

if(depth <= 100){
  vafs <- lapply(snvs, function(x){
    x[x$Depth>=10 & x$Depth<=300  & as.numeric(x$varCounts)>=3, ]$VAF
  })
}else{
 vafs <- lapply(snvs, function(x){
    x$VAF
  })
}



cumulated.vaf.count.mean <- lapply(vafs, function(y){
    sapply(vafs.of.interest, function(x){
    sum(y>=x)
  })})

## Counting error
cumulated.vaf.count.sd <- lapply(vafs, function(y){sapply(vafs.of.interest, function(x){
  sqrt(sum(y>=x)) 
  })
})


vafs = vafs.of.interest
mutation.count = cumulated.vaf.count.mean
sampled.sd <- lapply(cumulated.vaf.count.sd, function(x){
  x[x==0] <- min(x[x!=0])
  x
})

mySumStatData <- list(mutation.count=mutation.count, sampled.sd=sampled.sd, age = c(rep(age, length(vafs.of.interest))))



###### create the Bayesian setup


myModel <- function(parms){

  ## transform the parameters
  N <- round(10^parms$N)
  delta.exp <- parms$delta_exp
  lambda.ss <- 10^parms$lambda_ss
  mu <- parms$mu
  offset <- parms$offset
  size <- 10^unlist(parms[grep("size", names(parms))])
  ts <- c(0, unlist(parms[grep("ts", names(parms))]))
  s <- c(1, 1 - log(size*N)/(lambda.ss*(age - ts[-1]))) # compute the selective advantages from the size parameters
  s[s<0] <- 0
  s[s>1] <- 1

   
  modelResult <- list()
  
  # define the clone frequencies of interest: start from min.vaf/2 (where min.vaf is the smallest observable vaf) to account for stochastic fluxes from this interval into the measurable region
  clone.frequencies.for.fitting <- c(min.vaf/2, seq(min.vaf, 1, 0.01))   
  # don't take clone frequencies that would round to zero cells
  clone.frequencies.for.fitting <- clone.frequencies.for.fitting[round(clone.frequencies.for.fitting*N)>0]
  # simulate the mutation counts
  sim <- mutational.burden.multiclone(mu = mu, N = N, lambda.exp = 1, delta.exp = delta.exp, min.clone.size = min.clone.size, lambda.ss = lambda.ss, t.end=age, t.s = ts, s = s, mother.daughter = mother.daughter, b = round(clone.frequencies.for.fitting*N))

  # convert the cumulative counts into counts per interval
  sim <- sim - c(sim[-1], 0)
  
  ## add offset to clonal bin 
  sim[clone.frequencies.for.fitting==1] <- sim[clone.frequencies.for.fitting==1] + offset
  
  # simulate sequencing
  sim.vafs <- simulated.data(seq.type, clone.frequencies.for.fitting, sim, depth=depth, ncells=ncells, sensitivity= use.sensitivity, false.negative.per.vaf = false.negative.per.vaf)
  
  # compute cumulative counts
  sim.vafs <- sapply(vafs, function(x){
    sum(sim.vafs >=x)
  })
  
  modelResult[[1]] <- sim.vafs


  list(modelResult=modelResult)
 
}

mySummaryStatistics <- function(modelResult){
  modelResult
}

# compute distance between model and data
myDistance <- function(modelResult, mySumStatData){
     res <- sum((mySumStatData$mutation.count[[1]] - modelResult$modelResult[[1]])^2)

  res
}
