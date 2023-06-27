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

## compute the cumulative number of variants per VAF of interest
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
  x[x==0] <- min(x[x!=0]) # if there's no error, set to minimal error to avoid problems in cost function
  x
})

mySumStatData <- list(mutation.count=mutation.count, sampled.sd=sampled.sd, age = c(rep(age, length(vafs.of.interest))))


###### create the Bayesian setup

myModel <- function(parms){

  N <- round(10^parms$N) # prior is on log scale
  delta.exp <- parms$delta_exp
  lambda.ss <- 10^parms$lambda_ss # prior is on log scale
  mu <- parms$mu
  offset <- parms$offset

  min.t.s <- 0
  max.t.s <- age - log(max(1, min.prior.size*N))/lambda.ss # latest time point to reach the minimal prior size with unrestricted growth (delta = 0)
  
  t.s <- min.t.s + parms$t_s*(max.t.s - min.t.s)
  
  min.s <- (lambda.ss*(age-max(t.s)) - log(N))/(lambda.ss*(age-max(t.s))) # at this selective advantage, the clone will get fixed by the sampling time point. Smaller values of s hence make no sense for the observed subclones
  if(min.s < 0){
    min.s <- 0
  }
  max.s <- (lambda.ss*(age-min(t.s)) - log(max(1, min.prior.size*N)))/(lambda.ss*(age-min(t.s))) # at this selective advantage, the clone will reach at least the minimal size at the time point of observation. Larger values of s will not result in a visible clone. 
  
  s <- min.s + parms$s*(max.s-min.s)

  if(s < 0 & age >0){ # s must be at least 0
    modelResult <- rep(10^5, length(vafs))
    return( list(modelResult=modelResult))
  }else if(s<0 & age==0){
   ## if age = 0: no time for selection yet
   s <- 0.99
   t.s <- 0	
  }

   
  modelResult <- list()
  
  # define the clone sizes of interest: start from min.vaf/2 (where min.vaf is the smallest observable vaf) to account for stochastic fluxes from this interval into the measurable region
  clone.sizes.for.fitting <- c(min.vaf/2, seq(min.vaf, 1, 0.01))   

  # don't take clones that would round to zero cells
  clone.sizes.for.fitting <- clone.sizes.for.fitting[round(clone.sizes.for.fitting*N)>0]

  # simulate the mutation counts
  sim <-  mutational.burden.with.selection(mu=mu, N=N, lambda.exp=1, delta.exp=delta.exp, min.clone.size = min.clone.size,
                                       lambda.ss=lambda.ss, t.end=age, t.s=t.s[1], s=s, round(clone.sizes.for.fitting*N))
  
  # convert the cumulative counts into counts per interval
  sim <- sim - c(sim[-1], 0)
    
  ## add offset to clonal bin 
  sim[clone.sizes.for.fitting==1] <- sim[clone.sizes.for.fitting==1] + offset
  
  # simulate sequencing  
  sim.vafs <- simulated.data(seq.type, clone.sizes.for.fitting, sim, depth=depth, ncells=ncells, sensitivity= use.sensitivity, false.negative.per.vaf = false.negative.per.vaf)
  
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
