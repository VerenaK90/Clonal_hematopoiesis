import sys

tumor_id=sys.argv[1]

import math

from pyabc.external.r import R
print(tumor_id+"/Run_model.R")

r = R(tumor_id+"/Run_model.R")

model = r.model("myModel")
distance = r.distance("myDistance")
sum_stat = r.summary_statistics("mySummaryStatistics")

observation = r.observation("mySumStatData")

from pyabc import Distribution, RV, ABCSMC, sampler

prior = Distribution(mu=RV("uniform", 0.1, 9.9), offset=RV("uniform", 0, 10),  N=RV("uniform",2.5,5.5), delta_exp=RV("uniform",0,0.75), lambda_ss=RV("uniform", -3, 2), t_s=RV("uniform", 0, 1), s=RV("uniform", 0, 0.99))

sample_specs=sampler.MulticoreParticleParallelSampler(n_procs=8)

abc = ABCSMC(model, prior, distance, population_size = 1000, sampler=sample_specs)

import os
from tempfile import gettempdir

db = "sqlite:///" + "/home/koerber/pyABC/Blood/" + tumor_id + "/Model_fit.db"

abc.new(db, r.observation("mySumStatData"))

history = abc.run(minimum_epsilon=0.00005, max_nr_populations=25)