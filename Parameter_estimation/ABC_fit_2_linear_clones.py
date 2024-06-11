import sys

tumor_id=T2 # specify ID

age=77*365 # specify age (in days)

import math

import numpy as np

from pyabc.external.r import R

r = R(tumor_id+"Run_model_BM_WGS_data_2_linear_clones.R")

model = r.model("myModel")
distance = r.distance("myDistance")
sum_stat = r.summary_statistics("mySummaryStatistics")

observation = r.observation("mySumStatData")

from pyabc import Distribution, RV, ABCSMC, sampler, random_variables, parameters

class ConstrainedPrior(random_variables.DistributionBase):
        def __init__(self):
                self.N = RV("uniform", 5, 3)
                self.lambda_ss = RV("uniform", -2, 1)
                self.offset = RV("uniform", 0, 10)
                self.ts1 = RV("uniform", 0, age)
                self.ts2 = RV("uniform", 0, age)
                self.delta_exp = RV("uniform", 0, 0.75)
                self.size1 = RV("uniform", -0.6, 0.6)# specify based on results obtained with one-clone model
                self.size2 = RV("uniform", -2, 2)
                self.mu = RV("uniform", 1, 4)

        def rvs(self, *args, **kwargs):
                while True:
                  N, lambda_ss, offset, ts1, ts2, delta_exp, size1, size2, mu = self.N.rvs(), self.lambda_ss.rvs(), self.offset.rvs(), self.ts1.rvs(), self.ts2.rvs(), self.delta_exp.rvs(), self.size1.rvs(), self.size2.rvs(), self.mu.rvs()
                  if (ts2 > ts1) & (size1 > size2) & (np.log(10**size1*10**N)/(10**lambda_ss*(age - ts1)) <= 1) & (np.log(10**size2*10**N)/(10**lambda_ss*(age - ts2)) <= 1) & (np.log(10**size1*10**N)/(10**lambda_ss*(age-ts1)) >= 0) & (np.log(10**size2*10**N)/(10**lambda_ss*(age-ts2)) >= 0) :
                    return parameters.Parameter(ts1=ts1,ts2=ts2, size1=size1, size2=size2, delta_exp=delta_exp, mu=mu, offset=offset, N = N, lambda_ss=lambda_ss)

        def pdf(self, x):
          ts1, ts2, size1, size2, delta_exp, mu, offset, N, lambda_ss = x["ts1"], x["ts2"], x["size1"], x["size2"], x["delta_exp"], x["mu"], x["offset"], x["N"], x["lambda_ss"]
          if not((ts2 > ts1) & (size1 > size2) & (np.log(10**size1*10**N)/(10**lambda_ss*(age - ts1)) <= 1) & (np.log(10**size2*10**N)/(10**lambda_ss*(age - ts2)) <= 1)& (np.log(10**size1*10**N)/(10**lambda_ss*(age-ts1)) >= 0) & (np.log(10**size2*10**N)/(10**lambda_ss*(age-ts2)) >= 0)):
            return 0.0
          return self.ts1.pdf(ts1) * self.ts2.pdf(ts2) * self.size1.pdf(size1) * self.size2.pdf(size2) * self.delta_exp.pdf(delta_exp) * self.mu.pdf(mu) * self.offset.pdf(offset) * self.N.pdf(N) * self.lambda_ss.pdf(lambda_ss)

constrained_prior = ConstrainedPrior()

sample_specs=sampler.MulticoreParticleParallelSampler(n_procs=8)

abc = ABCSMC(model, constrained_prior, distance, population_size = 1000, sampler=sample_specs)

import os
from tempfile import gettempdir

db = "sqlite:///" + "/Model_fits_2_clones/WGS/" + tumor_id + "/Model_fit_linear.db"

abc.new(db, r.observation("mySumStatData"))

history = abc.run(minimum_epsilon=1, max_nr_populations=25)
