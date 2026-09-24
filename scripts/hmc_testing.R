rm(list=ls())
library(ContPhylo)

#HMC sampled trees
trees=hmc(n=5,niter=1e4,scale=2)
prior_test(trees,scale=2)
