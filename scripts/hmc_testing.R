rm(list=ls())
library(ContPhylo)

#HMC sampled trees
trees=hmc(n=5,niter=2e5,thin = 1e2,scale=2)
pdf('hmc.pdf',7,7)
prior_test(trees,scale=2)
dev.off()
browseURL('hmc.pdf')
