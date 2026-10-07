rm(list=ls())
library(ContPhylo)

#HMC sampled trees
trees=hmc(phyDatEmpty(5),niter=2e5,thin = 1e2,scale=2)
pdf('hmc.pdf',7,7)
prior_test(trees,model=2,scale=2)
dev.off()
browseURL('hmc.pdf')
