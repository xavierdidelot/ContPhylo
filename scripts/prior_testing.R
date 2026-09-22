rm(list=ls())
library(ape)
library(ContPhylo)

#Random independent trees
trees=lapply(1:1000,function(x) rtopology(4,rooted=T,br=rexp))
prior_test(trees)

#My own MCMC
trees=mcmc(n=4,niter=1e5)
prior_test(trees)
par(mfrow=c(2,2))
for (i in c(700,800,900,1000)) {
  plot(trees[[i]])
  axisPhylo(1,backward = F)
}
par(mfrow=c(2,2))
plot(sapply(trees,function(x) x$paccepted[1]),type='l',xlab='',ylab='paccepted[1]')
plot(sapply(trees,function(x) x$paccepted[2]),type='l',xlab='',ylab='paccepted[2]')
plot(sapply(trees,function(x) x$sd[1]),type='l',xlab='',ylab='sd[1]')
plot(sapply(trees,function(x) x$sd[2]),type='l',xlab='',ylab='sd[2]')

#MCMC from adaptMCMC package
trees=mcmc2(n=4,niter=1e5)
prior_test(trees)
