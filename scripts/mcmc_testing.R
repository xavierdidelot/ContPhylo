rm(list=ls())
library(ape)
library(ContPhylo)

#Random independent trees
trees=lapply(1:1000,function(x) rtopology(5,rooted=T,br=function (n) rexp(n,rate=1/2)))
prior_test(trees,2)

#My own MCMC
trees=mcmc(phyDatEmpty(5),niter=1e5,scale=2,model=3)
prior_test(trees,scale=2)
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
trees=mcmc2(phyDatEmpty(5),niter=2e5,thin = 1e2,scale=2)
mat=trees[[1]]$cov.jump
image(1:nrow(mat),1:ncol(mat),mat,xlab='',ylab='')
title('Covariance of last jump distribution')
pdf('mcmc.pdf',7,7)
trees2=trees[(length(trees)/2):length(trees)]
prior_test(trees2,scale=2)
dev.off()
browseURL('mcmc.pdf')

#Computing the approximate topological ESS from rwty
library(rwty)
filename=tempfile()
f=file(filename,'w')
for (i in 1:length(trees2)) write.tree(trees2[[i]],file=f,append=T)
close(f)
r=load.trees(filename,type = "newick")
print(topological.approx.ess(r,burnin=0))
