rm(list=ls())
library(ape)
library(ContPhylo)

n=3#number of leaves in tree
model=1#model of tree
scale=1#scaling factor on branch lengths
simu=2 #set to 1 for independent trees, 2 for my mcmc, 3 for adaptMCMC

if (simu==1) {
  #Random independent trees
  if (model==1) trees=lapply(1:1000,function(x) {y=rcoal(n);y$edge.length=y$edge.length*scale;return(y)})
  if (model==2) trees=lapply(1:1000,function(x) rtopology(n,rooted=T,br=function (n) rexp(n,rate=1/scale)))
  if (model==3) trees=lapply(1:1000,function(x) rtopology(n,rooted=F,br=function (n) rexp(n,rate=1/scale)))
}

if (simu==2) {
  #My own MCMC
  trees=mcmc(phyDatEmpty(n),model=model,scale=scale,niter=1e5)
  par(mfrow=c(2,2))
  for (i in c(1700,1800,1900,2000)) {
    plot(trees[[i]])
    axisPhylo(1,backward = F)
  }
  par(mfrow=c(2,2))
  plot(sapply(trees,function(x) x$paccepted[1]),type='l',xlab='',ylab='paccepted[1]')
  plot(sapply(trees,function(x) x$paccepted[2]),type='l',xlab='',ylab='paccepted[2]')
  plot(sapply(trees,function(x) x$sd[1]),type='l',xlab='',ylab='sd[1]')
  plot(sapply(trees,function(x) x$sd[2]),type='l',xlab='',ylab='sd[2]')
  trees=trees[(length(trees)/2):length(trees)]
}

if (simu==3) {
  #MCMC from adaptMCMC package
  trees=mcmc2(phyDatEmpty(n),model=model,scale=scale,niter=1e5)
  mat=trees[[1]]$cov.jump
  image(1:nrow(mat),1:ncol(mat),mat,xlab='',ylab='')
  title('Covariance of last jump distribution')
  trees=trees[(length(trees)/2):length(trees)]
}

pdf('mcmc.pdf',7,7)
prior_test(trees,model=model,scale=scale)
dev.off()
browseURL('mcmc.pdf')

if (F) {
  #Computing the approximate topological ESS from rwty
  library(rwty)
  filename=tempfile()
  f=file(filename,'w')
  for (i in 1:length(trees2)) write.tree(trees2[[i]],file=f,append=T)
  close(f)
  r=load.trees(filename,type = "newick")
  print(topological.approx.ess(r,burnin=0))
}
