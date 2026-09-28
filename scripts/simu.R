rm(list=ls())
library(ContPhylo)
library(ape)
library(phangorn)
set.seed(0)
tree=rtree(10)
scale=0.1
tree$edge.length=tree$edge.length*scale
seqs=phangorn::simSeq(tree,l=1000,rate=1)
njtree=NJ(dist.ml(seqs))
r=hmc(seqs,initTree = tree,niter = 1e4,scale=scale)

par(mfrow=c(2,2))
plot(sapply(r,function(x) x$posterior),type='l',xlab='',ylab='Posterior')

p=posterior(phangorn::pml(tree=tree,data=seqs,model='JC'),r=phylo2contphylo(tree),scale=scale)
plot(midpoint(tree),main=sprintf('Correct tree, lpost=%.2f',p))
axisPhylo(1,backward = F)

p=posterior(phangorn::pml(tree=njtree,data=seqs,model='JC'),r=phylo2contphylo(midpoint(njtree)),scale=scale)
plot(midpoint(njtree),main=sprintf('NJ tree, lpost=%.2f',p))
axisPhylo(1,backward = F)

p=posterior(phangorn::pml(tree=r[[length(r)]],data=seqs,model='JC'),r=phylo2contphylo(r[[length(r)]]),scale=scale)
plot(midpoint(r[[length(r)]]),main=sprintf('Final tree, lpost=%.2f',p))
axisPhylo(1,backward = F)

