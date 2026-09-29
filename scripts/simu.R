rm(list=ls())
library(ContPhylo)
library(ape)
library(phangorn)
set.seed(0)
tree=rtree(4)
scale=0.01
tree$edge.length=tree$edge.length*scale
seqs=phangorn::simSeq(tree,l=10000,rate=1)
#seqs=phyDatEmpty(Ntip(tree));names(seqs)<-tree$tip.label
r=hmc(seqs,initTree = tree, niter = 1e3,scale=scale)

par(mfrow=c(2,2))
plot(sapply(r,function(x) x$posterior),type='l',xlab='',ylab='Posterior')

p=posterior(phangorn::pml(tree=tree,data=seqs,model='JC'),r=phylo2contphylo(tree),scale=scale)
plot(midpoint(tree),main=sprintf('Correct tree, lpost=%.2f',p))
axisPhylo(1,backward = F)

if (length(seqs)>2 && all(!is.na(dist.ml(seqs)))) {
  njtree=NJ(dist.ml(seqs))
  p=posterior(phangorn::pml(tree=njtree,data=seqs,model='JC'),r=phylo2contphylo(midpoint(njtree)),scale=scale)
  plot(midpoint(njtree),main=sprintf('NJ tree, lpost=%.2f',p))
  if (sum(njtree$edge.length)>0) axisPhylo(1,backward = F)
} else {
  plot(0,0,main=sprintf('NJ tree not available'),xlab='',ylab='')
}

p=posterior(phangorn::pml(tree=r[[length(r)]],data=seqs,model='JC'),r=phylo2contphylo(r[[length(r)]]),scale=scale)
plot(midpoint(r[[length(r)]]),main=sprintf('Final tree, lpost=%.2f',p))
axisPhylo(1,backward = F)

