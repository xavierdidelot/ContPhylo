rm(list=ls())
library(ContPhylo)
library(ape)
library(phangorn)
set.seed(0)

n=5
model=1
scale=0.01
l=1000

tree=randtree(n,model,scale)
if (l==0) {seqs=phyDatEmpty(Ntip(tree));names(seqs)<-tree$tip.label} else
  seqs=phangorn::simSeq(tree,l=l,rate=1)

n=length(seqs)
initTree=randtree(n,model,0.01)
initTree$tip.label=names(seqs)

r=mle(seqs,model=model,initTree=initTree)
par(mfrow=c(1,3))
plot(tree,cex=1.5)
mypml=phangorn::pml(tree=tree,data=seqs,model='JC')
title(sprintf('Truth ll=%e',mypml$logLik))
axisPhylo(1,backward = F)
plot(initTree,cex=1.5)
mypml=phangorn::pml(tree=initTree,data=seqs,model='JC')
title(sprintf('Init ll=%e',mypml$logLik))
axisPhylo(1,backward = F)
plot(r,cex=1.5)
title(sprintf('MLE ll=%e',r$ll))
axisPhylo(1,backward = F)
