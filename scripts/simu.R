rm(list=ls())
library(ContPhylo)
library(ape)
library(phangorn)
set.seed(2)

n=10
model=1
scale=0.01
l=10000

tree=randtree(n,model,scale)
if (l==0) {seqs=phyDatEmpty(Ntip(tree));names(seqs)<-tree$tip.label} else
seqs=phangorn::simSeq(tree,l=l,rate=1)

pdf('scripts/simu.pdf',7,7)
par(mfrow=c(3,2),xpd=NA)
p=posterior(phangorn::pml(tree=tree,data=seqs,model='JC'),r=phylo2contphylo(tree,model=model),scale=scale)
plot(tree,main=sprintf('Correct tree, lpost=%.2f',p),cex=1.5)
axisPhylo(1,backward = F)

dist=dist.ml(seqs)
if ((length(seqs)>2||model==1) && all(!is.na(dist))) {
    if (model==1) disttree=phangorn::upgma(dist)
    if (model==2) disttree=midpoint(NJ(dist))
    if (model==3) disttree=NJ(dist)
  p=NA
  tryCatch(suppressWarnings(
    p<-posterior(phangorn::pml(tree=disttree,data=seqs,model='JC'),r=phylo2contphylo(disttree,model=model),scale=scale)
    ),error = function(e) e)
  plot(disttree,main=sprintf('Distance-based tree, lpost=%.2f',p),cex=1.5)
  if (sum(disttree$edge.length)>0) axisPhylo(1,backward = F)
} else {
  plot(0,0,main=sprintf('Distance-based tree not available'),xlab='',ylab='')
}

r=mcmc(seqs,model=model,niter=2e5,scale=scale,initTree = tree)
#r=r[(length(r)/2+1):length(r)]

plot(sapply(r,function(x) x$posterior),type='l',xlab='',ylab='Posterior')

p=posterior(phangorn::pml(tree=r[[length(r)]],data=seqs,model='JC'),r=phylo2contphylo(r[[length(r)]],model=model),scale=scale)
plot(r[[length(r)]],main=sprintf('Final tree, lpost=%.2f',p),cex=1.5)
axisPhylo(1,backward = F)

clusters=topologies(r)
plot(clusters,xlab='',ylab='')
counts=table(clusters)
names(counts)<-1:length(counts)
barplot(counts)
dev.off()
browseURL('scripts/simu.pdf')
