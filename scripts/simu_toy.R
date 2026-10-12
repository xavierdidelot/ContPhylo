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
dist=dist.ml(seqs)
disttree=phangorn::upgma(dist)
results=vector('list',2)
clusters=vector('list',2)
for (i in 1:2) {
  results[[i]]=mcmc2(seqs,model=model,niter=2e5,scale=scale,initTree = tree)
  clusters[[i]]=topologies(results[[i]])
}
save.image('scripts/simu_toy.RData')

rm(list=ls())
load('scripts/simu_toy.RData')
set.seed(0)




pdf('scripts/simu_toy.pdf',7,7)
r=results[[1]]
#r=r[(length(r)/2+1):length(r)]
par(mfrow=c(3,2),xpd=NA)
p=posterior(phangorn::pml(tree=tree,data=seqs,model='JC'),r=phylo2contphylo(tree,model=model),scale=scale)
plot(tree,main=sprintf('Correct tree, lpost=%.2f',p),cex=1.5)
axisPhylo(1,backward = F)

plot(sapply(r,function(x) x$posterior),type='l',xlab='',ylab='Posterior')

p=posterior(phangorn::pml(tree=r[[length(r)]],data=seqs,model='JC'),r=phylo2contphylo(r[[length(r)]],model=model),scale=scale)
plot(r[[length(r)]],main=sprintf('Final tree, lpost=%.2f',p),cex=1.5)
axisPhylo(1,backward = F)

plot(clusters[[1]],xlab='',ylab='')
counts=table(clusters[[1]])
names(counts)<-1:length(counts)
barplot(counts)
dev.off()
browseURL('scripts/simu_toy.pdf')
