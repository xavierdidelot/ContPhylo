rm(list=ls())
library(ContPhylo)
library(ape)
library(phangorn)
set.seed(0)
tree=rtree(6)
scale=0.01
tree$edge.length=tree$edge.length*scale
#the next line moves the root to half-way above first leaf
#tree=unroot(tree);tree=phytools::reroot(tree,1,tree$edge.length[which(tree$edge[,2]==1)]/2)
seqs=phangorn::simSeq(tree,l=10000,rate=1)
seqs=phyDatEmpty(Ntip(tree));names(seqs)<-tree$tip.label
r=mcmc2(seqs,initTree = tree, niter = 1e3,scale=scale)

par(mfrow=c(3,2))
plot(sapply(r,function(x) x$posterior),type='l',xlab='',ylab='Posterior')

p=posterior(phangorn::pml(tree=tree,data=seqs,model='JC'),r=phylo2contphylo(tree),scale=scale)
plot(tree,main=sprintf('Correct tree, lpost=%.2f',p))
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
plot(r[[length(r)]],main=sprintf('Final tree, lpost=%.2f',p))
axisPhylo(1,backward = F)

clusters=rep(1,length(r))
reps=vector('list',length(r))
k=0
for (i in seq_along(r)) {
  assigned=FALSE
  if (k>0) for (j in 1:k)
    if (all.equal(r[[i]], reps[[j]],use.edge.length = F)) {clusters[i]=j;assigned=TRUE;break}
  if (!assigned) {k=k+1;reps[[k]]=r[[i]];clusters[i]=k}
}
ma=max(clusters)
map=sample(1:ma,ma)
plot(map[clusters],xlab='',ylab='')
