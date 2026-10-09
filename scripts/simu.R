rm(list=ls())
library(ContPhylo)
library(ape)
library(phangorn)
set.seed(1)

n=6
model=1
scale=0.01
l=1000

tree=randtree(n,model,scale)
if (l==0) {seqs=phyDatEmpty(Ntip(tree));names(seqs)<-tree$tip.label} else
seqs=phangorn::simSeq(tree,l=l,rate=1)
r=mcmc2(seqs,model=model,niter=1e4,scale=scale)

par(mfrow=c(3,2))
plot(sapply(r,function(x) x$posterior),type='l',xlab='',ylab='Posterior')

p=posterior(phangorn::pml(tree=tree,data=seqs,model='JC'),r=phylo2contphylo(tree,model=model),scale=scale)
plot(tree,main=sprintf('Correct tree, lpost=%.2f',p))
axisPhylo(1,backward = F)

dist=dist.ml(seqs)
if ((length(seqs)>2||model==1) && all(!is.na(dist))) {
    if (model==1) disttree=phangorn::upgma(dist)
    if (model==2) disttree=midpoint(NJ(dist))
    if (model==3) disttree=NJ(dist)
  p=NA
  tryCatch(suppressWarnings(
    p<-posterior(phangorn::pml(tree=disttree,data=seqs,model='JC'),r=phylo2contphylo(disttree,model=model),scale=scale)),
    error = function(e) e)
  plot(disttree,main=sprintf('Distance-based tree, lpost=%.2f',p))
  if (sum(disttree$edge.length)>0) axisPhylo(1,backward = F)
} else {
  plot(0,0,main=sprintf('Distance-based tree not available'),xlab='',ylab='')
}

p=posterior(phangorn::pml(tree=r[[length(r)]],data=seqs,model='JC'),r=phylo2contphylo(r[[length(r)]],model=model),scale=scale)
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
counts=table(clusters)
b=counts[map]
names(b)<-1:length(b)
barplot(b)
