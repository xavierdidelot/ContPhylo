rm(list=ls())
library(ContPhylo)
library(ape)
library(phangorn)

n=10
model=1
scale=0.004
l=10000

set.seed(4)
tree=randtree(n,model,scale)#nice balanced tree

set.seed(2)
if (l==0) {seqs=phyDatEmpty(Ntip(tree));names(seqs)<-tree$tip.label} else
seqs=phangorn::simSeq(tree,l=l,rate=1)
dist=dist.ml(seqs)
disttree=phangorn::upgma(dist)
results=vector('list',2)
#clusters=vector('list',2)
for (i in 1:2) {
  x=mcmc2(seqs,model=model,niter=1e5,scale=scale,initTree = disttree)
  x=x[(floor(length(x) / 2) + 1):length(x)]
  results[[i]]=x
  #clusters[[i]]=topologies(x)
}
clusters=topologies(c(results[[1]],results[[2]]),ranked = F)
save.image('scripts/simu_toy.RData')

rm(list=ls())
load('scripts/simu_toy.RData')
source('scripts/topo_ess.R')
set.seed(0)

mcmcs=vector('list',2)
for (i in 1:2) {
  trees=results[[i]]
  data=sapply(trees,function(x) return(sum(x$edge.length)))
  for (wi in 1:n) {
    data=cbind(data,sapply(trees,function(x) {w=which(x$edge[,2]==wi);return(x$edge.length[w])}))
  }
  mcmcs[[i]]=coda::as.mcmc(data)
  ess=coda::effectiveSize(mcmcs[[i]])
  print(sprintf('ESS values: %s',paste(round(ess),sep=',',collapse = ',')))
  #tess=topo_ess(trees)
  #print(sprintf('Approx ESS for topology: %f',tess))
}
coda::gelman.diag(coda::as.mcmc.list(mcmcs),multivariate = F)

source('scripts/fig_label.R')
let=1
pdf('scripts/simu_toy.pdf',10,7)
par(mfrow=c(2,3),xpd=NA)
plot(tree,cex=1.5)
axisPhylo(1,backward = F)
fig_label(LETTERS[let],cex=2);let=let+1

#plot(disttree,cex=1.5)
#axisPhylo(1,backward = F)

r=results[[1]]
r=r[(length(r)/2+1):length(r)]
class(r)<-'multiPhylo'
densiTree(r,type='cladogram',col='blue',alpha = 0.01)
fig_label(LETTERS[let],cex=2);let=let+1

r=results[[2]]
r=r[(length(r)/2+1):length(r)]
class(r)<-'multiPhylo'
densiTree(r,type='cladogram',col='blue',alpha = 0.01)
fig_label(LETTERS[let],cex=2);let=let+1

col1 <- rgb(1, 0, 0, 0.4)   # red, 40% opaque
col2 <- rgb(0, 0, 1, 0.4)   # blue, 40% opaque
s1=sapply(results[[1]],function(x) return(sum(x$edge.length)))
s2=sapply(results[[2]],function(x) return(sum(x$edge.length)))
plot (s1,type='l',xlab='',ylab='',col=col1,ylim=c(min(c(s1,s2)),max(c(s1,s2))))
lines(s2,col=col2)
fig_label(LETTERS[let],cex=2);let=let+1

counts=table(clusters[1:1000])
#names(counts)<-1:length(counts)
barplot(counts)
fig_label(LETTERS[let],cex=2);let=let+1

counts=table(clusters[1001:2000])
#names(counts)<-1:length(counts)
barplot(counts)
fig_label(LETTERS[let],cex=2);let=let+1

dev.off()
browseURL('scripts/simu_toy.pdf')
