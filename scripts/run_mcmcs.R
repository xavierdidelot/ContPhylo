rm(list=ls())
library(ape)
library(ContPhylo)
source('scripts/fig_label.R')

set.seed(0)

n=5#number of leaves in tree
scale=0.01#scaling factor on branch lengths

alltrees=vector('list',3)
allcounts=vector('list',3)
for (model in 1:3) {
  trees=mcmc(phyDatEmpty(n),model=model,scale=scale,niter=2e5,thin=10)
  trees=trees[(length(trees)/2+1):length(trees)]
  alltrees[[model]]=trees

  clusters=topologies(trees,ranked=(model==1))
  allcounts[[model]]=table(clusters)
}
save.image('scripts/run_mcmcs.RData')

rm(list=ls())
load('scripts/run_mcmcs.RData')
set.seed(0)
pdf('scripts/priors.pdf',7,7)
let=1
par(mfrow=c(3,3))
for (model in 1:3) {
  trees=alltrees[[model]]
  counts=allcounts[[model]]
  ma=length(counts)
  map=sample(1:ma,ma)
  b=counts[map]
  names(b)<-1:length(b)
  barplot(b)
  lines(c(0,2*length(b)),c(1,1)*sum(b)/length(b),col='red')
  fig_label(LETTERS[let],cex=2);let=let+1

  b1=sapply(trees,function(x) {wi=which(x$tip.label==trees[[1]]$tip.label[1]);w=which(x$edge[,2]==wi);return(x$edge.length[w])})
  #b1=as.vector(sapply(trees,function(x) {w=which(x$edge[,2]<=n);return(x$edge.length[w])}))
  hist(b1[b1<=0.05],main='',xlab='',probability=T,breaks=seq(0,0.05,length.out=20))
  xs=seq(0,0.05,length.out=1000)
  if (model==1) lines(xs,sapply(xs,function (t) ContPhylo:::p_t(t,n,scale)),col='red')
  if (model>1)  lines(xs,dexp(xs,rate=1/scale),col='red')
  fig_label(LETTERS[let],cex=2);let=let+1

  sbl=sapply(trees,function(x) return(sum(x$edge.length)))
  hist(sbl[sbl<=0.2],main='',xlab='',probability=T,breaks=seq(0,0.2,length.out=20))
  xs=seq(0,0.2,length.out=1000)
  dcoal=function(x) {(n-1)/(2*scale)*exp(-x/(2*scale))*(1-exp(-x/(2*scale)))^(n-2)}
  par(xpd=NA)
  if (model==1) lines(xs,dcoal(xs),col='red')
  if (model==2) lines(xs,dgamma(xs,shape=(2*n-2),scale=scale),col='red')
  if (model==3) lines(xs,dgamma(xs,shape=(2*n-3),scale=scale),col='red')
  par(xpd=F)
  fig_label(LETTERS[let],cex=2);let=let+1
}
dev.off()
browseURL('scripts/priors.pdf')
