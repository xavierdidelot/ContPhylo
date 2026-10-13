rm(list=ls())
library(ape)
library(ContPhylo)

set.seed(0)
n=6
t=rcontphylo(n)
t$tip.label=1:n
as=c(-1,1.5,-1,2.75,2.25)
bs=log(c(1,0.5,2,0.75,0.25))
t$coordinates=matrix(c(0,as,0,bs),nrow=6,ncol=2,byrow=F)
pdf('scripts/figUltra.pdf',10,3)
par(mfrow=c(1,n-1),xpd=NA)
for (i in 2:n) {
  t2=t
  t2$coordinates=t2$coordinates[1:i,]
  t2$tip.label=t2$tip.label[1:i]
  p=contphylo2phylo(t2)
  p=rotateConstr(p,c(4,6,5,2,3,1))
  plot(p,direction='downwards',cex=1.5,x.lim=c(-1,Ntip(p)))
  axisPhylo(2,backward = F)
}
dev.off()
browseURL('scripts/figUltra.pdf')

set.seed(0)
n=6
t=rcontphylo(n)
t$tip.label=1:n
as=c(-1,1.5,-1,2.75,3)
bs=log(c(0.5,1,0.5,1,0.5))
t$coordinates=matrix(c(0,as,0,bs),nrow=6,ncol=2,byrow=F)
pdf('scripts/figRooted.pdf',10,3)
par(mfrow=c(1,n-1),xpd=NA)
for (i in 2:n) {
  t2=t
  t2$coordinates=t2$coordinates[1:i,]
  t2$tip.label=t2$tip.label[1:i]
  p=contphylo2phylo(t2)
  p=rotateConstr(p,c(4,2,6,5,3,1))
  plot(p,direction='downwards',cex=1.5,x.lim=c(-1,Ntip(p)))
  axisPhylo(2,backward = F)
}
dev.off()
browseURL('scripts/figRooted.pdf')

set.seed(0)
n=6
t=rcontphylo(n)
t$tip.label=1:n
as=c(-1,1.5,-1,2.75,0)
bs=log(c(0.5,1,0.5,1,0.5))
t$coordinates=matrix(c(0,as,0,bs),nrow=6,ncol=2,byrow=F)
pdf('scripts/figUnrooted.pdf',10,5)
par(mfrow=c(2,3),xpd=NA)
for (i in 2:n) {
  t2=t
  t2$coordinates=t2$coordinates[1:i,]
  t2$tip.label=t2$tip.label[1:i]
  p=contphylo2phylo(t2)
  p=rotateConstr(p,c(4,2,6,5,3,1))
  plot(p,direction='downwards',cex=1.5)
  axisPhylo(2,backward = F)
}
plot(p,type='unrooted',cex=1.5)
dev.off()
browseURL('scripts/figUnrooted.pdf')
