#Testing basic package functions
context("Testing convertion package functions")

test_that("Testing convertion package functions.", {
  for (n in 2:10) {
  r=list()
  r$coordinates=rbind(c(0,0),matrix(c(seq(from=1,length.out=n-1,by=4),rep(1,n-1)),nrow=n-1,ncol=2))
  r$tip.label=paste0('t',1:n)
  class(r)<-'contphylo'
  r$coordinates=log(r$coordinates)
  p<-contphylo2phylo(r)
  mat1=cophenetic.phylo(p)
  mat2=outer(1:n,1:n,\(i, j) ifelse(i == j, 0, abs(i - j) + pmin(i, j, 2)))
  expect_equal(max(abs(as.vector(mat1)-as.vector(mat2))),0)
  }

  for (n in 2:10) {
    r=list()
    r$coordinates=rbind(c(0,0),matrix(c(1/2^(0:(n-2)),rep(1,n-1)),nrow=n-1,ncol=2))
    r$tip.label=paste0('t',1:n)
    class(r)<-'contphylo'
    r$coordinates=log(r$coordinates)
    p<-contphylo2phylo(r)
    expect_equal(sum(p$edge.length),n)
  }

  for (n in 3:10) {
    r=list()
    r$coordinates=rbind(c(0,0),rexp(2),matrix(0,nrow=n-2,ncol=2))
    s=r$coordinates[2,1]+r$coordinates[2,2]
    for (i in 3:n) {
      r$coordinates[i,1]=runif(1,0,s)
      r$coordinates[i,2]=rexp(1)
      s=s+r$coordinates[i,2]
    }
    r$tip.label=paste0('t',1:n)
    class(r)<-'contphylo'
    r$coordinates=log(r$coordinates)
    p<-contphylo2phylo(r)
    expect_equal(sum(p$edge.length),s)
  }
})


