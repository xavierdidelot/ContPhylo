#Testing basic package functions
context("Testing convertion package functions")

test_that("Testing convertion package functions.", {
  for (n in 2:10) {
    r=list()
    r$coordinates=rbind(c(0,0),matrix(c(seq(from=1,length.out=n-1,by=4),rep(0,n-1)),nrow=n-1,ncol=2))
    r$tip.label=paste0('t',1:n)
    r$model=2
    class(r)<-'contphylo'
    p<-contphylo2phylo(r)
    mat1=cophenetic.phylo(p)
    mat2=outer(1:n,1:n,\(i, j) ifelse(i == j, 0, abs(i - j) + pmin(i, j, 2)))
    expect_equal(max(abs(as.vector(mat1)-as.vector(mat2))),0)
  }

  for (n in 2:10) {
    r=list()
    r$coordinates=rbind(c(0,0),matrix(c(1/2^(0:(n-2)),rep(0,n-1)),nrow=n-1,ncol=2))
    r$tip.label=paste0('t',1:n)
    r$model=2
    class(r)<-'contphylo'
    p<-contphylo2phylo(r)
    expect_equal(sum(p$edge.length),n)
  }

  for (n in 3:10) {
    r=list()
    r$coordinates=rbind(c(0,0),log(rexp(2)),matrix(0,nrow=n-2,ncol=2))
    s=abs(r$coordinates[2,1])+exp(r$coordinates[2,2])
    for (i in 3:n) {
      r$coordinates[i,1]=runif(1,0,s)
      r$coordinates[i,2]=log(rexp(1))
      s=s+exp(r$coordinates[i,2])
    }
    r$tip.label=paste0('t',1:n)
    r$model=2
    class(r)<-'contphylo'
    p<-contphylo2phylo(r)
    expect_equal(sum(p$edge.length),s)
  }

  #Test properties after decode
  for (n in 3:10) {
    r=rcontphylo(n,1)
    p=contphylo2phylo(r)
    expect_equal(is.ultrametric(p),T)
    expect_equal(is.rooted(p),T)
    r=rcontphylo(n,2)
    p=contphylo2phylo(r)
    expect_equal(is.ultrametric(p),F)
    expect_equal(is.rooted(p),T)
    r=rcontphylo(n,3)
    p=contphylo2phylo(r)
    expect_equal(is.ultrametric(p),F)
    expect_equal(is.rooted(p),F)
  }

  #Test model 2 forward and backward conversions
  for (n in 2:10) {
    r=rtree(n)
    r2=contphylo2phylo(phylo2contphylo(r))
    mat1=cophenetic.phylo(r)
    mat2=cophenetic.phylo(r2)
    expect_equal(max(abs(as.vector(mat1)-as.vector(mat2))),0)
  }

  #Test model 1 forward and backward conversions
  for (n in 2:10) {
    r=rcoal(n)
    r2=contphylo2phylo(phylo2contphylo(r,model=1))
    mat1=cophenetic.phylo(r)
    mat2=cophenetic.phylo(r2)
    expect_equal(max(abs(as.vector(mat1)-as.vector(mat2))),0)
    r3=phangorn::upgma(mat1)
    mat3=cophenetic.phylo(r3)
    expect_equal(max(abs(as.vector(mat1)-as.vector(mat3))),0)
    r4=contphylo2phylo(phylo2contphylo(r3,model=1))
    mat4=cophenetic.phylo(r4)
    expect_equal(max(abs(as.vector(mat1)-as.vector(mat4))),0)
  }

  #Test model 3 forward and backward conversions
  for (n in 3:10) {
    r=rtree(n)
    r=unroot(r)
    r2=contphylo2phylo(phylo2contphylo(r,model=3))
    mat1=cophenetic.phylo(r)
    mat2=cophenetic.phylo(r2)
    expect_equal(max(abs(as.vector(mat1)-as.vector(mat2))),0)
  }

})


