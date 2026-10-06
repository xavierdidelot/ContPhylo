#Testing  prior function
context("Testing prior function")

test_that("Testing prior function.", {
  t=rtree(2)
  r=phylo2contphylo(t)
  expect_equal(prior(r,2),sum(dexp(t$edge.length,rate=1/2,log=T))+sum(r$coordinates[,2])-log(2))
  t=rtree(3)
  r=phylo2contphylo(t)
  expect_equal(prior(r,2),log(1/3)+sum(dexp(t$edge.length,rate=1/2,log=T))+sum(r$coordinates[,2])-2*log(2))

  t=rtree(10)
  r=phylo2contphylo(t)
  data=phyDatEmpty(10)
  names(data)=t$tip.label
  mypml=phangorn::pml(tree=t,data,model='JC')
  expect_equal(prior(r,scale=2),posterior(mypml,r,scale=2))

  t=rcoal(3,br=c(1,2))
  a=priorT(t,model=1,scale=1)
  b=dexp(1,3,T)+dexp(2,1,T)+log(1/3)
  expect_equal(a,b)
})


