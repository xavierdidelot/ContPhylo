#Testing  stan functions
context("Testing stan functions")

test_that("Testing stan functions.", {
  t=rtree(2)
  r=phylo2contphylo(t)
  expected=sum(dexp(t$edge.length,rate=1/2,log=T))+sum(r$coordinates[,2])
  expect_equal(prior_stan(r,2),expected,tolerance=1e-6)

  t=rtree(3)
  r=phylo2contphylo(t)
  expected=log(1/3)+sum(dexp(t$edge.length,rate=1/2,log=T))+sum(r$coordinates[,2])
  expect_equal(prior_stan(r,2),expected,tolerance=1e-6)

  n=10
  r=rcontphylo(n)
  p1=prior(r,2)
  p2=prior_stan(r,2)
  expect_equal(p1,p2,tolerance=1e-6)

})


