#Testing  prior function
context("Testing prior function")

test_that("Testing prior function.", {
  t=rtree(2)
  r=phylo2contphylo(t)
  expect_equal(prior(r,2),sum(dexp(t$edge.length,rate=1/2,log=T))+sum(r$coordinates[,2]))
  t=rtree(3)
  r=phylo2contphylo(t)
  expect_equal(prior(r,2),log(1/3)+sum(dexp(t$edge.length,rate=1/2,log=T))+sum(r$coordinates[,2]))
})


