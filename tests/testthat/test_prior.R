#Testing  prior function
context("Testing prior function")

test_that("Testing prior function.", {
  t=rtree(2)
  r=phylo2contphylo(t)
  expect_equal(prior(r,2,2),sum(dgamma(t$edge.length,shape=2,scale=2,log=T))+sum(r$coordinates[,2]))
  t=rtree(3)
  r=phylo2contphylo(t)
  expect_equal(prior(r,2,2),log(1/3)+sum(dgamma(t$edge.length,shape=2,scale=2,log=T))+sum(r$coordinates[,2]))
})


