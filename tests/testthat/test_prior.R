#Testing  prior function
context("Testing prior function")

test_that("Testing prior function.", {
  t=rtree(2)
  expect_equal(prior(t,2,2),sum(dgamma(t$edge.length,shape=2,scale=2,log=T)))
  t=rtree(3)
  expect_equal(prior(t,2,2),log(1/3)+sum(dgamma(t$edge.length,shape=2,scale=2,log=T)))
})


