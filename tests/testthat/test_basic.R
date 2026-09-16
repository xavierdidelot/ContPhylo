#Testing basic package functions
context("Testing basic package functions")

test_that("Testing basic package functions.", {
  set.seed(0)
  expect_silent(r<-rcontphylo(10))
  expect_is(r,'contphylo')
  expect_silent(p<-contphylo2phylo(r))
  expect_is(p,'phylo')
  expect_silent(r2<-phylo2contphylo(p))
  expect_is(r2,'contphylo')
  expect_equal(treediff(p,p),0)
})


