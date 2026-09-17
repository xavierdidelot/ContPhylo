#Testing basic package functions
context("Testing basic package functions")

test_that("Testing basic package functions.", {
  for (i in 1:10) {
    set.seed(i)
    p0=rtree(10)
    expect_silent(r<-phylo2contphylo(p0))
    expect_is(r,'contphylo')
    expect_silent(p<-contphylo2phylo(r))
    expect_is(p,'phylo')
    expect_silent(r2<-phylo2contphylo(p))
    expect_is(r2,'contphylo')
    expect_equal(treediff(p0,p),0)
  }
})


