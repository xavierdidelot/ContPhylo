#Testing basic package functions
context("Testing sum of branch lengths equation")

test_that("Testing sum of branch lengths equation.", {
  for (i in 1:10) for (model in 1:3) {
    set.seed(i)
    r=rcontphylo(10,model=model)
    s=sum_bra_len(r)
    p=contphylo2phylo(r)
    expect_equal(s,sum(p$edge.length))
  }
})


