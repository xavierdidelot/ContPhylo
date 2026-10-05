# #Testing  stan functions
# context("Testing stan functions")
#
# test_that("Testing stan functions.", {
#
#   #Check stan prior is as expected on a small example
#   t=rtree(2)
#   r=phylo2contphylo(t)
#   expected=sum(dexp(t$edge.length,rate=1/2,log=T))+sum(r$coordinates[,2])
#   expect_equal(prior_stan(r,2),expected,tolerance=1e-6)
#
#   #Check stan prior is as expected on another small example
#   t=rtree(3)
#   r=phylo2contphylo(t)
#   expected=log(1/3)+sum(dexp(t$edge.length,rate=1/2,log=T))+sum(r$coordinates[,2])
#   expect_equal(prior_stan(r,2),expected,tolerance=1e-6)
#
#   #Check stan prior is same as R prior
#   n=10
#   r=rcontphylo(n)
#   p1=prior(r,2)
#   p2=prior_stan(r,2)
#   expect_equal(p1,p2,tolerance=1e-6)
#
#   #Check stan posterior is same as R posterior when there is no data
#   n=10
#   tree=rtree(n)
#   seqs=phyDatEmpty(Ntip(tree));names(seqs)<-tree$tip.label
#   mypml=phangorn::pml(tree=tree,data=seqs,model='JC')
#   r=phylo2contphylo(tree)
#   p1=posterior(mypml,r,2)
#   p2=posterior_stan(mypml,r,2)
#   expect_equal(p1,p2,tolerance=1e-6)
#
#   #Check stan posterior is same as R posterior when there is data
#   n=10
#   tree=rtree(n)
#   scale=0.01
#   tree$edge.length=tree$edge.length*scale
#   seqs=phangorn::simSeq(tree,l=100,rate=1)
#   mypml=phangorn::pml(tree=tree,data=seqs,model='JC')
#   r=phylo2contphylo(tree)
#   p1=posterior(mypml,r,scale)
#   p2=posterior_stan(mypml,r,scale)
#   expect_equal(p1,p2,tolerance=1e-6)
# })
#
#
