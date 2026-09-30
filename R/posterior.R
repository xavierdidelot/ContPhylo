#' Posterior function
#' @param mypml Object of class phangorn::pml containing the data
#' @param r Object of class ContPhylo
#' @param scale Scale of the exponential distribution for branch lengths
#' @returns Log of the posterior
#' @export
posterior=function(mypml,r,scale=1) {
  tree=contphylo2phylo(r)
  n=nrow(r$coordinates)
  p=(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)
  p=p+sum(dexp(tree$edge.length,rate=1/scale,log=T))
  p=p+sum(r$coordinates[,2])
  mypml=stats::update(mypml,tree=tree)
  p=p+mypml$logLik
  #p=p+dnorm(dist.nodes(tree)[1,n+1],scale/2,0.1,log=T)#Encourages root position
  return(p)
}
