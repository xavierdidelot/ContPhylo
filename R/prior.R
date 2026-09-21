#' Prior function
#' @param tree Tree in ape::phylo format
#' @param shape Shape of the gamma distribution for branch lengths
#' @param scale Scale of the gamma distribution for branch lengths
#' @returns Log of the prior
#' @export
prior=function(tree,shape=1,scale=1) {
  n=Ntip(tree)
  p=(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)
  p=p+sum(dgamma(tree$edge.length,shape=shape,scale=scale,log=T))
  return(p)
}
