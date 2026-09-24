#' Prior function
#' @param r Object of class ContPhylo
#' @param shape Shape of the gamma distribution for branch lengths
#' @param scale Scale of the gamma distribution for branch lengths
#' @returns Log of the prior
#' @export
prior=function(r,shape=1,scale=1) {
  tree=contphylo2phylo(r)
  n=nrow(r$coordinates)
  p=(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)
  p=p+sum(dgamma(tree$edge.length,shape=shape,scale=scale,log=T))
  p=p+sum(r$coordinates[,2])
  return(p)
}
