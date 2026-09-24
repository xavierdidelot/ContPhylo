#' Prior function
#' @param r Object of class ContPhylo
#' @param scale Scale of the exponential distribution for branch lengths
#' @returns Log of the prior
#' @export
prior=function(r,scale=1) {
  tree=contphylo2phylo(r)
  n=nrow(r$coordinates)
  p=(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)
  p=p+sum(dexp(tree$edge.length,rate=1/scale,log=T))
  p=p+sum(r$coordinates[,2])
  return(p)
}
