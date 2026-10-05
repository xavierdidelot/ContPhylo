#' Prior function for a contphylo object
#' @param r Object of class contphylo
#' @param scale Scale of the branch lengths
#' @returns Log of the prior
#' @export
prior=function(r,scale=1) {
  tree=contphylo2phylo(r)
  n=nrow(r$coordinates)
  p=priorT(tree,r$type,scale)
  if (r$type==1)  p=p-(n-1)*log(2)
  if (r$type==2) p=p-(n-1)*log(2)+sum(r$coordinates[,2])
  if (r$type==3) p=p-(n-2)*log(2)+sum(r$coordinates[,2])
  return(p)
}

#' Prior function for a phylo object
#' @param r Object of class phylo
#' @param type Type of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted.
#' @param scale Scale of the branch lengths
#' @returns Log of the prior
#' @export
priorT=function(tree,type=2,scale=1) {
  n=Ntip(tree)
  if (type==1) {
    p=(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)#topology
    heights=unname(dist.nodes(tree)[1,(n+1):(2*n-1)])
    heights=sort(heights)
    difs=diff(c(0,heights))
    p=p+sum(dexp(difs,rate=choose(n:2,2)/scale,log=T))#coalescent intervals
  }

  if (type==2) {
    p=(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)#topology
    p=p+sum(dexp(tree$edge.length,rate=1/scale,log=T))
  }

  if (type==3) {
    p=(n-3)*log(2)+lfactorial(n-3)-lfactorial(2*n-5)#topology
    p=p+sum(dexp(tree$edge.length,rate=1/scale,log=T))
  }

  return(p)
}
