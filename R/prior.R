#' Prior function for a contphylo object
#' @param r Object of class contphylo
#' @param scale Scale of the branch lengths
#' @returns Log of the prior
#' @export
prior=function(r,scale=1) {
  tree=contphylo2phylo(r)
  n=nrow(r$coordinates)
  p=priorT(tree,r$model,scale)
  if (r$model==1)  p=p-(n-1)*log(2)
  if (r$model==2) p=p-(n-1)*log(2)+sum(r$coordinates[,2])
  if (r$model==3) p=p-(n-2)*log(2)+sum(r$coordinates[,2])
  return(p)
}

#' Prior function for a phylo object
#' @param r Object of class phylo
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted.
#' @param scale Scale of the branch lengths
#' @returns Log of the prior
#' @export
priorT=function(tree,model=2,scale=1) {
  n=Ntip(tree)
  if (model==1) {
    heights=unname(dist.nodes(tree)[n+1,c(1,(n+1):(2*n-1))])
    heights=sort(max(heights)-heights)
    difs=diff(heights)
    p=sum(-difs*choose(n:2,2)/scale)#coalescent intervals
  }

  if (model==2) {
    p=(n-1)*log(2)+lfactorial(n-1)-lfactorial(2*n-2)#topology
    p=p+sum(dexp(tree$edge.length,rate=1/scale,log=T))
  }

  if (model==3) {
    p=(n-3)*log(2)+lfactorial(n-3)-lfactorial(2*n-5)#topology
    p=p+sum(dexp(tree$edge.length,rate=1/scale,log=T))
  }

  return(p)
}
