#' Importance sampler
#' @param n Number of leaves
#' @param model Model of tree: 1 for ultrametric, 2 for rooted, 3 for unrooted.
#' @param scale Scale parameter of exponential distribution for branch lengths
#' @param nrep Number of repeated simulations
#' @returns List of objects of class contphylo
#' @export
importance=function(n=5,model=2,scale=1,nrep=1000) {
  res=vector(mode='list',nrep)
  for (rep in 1:nrep) {
    r=rcontphylo(n,model)
    r$w=exp(prior(r,scale=scale)-r$p)
    res[[rep]]=r
  }
  weights=unlist(lapply(res,function(x) {x$w}))
  ess=sum(weights)^2/sum(weights^2)
  print(ess)
  return(res)
}
