#' Posterior function
#' @param mypml Object of class phangorn::pml containing the data
#' @param r Object of class ContPhylo
#' @param scale Scale of branch lengths
#' @returns Log of the posterior
#' @export
posterior=function(mypml,r,scale=1) {
  tree=contphylo2phylo(r)
  n=nrow(r$coordinates)
  p=priorT(tree,model=r$model,scale)
  if (r$model==1)  p=p-(n-1)*log(2)
  if (r$model==2) p=p-(n-1)*log(2)+sum(r$coordinates[,2])
  if (r$model==3) p=p-(n-2)*log(2)+sum(r$coordinates[,2])
  if (sum(attr(mypml$data, "weight"))>0) {
    mypml=stats::update(mypml,tree=tree)
    p=p+mypml$logLik
  }
  return(p)
}
